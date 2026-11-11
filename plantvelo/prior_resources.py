"""Shared validation and directory locks for managed prior resources.

Readers and writers lock the existing directory inode, requiring no writable
lock file. Callers hold this lock across their complete resource transaction.
"""
from contextlib import contextmanager
import fcntl
import json
import os
from pathlib import Path
import re
import stat

PENDING_NAME = '.prior-sync.pending.json'
REGISTRY_NAME = 'priors.json'
_SHA256 = re.compile(r'^[0-9a-f]{64}$')
_SPECIES = re.compile(r'^[a-z0-9][a-z0-9_-]*$')
_FIELDS = {'reference_id', 'annotation_id', 'reference_sha256', 'release_id', 'sha256', 'rows'}


@contextmanager
def resource_lock(directory, exclusive=False):
    """Lock the directory itself; never create or modify resource files."""
    fd = os.open(str(directory), os.O_RDONLY | os.O_DIRECTORY)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX if exclusive else fcntl.LOCK_SH)
        yield
    finally:
        fcntl.flock(fd, fcntl.LOCK_UN)
        os.close(fd)


def _identifier(value):
    return (isinstance(value, str) and bool(value) and value not in ('.', '..')
            and not any(c in value for c in ('/', '\\', '\x00', '\n', '\r')))


def validate_entry(species, entry):
    """Validate one exact registry entry, raising ValueError on invalid data."""
    if not isinstance(species, str) or _SPECIES.fullmatch(species) is None:
        raise ValueError('invalid prior species: {}'.format(species))
    if not isinstance(entry, dict) or set(entry) != _FIELDS:
        raise ValueError('invalid prior registry entry fields for {}'.format(species))
    for name in ('reference_id', 'annotation_id', 'release_id'):
        if not _identifier(entry[name]):
            raise ValueError('invalid {} for {}'.format(name, species))
    hashes = entry['reference_sha256']
    if not isinstance(hashes, dict) or set(hashes) != {'genome', 'annotation'}:
        raise ValueError('invalid reference_sha256 for {}'.format(species))
    for value in [entry['sha256']] + list(hashes.values()):
        if not isinstance(value, str) or _SHA256.fullmatch(value) is None:
            raise ValueError('invalid sha256 for {}'.format(species))
    if type(entry['rows']) is not int or entry['rows'] <= 0:
        raise ValueError('invalid rows for {}'.format(species))


def regular_resource(path, required=True):
    """Reject symlinks (including dangling ones) and special resource files."""
    path = Path(path)
    try:
        mode = path.lstat().st_mode
    except FileNotFoundError:
        if not required:
            return False
        raise ValueError('missing prior resource: {}'.format(path))
    if stat.S_ISLNK(mode):
        raise ValueError('symlink prior resource prohibited: {}'.format(path))
    if not stat.S_ISREG(mode):
        raise ValueError('prior resource must be a regular file: {}'.format(path))
    return True


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError('duplicate JSON key in prior registry: {}'.format(key))
        result[key] = value
    return result


def read_registry(directory, required=False):
    """Read and validate registry only; callers manage locking and recovery."""
    path = Path(directory) / REGISTRY_NAME
    if not regular_resource(path, required):
        return None
    try:
        with path.open(encoding='utf-8') as handle:
            registry = json.load(handle, object_pairs_hook=_unique_object)
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError('invalid prior registry {}: {}'.format(path, exc)) from exc
    if (not isinstance(registry, dict) or set(registry) != {'schema_version', 'species'}
            or type(registry['schema_version']) is not int or registry['schema_version'] != 1
            or not isinstance(registry['species'], dict)):
        raise ValueError('invalid prior registry schema: {}'.format(path))
    for species, entry in registry['species'].items():
        validate_entry(species, entry)
    return registry


def load_managed_ir_prior(path, species, required=False, source_kind='builtin'):
    """Load a verified snapshot while holding the shared resource lock."""
    from dataclasses import replace
    from types import MappingProxyType
    from plantvelo.ir_prior import load_builtin_ir_prior, sha256_file
    directory = Path(path).parent
    with resource_lock(directory):
        pending = directory / PENDING_NAME
        if regular_resource(pending, required=False):
            raise ValueError('prior sync is pending; recover the resource sync before loading')
        registry = read_registry(directory, required=required)
        regular_resource(path)
        if registry is None:
            return replace(load_builtin_ir_prior(path, species), source_kind=source_kind)
        if species not in registry['species']:
            raise ValueError('species {} missing from prior registry'.format(species))
        entry = registry['species'][species]
        digest = sha256_file(path)
        if digest != entry['sha256']:
            raise ValueError('prior sha256 mismatch for {}'.format(species))
        prior = load_builtin_ir_prior(path, species)
        if len(prior.fieldnames) != 5:
            raise ValueError('managed prior must contain exactly five columns')
        if (prior.sha256 != digest or prior.input_rows != entry['rows']
                or len(prior.records) != entry['rows']
                or len({record.key for record in prior.records}) != entry['rows']):
            raise ValueError('prior rows/unique records/hash mismatch for {}'.format(species))
        return replace(prior, source_kind=source_kind,
                       reference_id=entry['reference_id'], annotation_id=entry['annotation_id'],
                       release_id=entry['release_id'],
                       reference_sha256=MappingProxyType(dict(entry['reference_sha256'])))
