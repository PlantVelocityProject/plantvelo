"""Lazy public command exports."""
from importlib import import_module

__all__ = ["cli", "run", "merge", "irprior"]


def __getattr__(name):
    if name in __all__:
        module = import_module("plantvelo.commands." + ("plantvelo" if name == "cli" else name))
        command = getattr(module, name)
        # Importing a same-named submodule binds a module on this package.
        # Restore the documented public command object before from-import
        # checks the package attribute a second time.
        globals()[name] = command
        return command
    raise AttributeError(name)
