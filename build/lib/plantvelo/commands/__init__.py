"""Lazy public command exports."""

__all__ = ["cli", "run", "merge"]


def __getattr__(name):
    if name == "cli":
        from plantvelo.commands.plantvelo import cli

        return cli
    if name == "run":
        from plantvelo.commands.run import run

        return run
    if name == "merge":
        from plantvelo.commands.merge import merge

        return merge
    raise AttributeError(name)
