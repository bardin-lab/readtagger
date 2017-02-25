from itertools import chain


def namedtuple_to_argv(nt, prog='prog.py'):
    """Convert a namedtuple into an argv list for testing argparser arguments."""
    d = nt._asdict()
    keys = ["--%s" % k for k in d]
    argv = list(chain.from_iterable(zip(keys, d.values())))
    argv = [arg for arg in argv if not isinstance(arg, bool)]
    argv = [arg if not isinstance(arg, list) else " ".join(arg) for arg in argv]
    argv.insert(0, prog)
    return argv
