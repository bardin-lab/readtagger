import functools
import traceback


def reraise_with_stack(func):
    """Include original stacktrace in multiprocessing context."""
    # Adapted from https://stackoverflow.com/a/29357032
    @functools.wraps(func)
    def wrapped(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            traceback_str = traceback.format_exc(e)
            raise Exception("Error occurred. Original traceback is\n%s\n" % traceback_str)

    return wrapped
