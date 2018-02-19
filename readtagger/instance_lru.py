from functools import wraps
try:
    from functools import lru_cache
except ImportError:
    from backports.functools_lru_cache import lru_cache

lru_cache = lru_cache


def instance_method_lru_cache(*cache_args, **cache_kwargs):
    """
    Wrap instance methods with a LRU cache.

    Just like functools.lru_cache, but a new cache is created for each instance
    of the class that owns the method this is applied to.
    """
    # from https://gist.github.com/z0u/9df24dda2b1fe0613a85e7349d5f7d62
    def cache_decorator(func):
        @wraps(func)
        def cache_factory(self, *args, **kwargs):
            # Wrap the function in a cache by calling the decorator
            instance_cache = lru_cache(*cache_args, **cache_kwargs)(func)
            # Bind the decorated function to the instance to make it a method
            instance_cache = instance_cache.__get__(self, self.__class__)
            setattr(self, func.__name__, instance_cache)
            # Call the instance cache now. Next time the method is called, the
            # call will go directly to the instance cache and not via this
            # decorator.
            return instance_cache(*args, **kwargs)
        return cache_factory
    return cache_decorator
