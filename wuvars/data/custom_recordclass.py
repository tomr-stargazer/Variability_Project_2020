"""
Provides a variant of `recordclass` that is indexable by strings.
"""

from recordclass import recordclass as _recordclass


def recordclass(name, fields, **kwargs):
    """
    A wrapper around recordclass() that restores dict-like string indexing.
    """
    Base = _recordclass(name, fields, **kwargs)

    class Wrapped:
        __slots__ = ("_base",)

        def __init__(self, *args, **kw):
            self._base = Base(*args, **kw)

        def __getitem__(self, key):
            if isinstance(key, str):
                return getattr(self._base, key)
            return self._base[key]

        def __setitem__(self, key, value):
            if isinstance(key, str):
                return setattr(self._base, key, value)
            self._base[key] = value

        def __getattr__(self, name):
            # forward attributes to the underlying recordclass
            return getattr(self._base, name)

        def __setattr__(self, name, value):
            if name == "_base":
                object.__setattr__(self, name, value)
            else:
                setattr(self._base, name, value)

        def __repr__(self):
            return repr(self._base)

    Wrapped.__name__ = name
    return Wrapped
