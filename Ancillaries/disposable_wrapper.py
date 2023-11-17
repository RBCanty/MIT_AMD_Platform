"""
Wrapper function for deploying API-based classes of dubious context management

Intended Usage::

| with disposable(disposableObject):
|   blah...
| # or
| with disposable(disposableObject) as temp:
|   blah...
"""

from inspect import isclass, isfunction
from contextlib import contextmanager


@contextmanager
def disposable(obj_or_class, *args, **kwargs):
    """
    Creates a context manager which implements a 'Dispose' method

    Uses: For Spark and HPLC's IDisposables

    Note about contextmanager and WITH: If this method exits via RETURN the contents of the WITH block are not
    executed, it is only if YIELD is invoked that the WITH block will be executed.

    :param obj_or_class: An object or class being managed
    :param args: The positional arguments for the object's constructor or the method
    :param kwargs: The keyword arguments for the object's constructor or the method
    :return:
    """
    # If it takes arguments, pass them in; otherwise, don't
    if isclass(obj_or_class) or isfunction(obj_or_class):
        obj = obj_or_class(*args, **kwargs)
    else:
        obj = obj_or_class

    # Is it already context-managed?
    if hasattr(obj, '__enter__') and hasattr(obj, '__exit__'):
        return obj
    # If obj does not have a Dispose, then it should not be managed by Disposable
    elif (not hasattr(obj, 'Dispose')) or (not callable(obj.Dispose)):
        return obj
    # Yield for @contextmanager
    try:
        yield obj
    finally:
        obj.Dispose()
