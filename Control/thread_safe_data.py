""" ThreadDataContainer
Is subclassed by State\n
Uses a Lock objects to make thread-safe.\n
Note: Can only protect Python Objects (most data types in Python)

Module Methods
  - thread_protected() : Provided a Python Object it creates an anonymous queue.Queue
                         object and returns a ThreadDataContainer for thread-safe access to the Object

Class Fields
  - _memory : Stores a queue.Queue object whose sole member is the protected data
  - _lock   : Used to temporarily store the protected data during enter/exit use

Class Methods
  - __contains__() : Allows use of 'in' keyword for top-level of the protected object
  - __getitem__()  : Raises NotImplementedError
                     Python's get-[] method does not allow safe access when multiple []'s are involved
  - __setitem__()  : Raises NotImplementedError
                     Python's set-[] method does not allow safe access when multiple []'s are involved
  - callattr()     : Allows use of top-level items' attributes as functions (better to use With-Statement)
  - __enter__()    : Implemented to allow With-statement use
  - __exit__()     : Implemented to allow With-statement use

@author: Ben C
"""

import sys
import traceback
from threading import Lock


class ThreadSafeDataContainer:
    """ Safety wrapper around a Python Object

    Locks object when accessing an always* releases lock when done (* done via finally statements)

    Use within a WITH statement (e.g. 'with ThreadSafeDataContainer as data: ...access data according to its data type')
    or use a built-in method.  Does not support direct [] access/mutation
    """
    def __init__(self, data):
        """ Makes a thread-safe data container

        :param data: The protected object
        """
        self._memory = data
        self._lock = Lock()

    def __contains__(self, item):
        """ Allows use of 'in' keyword on protected object from ThreadDataContainer level

        Top-level only

        :param item: item being searched for
        :return: True/False if present at top level
        """
        self._lock.acquire()
        try:
            return item in self._memory
        finally:
            self._lock.release()

    def __getitem__(self, key):
        """ Forbidden, [] access not thread-safe due to Python internal implementation

        Technically could be allowed but since [] absolutely can't, to prevent confusion
        with [] for get but set() for set, both are blocked
        """
        raise NotImplementedError(f"For thread-safe access, use get() or call from WITH statement")

    def __setitem__(self, key, value):
        """ Forbidden, [] mutation not thread-safe due to Python internal implementation """
        raise NotImplementedError(f"For thread-safe mutation, use set() or call from WITH statement")

    def __dict__(self):
        return {"ThreadSafeDataContainer": self._memory}

    def callattr(self, attr, *args, **kwargs):
        """ Allows the protected object to call one of its attributes from the ThreadDataContainer level

        :param attr: Attribute name
        :param args: Positional arguments for the attribute
        :param kwargs: Keyword arguements for the attribute
        :return: The return value of the attribute being called
        """
        self._lock.acquire()
        try:
            return self._memory.__getattribute__(attr)(*args, **kwargs)
        finally:
            self._lock.release()

    def __enter__(self):
        """ Implemented to allow WITH statement use

        :return: The protected object
        """
        self._lock.acquire()
        return self._memory

    def __exit__(self, exc_type, exc_value, tb):
        """ Implemented to allow WITH statement use

        :param exc_type: Exception type
        :param exc_value: Exception value
        :param tb: Traceback object
        :return: None
        """
        if exc_type is not None:
            print(f"ThreadDataContainer:: \n\t- type = {exc_type}\n\t- value = {exc_value}\n\t- traceback = {tb}")
            traceback.print_tb(tb, file=sys.stdout)
        if self._memory is None:
            raise ValueError(f"If self.memory is None from an __exit__, caller destroyed object?")
        self._lock.release()
