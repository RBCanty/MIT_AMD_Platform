""" Operation

Proxy of a Python Dictionary

Stores requisite information for a job on the platform

Class Fields
  - FUNC   : A keyword for a function (looked up on command list)
  - ARGS   : A tuple or positional arguments
  - KWARGS : A dictionary of keyword arguments
  - AGENT  : The name of the system which should perform the operation
  - RETVAL : Completed (True), Encountered a problem or error (False), Not yet run (None)

Class Methods
  - __dict__()     : (Implemented to bypass json not knowing how to serialize nested dictionary-proxies)
  - __deepcopy__() : (Implemented to bypass json not knowing how to serialize nested dictionary-proxies)
  - keys()         : Returns a dictionary list of keys
  - prints()       : Returns a string containing the agent of the function and the keyword of the function
  - package()      : Returns a json string representation of the
  - __setitem__()  : Allows Python-Dictionary-like mutation of fields
  - __getitem__()  : Allows for []-keyed access of fields

Static Methods
  - build_from_json()         : Helper function for constructing an Operation from a dictionary made by json
  - build_from_args()         : Helper function for constructing an Operation
  - generate_item()           : Helper function for constructing an Operation
  - generate_item_from_list() : Helper function for constructing an Operation
  - generate_item_from_dict() : Helper function for constructing an Operation

@author: Ben C
"""

import json
import copy
from typing import Union

FUNC = 'function'
ARGS = 'args'
KWARGS = 'kwargs'
AGENT = 'agent'
RETVAL = 'return value'
OPERATION_KEYS = [FUNC, ARGS, KWARGS, AGENT, RETVAL]
PID = 'plate_id'


class Operation(object):
    """ A single mid-level command

    AGENT - The system which shall be performing the function
    FUNC - A function name or command string
    ARGS - Positional arguments (args[0] = PlateID)
    KWARGS - Keyword Arguments
    RETVAL - (Book-keeping) stores the return code from the function
    """
    def __init__(self, dictionary):
        """ Creates a new operation from a dictionary

        Gives missing values default values, does not remove nonstandard values

        :param dictionary: A dictionary
        """
        self._dictionary = dict()
        for k in dictionary.keys():
            try:
                self[k] = dictionary[k]
            except KeyError:
                pass
        for k in OPERATION_KEYS:
            self[k] = self[k] if k in self.keys() else None

    def __dict__(self):
        """ Provides a dictionary

        Implemented for json

        :return: The underlying dictionary object
        """
        return self._dictionary

    def __deepcopy__(self, memodict={}):  # noqa # This is the recommended syntax for deploying deepcopy
        """ Provides a dictionary

        Implemented for json

        :param memodict: Needed by deepcopy
        :return: A deep copy of the underlying dictionary
        """
        return copy.deepcopy(self._dictionary)

    def keys(self):
        """ Key access

        :return: The keys of the underlying dictionary
        """
        return self._dictionary.keys()

    def prints(self):
        """ A print-string method

        :return: A summary of the operation "agent>function"
        """
        return f"'{self[AGENT]}'>{self[FUNC]}"

    def package(self):
        """ Converts the Operation into a string

        :return: A json-generated string
        """
        return json.dumps(self._dictionary)

    def __setitem__(self, key, item):
        """ Provides [] mutation

        :param key: Dictionary keyword
        :param item: Value overwriting the contents of [key]
        :return: None
        """
        if key not in OPERATION_KEYS:
            raise KeyError(f"Operation cannot have field {key}")
        if item is None:
            if key == FUNC:
                self._dictionary[key] = str()
            if key == ARGS:
                self._dictionary[key] = tuple()
            if key == KWARGS:
                self._dictionary[key] = dict()
            if key == AGENT:
                self._dictionary[key] = str()
            if key == RETVAL:
                self._dictionary[key] = None
        else:
            if (key == FUNC) and (not isinstance(item, str)):
                raise ValueError(f"Operation function must be a string")
            if (key == ARGS) and isinstance(item, list):
                item = tuple(item)
            if (key == ARGS) and (not isinstance(item, tuple)):

                raise ValueError(f"Operation positional arguments must be a tuple")
            if (key == KWARGS) and (not isinstance(item, dict)):
                raise ValueError(f"Operation keyword arguments must be a dictionary")
            if (key == AGENT) and not isinstance(item, str):
                raise ValueError(f"Operation metadata must be a string")
            if (key == RETVAL) and (not isinstance(item, bool)):
                raise ValueError(f"Operation's return value must be a boolean")
            self._dictionary[key] = item

    def __getitem__(self, key):
        """ Allows [] access

        :param key: dictionary key of accessed item
        :return: The item (or raises KeyError if not present)
        """
        return self._dictionary[key]

    @staticmethod
    def build_from_json(dictionary):
        """ Builds an object from a json-made dictionary

        Originally different from the class constructor because json has no list/tuple distinction
        But with a more detailed __setitem__() method this became redundant, kept for semantics

        :param dictionary: json-made dictionary of an Operation
        :return: An Operation object
        """
        init = dict()
        init[AGENT] = dictionary.get(AGENT, None)
        init[FUNC] = dictionary.get(FUNC, None)
        init[ARGS] = dictionary.get(ARGS, None)
        init[KWARGS] = dictionary.get(KWARGS, None)
        init[RETVAL] = dictionary.get(RETVAL, None)
        return Operation(init)

    @staticmethod
    def build_from_args(func, args=None, kwargs=None, agent=None, retval=None):
        """ Allows an operation to be built from argumetns rather than a dictionary

        :param func: String containing the function name
        :param args: Positional arguments in a tuple
        :param kwargs: Keyword arguments in a dictionary
        :param agent: Name of the system which will perform func
        :param retval: Book-keeping the return value of the function
        :return: Operation object
        """
        self = dict()
        self[FUNC] = func
        self[ARGS] = args
        self[KWARGS] = kwargs
        self[AGENT] = agent
        self[RETVAL] = retval
        return Operation(self)

    @staticmethod
    def generate_item(_agent, _function, _args, _kwargs):
        """ Generates operation from arguments but less restrictive

        :param _agent: Name of the system which will perform func
        :param _function: String containing the function name
        :param _args: Positional arguments in a tuple
        :param _kwargs: Keyword arguments in a dictionary
        :return: Operation object
        """
        item = dict()
        item[AGENT] = _agent if isinstance(_agent, str) else str()
        item[FUNC] = _function if isinstance(_function, str) else str()
        item[ARGS] = _args if isinstance(_function, tuple) else tuple()
        item[KWARGS] = _kwargs if isinstance(_function, dict) else dict()
        return Operation(item)

    @staticmethod
    def generate_item_from_list(input_list: Union[list, tuple]):
        """ Generates an operation from an ordered list

        :param input_list: [agent, func, args, kwargs]
        :return: Operation object
        """
        if len(input_list) != 4:
            return None
        item = {AGENT: input_list[0], FUNC: input_list[1], ARGS: input_list[2], KWARGS: input_list[3]}
        return Operation(item)

    @staticmethod
    def generate_item_from_dict(input_dict: dict):
        """ Generates an operation from a dictionary

        :param input_dict: python dictionary
        :return: Operation object
        """
        return Operation(input_dict)
