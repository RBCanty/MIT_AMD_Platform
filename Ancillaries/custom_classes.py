"""
Special classes used on the platform as well as some helper methods for file IO

@author: Ben Canty
"""

import os.path
import traceback
from pprint import pprint, pformat
import time
import json
from os import listdir
from os.path import join, isdir, dirname, basename
from threading import Event
import schema as sch
from typing import Iterable, Callable


class SliceDict(dict):
    """ Mutation of dict (keyed by numbers [int] only)

    Allows multiple keys to be accessed like slices of a list
    """
    def __setitem__(self, key, value):
        try:
            int(key)
        except (TypeError, ValueError) as err:
            raise ValueError("Sliceable Dictionary must use integer indices") from err
        dict.__setitem__(self, int(key), value)

    def __getitem__(self, item):
        if isinstance(item, slice):
            start, stop = item.start, item.stop
            if start is None:
                start = min(list(self.keys()))
            if stop is None:
                stop = max(list(self.keys())) + 1
            return [dict.get(self, value, None) for value in range(start, stop)]
        else:
            item = int(item)
            return dict.__getitem__(self, item)


class CommandLineEnvironment:
    """
    For Debugging and Troubleshooting, emulates a command-line-like environment for controlling a system via its class
    """
    def __init__(self, child_class, extra_keywords=None):
        """
        Creates a rough command line emulator

        :param child_class: The Class being managed via command line
        :param extra_keywords: A list of function names to be ignored (defaults to dunders and base Object methods)
        """
        if extra_keywords is None:
            extra_keywords = []
        self.child = child_class
        self.extra = ['__del__', '__dict__', '__enter__', '__exit__', '__module__', '__weakref__']
        self.extra.extend(extra_keywords)
        self.available_functions = [x for x in dir(self.child) if x not in dir(object) + self.extra]
        pprint(self.available_functions)

    def run(self):
        """
        Creates a rough command line emulation.

        Formatting:
         - '-exit' to exit the emulation
         - '-help' to list the methods which can be called
         - 'method -help' to print the docstring of method
         - 'method kwargs' - kwargs are taken in json format (positional arguments not accepted)

        :return: None
        """
        while True:
            time.sleep(0.5)
            command = input("Enter Command (-exit to exit, -help for command list)\n>> ")
            command = command.strip()
            if command == "-exit":
                break
            elif command == "-help":
                pprint(self.available_functions)
                continue

            kwargs = dict()
            if " " in command:
                try:
                    command, args = command.split(" ", 1)
                except ValueError:
                    print("Arguments could not be separated from command")
                    continue

                command = command.strip()
                args = args.strip()

                if args == "-help":
                    kwargs = {"help": True}
                else:
                    try:
                        kwargs = json.loads(args)
                    except json.decoder.JSONDecodeError as jde:
                        traceback.print_tb(jde.__traceback__)
                        continue

            func = getattr(self.child, command, None)
            if func is None:
                continue

            if kwargs.get("help", False):
                print(getattr(func, "__doc__", "no __doc__ found"))
                try:
                    kwargs = json.loads(input("Enter kwargs as a json dictionary\n>> "))
                except json.decoder.JSONDecodeError as jde:
                    traceback.print_tb(jde.__traceback__)
                    continue
            print(f"Running {command} with arguments: {pformat(kwargs)}")
            try:
                ret_val = func(**kwargs)
                print(f"{command} returned:\n{pformat(ret_val)}")
            except Exception as e:
                print(f"Exception\n{repr(e)}\n")
                traceback.print_tb(e.__traceback__)


class Nop(object):
    """
    Dummy class that will do nothing

    But shouldn't throw syntax errors when used
    """
    def __call__(self, *_, **__): pass
    def __getattr__(self, *_, **__): return Nop()
    def __get__(self, *_, **__): return Nop()
    def __getitem__(self, *_, **__): return Nop()
    def __setattr__(self, *_, **__): return
    def __set__(self, *_, **__): return
    def __setitem__(self, *_, **__): return
    def __enter__(self, *_, **__): return Nop()
    def __exit__(self, *_, **__): pass


class Narg(Nop):
    """ A default argument class

    For when the default can't be None (the modern '...' object fulfils this role but is not compatible with all
    versions of python running MCN code)
    """
    pass


class Generic:
    """
    A generic object which can be monkey-patched (not recommended for main code, but useful for testing)
    """
    pass


class ParcelEvent(Event):
    """ Extension of Event class

    Used such that the Event can hold a value which can be read once the Event is raised (trying to read before being
    set will raise a ValueError).
    """
    def __init__(self):
        super(ParcelEvent, self).__init__()
        self.parcel = None

    def set_parcel(self, value):
        """ Sets the internal memory (self.parcel) to a value.

        :param value: The value stored (the Event does not have to be (un)set to assign this value, setting the value
          does not set the event, and it can be reassigned freely)
        :return: None
        """
        self.parcel = value

    def get_parcel(self):
        """ Gets the internal memory (self.parcel) if the event is set

        :return: (if Event is set) the value stored, (else) raises ValueError
        """
        if self.is_set():
            return self.parcel
        else:
            raise ValueError("Parcel cannot be retrieved until the Event is set")


class ResourceCapacityManager:
    """
    Contingency plan: If knowledge of how occupied certain resources are be required by scheduler, then this class would
    help keep the scheduler informed.

    At present, an overload of an instrument will return a Busy response, and the scheduler will devalue trying the same
    operation again for 5 minutes.

    While not used, it is referenced in commented-out code in the mcn_status.py file for the system configuration.

    In implementation, it is a counted semaphore where on/off calls increment/decrement a register which has a max and
    a minimum value.
    """
    def __init__(self, capacity: int, context_on: dict, context_off: dict,
                 parser: Callable = Narg(), init_register: int = 0):
        """ Creates a manager for the resource utilization of the platform.

        Provides a basic functionality for the Master Controller to be aware of a system having a limited capacity for
        an operation (e.g. number of heater shakers) such that it does not have to become aware of running out of a
        resource by the corresponding Agent giving a 'busy' reply.

        :param capacity: The maximum number of 'on' calls allowed
        :param context_on: A Scheme for when 'on' is called
        :param context_off: A Scheme for when 'off' is called
        :param parser: Logic function which takes (bool if context_on scheme is met, bool if context_off scheme is met)
          and returns a boolean (true - on, false - off), or None (on/off not applicabe)
          or raises (ValueError - both on and off)
        :param init_register: The value the register should be initialized to (default: 0)
        """
        self.capacity = capacity
        self.context_on = sch.Schema(context_on, ignore_extra_keys=True)
        self.context_off = sch.Schema(context_off, ignore_extra_keys=True)
        if isinstance(parser, Narg):
            self.parser = self.default_parser
        else:
            self.parser = parser
        self.register = init_register

    @staticmethod
    def str_or(values: Iterable):
        """ Creates an AND schema which will validate strings into lowercase versions found in values

        :param values: Acceptable values for a string to be
        :return: And(str, Use(str.lower), lambda s: s in values)
        """
        return sch.And(str, sch.Use(str.lower), lambda s: s in values)

    @staticmethod
    def int_range(low: int, high: int):
        """ Creates an AND schema which will check integers between lower (inclusive) and upper bounds (exclusive)

        :param low: Lower bound (inclusive)
        :param high: Upper bound (exclusive)
        :return: And(Use(int), lambda n: low <= n < high)
        """
        return sch.And(sch.Use(int), lambda n: low <= n < high)

    @staticmethod
    def optional(*args, **kwargs):
        """ Creates an OPTIONAL schema

        Wraps schema directly for Optional, mostly to avoid redundant import statements

        :param args: Given to Optional()
        :param kwargs: Given to Optional()
        :return: Optional(*args, **kwargs)
        """
        return sch.Optional(*args, **kwargs)

    def validate(self, operation):
        """ Invokes schema.validate() but catches SchemaError

        :param operation: The first positional argument to the context_on and context_off schema
        :return: (is_on, is_off) - booleans for if the respective schema were matched
        """
        is_on, is_off = False, False

        try:
            self.context_on.validate(operation)
        except sch.SchemaError:
            pass
        else:
            is_on = True

        try:
            self.context_off.validate(operation)
        except sch.SchemaError:
            pass
        else:
            is_off = True

        return is_on, is_off

    def manageable(self, operation, check_only=True):
        """ Checks and modifies the register based on the given operation matching (or nor matching) the on and off
        criteria.

        :param operation: An operation which is checked against the schema
        :param check_only: True - does not modify the register (is the operation allowed?), False - modifies the
          register (perform the operation.)
        :return: True - the operation is manageable, False - the operation is not manageable, None - the operation is
          not relevant to this resource manager
        :raises ValueError: The criteria are bad (both were met)
        """
        is_on, is_off = self.validate(operation)

        try:
            code = self.parser(is_on, is_off)
            if code is None:
                return None
            elif code:
                if self.register + 1 > self.capacity:
                    return False
                if not check_only:
                    self.register += 1
                return True
            else:
                if not check_only:
                    self.register -= 1
                    self.register = max(0, self.register)
                return True
        except ValueError as ve:
            raise ve

    @staticmethod
    def default_parser(is_on: bool, is_off: bool):
        """
        Together with self.manageable(), determines the behavior of the resource manager, converts the (is_on, if_off)
        tuple into a single, usable value.

        :param is_on: A boolean for if the context_on scheme is met
        :param is_off: A boolean for is the context_off scheme is met
        :return: True - on, False - off, None - neither
        :raises ValueError: When both are met
        """
        if is_on and is_off:
            raise ValueError("Both criteria are met, contexts poorly defined")
        if is_on:
            return True
        if is_off:
            return False
        return None


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def pathfinder(target: str, filename=None, max_iter=4) -> str:
    """ Searches for a path starting from the working directory and moving up.

    :param target: The directory being searched for
    :param filename: A file name that can be joined to the path once the path is found
    :param max_iter: Height of search
    :return: The directory's path
    """
    filepath = str()
    par_dir = dirname(__file__)
    _iter = 0
    while _iter <= max_iter:
        dirs = [basename(join(par_dir, d)) for d in listdir(par_dir) if isdir(join(par_dir, d))]
        if target in dirs:
            filepath = par_dir
            break
        par_dir = dirname(par_dir)
        _iter += 1
    else:
        pass  # Should there be special behavior for this case--where the iter maxes out
        # like it should return None right?
    if filename:
        return os.path.join(filepath, target, filename)
    else:
        return os.path.join(filepath, target)


def safe_open(path, mode):
    """
    Helper to replace open(path, mode) such that if the directory does not yet exist, it is created.

    :param path: The path to the file
    :param mode: the mode of basic IO's open() method
    :return:
    """
    if isinstance(path, list) or isinstance(path, tuple):
        path = os.path.join(*path)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    return open(path, mode)


def directory_check(path_to_check):
    """
    Helper to enforce a directory exists (for cases where safe_open may not be applicable,
    such as a third-party function that requests a path).

    :param path_to_check: The path to the directory
    :return: None
    """
    if not os.path.exists(path_to_check):
        os.makedirs(path_to_check)
        print('Made a new directory: ' + path_to_check)
    return


def grab_newest_folder(directory):
    """
    Searches a directory for the most recent (determined by os.path.getmtime) item

    :param directory: A path explorable by os.listdir
    :return: max([os.path.join(directory, d) for d in os.listdir(directory)], key=os.path.getmtime)
    """
    return max([os.path.join(directory, d) for d in os.listdir(directory)], key=os.path.getmtime)


if __name__ == '__main__':
    # print(grab_newest_folder(r"C:\Users\chemegrad2018\PycharmProjects\AMD_Control_Platform\venv\Shell_AH\Spark_API"))
    data = {'agent': "Lh",
            'operation': "start_stop_heater_shaker",
            'container': "reaction_plate",
            'details': {
                'power': "on",
                'temperature': 84.0,
                'rpms': 600,
                'is_paired': "yes",
                'schedule_time': 7200,
            }}
    schema = sch.Schema(
        {
            'agent': "Lh",
            'operation': "start_stop_heater_shaker",
            'details': {
                'power': "off"
            }
        },
        ignore_extra_keys=True
    )

    validated = schema.validate(data, )

    print(validated)
