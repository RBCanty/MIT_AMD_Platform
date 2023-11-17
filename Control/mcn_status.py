""" mcn_status module

Contains:
 - RetObj : Return Objects (all methods called from the MCN must return a RetObj)
 - Fault : When RetObj indicate failure or a fault occurs on the platform, Faults are generated to handle reporting
 - Status : Stores system- and platform-level information on the MCN

 @author: Ben C
"""
import json
import time
from copy import deepcopy
from datetime import datetime
from functools import wraps
from pprint import pformat
from threading import Lock, Event
from typing import Union, Any

from custom_classes import Narg
from custom_exceptions import StatusAccessTimeout, FaultValueError
from mcn_logging_manager import system_log

S_NETWORK = "network"
S_MAJ_SYS = "major_systems"
S_ALL_SYS = "all_systems"
S_STRUCTR = "structure"
S_AUXILIS = "auxiliaries"
S_DATABASES = "databases"
S_CAPACITY = "capacities"
S_SYSTEM_ = "system"
S_NAME = "name"
S_FAULTS = "error_and_warnings"
S_CHECKPOINTS = "checkpoints"
DS_NAME = [S_SYSTEM_, S_NAME]
DS_FAULTS = [S_SYSTEM_, S_FAULTS]
DS_CHECKPOINTS = [S_SYSTEM_, S_CHECKPOINTS]
N_CONNECTIVITY = 'Connectivity'
N_MODE = 'Mode'
N_MEMBERS = 'Members'
DN_MEMBERS = [S_SYSTEM_, S_NETWORK, N_MEMBERS]
DN_CONNECTIVITY = [S_SYSTEM_, S_NETWORK, N_CONNECTIVITY]
DN_MODE = [S_SYSTEM_, S_NETWORK, N_MODE]
V_ONLINE = 'Online'
V_OFFLINE = 'Offline'
V_GOOD = 'Good'
V_BUSY = 'Busy'
V_PROBLEM = 'Problem'
V_FATAL = 'Fatal'
TIME_FORMAT = "%m/%d/%Y %H:%M:%S"
LEVEL_MAP = {
    V_GOOD: 4,
    V_BUSY: 3,
    V_PROBLEM: 2,
    V_FATAL: 1
}

C_ACTIVE_MODE = 'active_mode'
C_QUEUE_BLOCK = 'queue_blacklist'
C_SAFE_MODE = 'safe_mode'
C_EXPIRED = 'expired_flag'
C_TASK_QUEUE = 'tasqueue'

# If Slack should be used and when
USE_SLACK = False
SLACK_LEVELS = [V_FATAL, ]

# Used for Message passing (all systems which are nodes in the network)
MEMBER_LIST = ["AH", "LC", "SP", "MC", "_A"]

# Numerical Map
SYS_ENUM = {
    "MC": 0,
    "Ra": 1,
    "AH": 2,
    "Lh": 3,
    "Pr": 4,
    "LC": 5,
    "Sw": 6,
    "SP": 7,
    "Fs": 8,
    "Ss": 9,
    "Th": 10,
    "MC.__": 11,
}

# Which subsystem is the silent/bypass subsystem
SCH_BYPASS = "MC.__"

# Master Configuration "file" for the MCN
MCN_CFG = {
    S_NETWORK: {"IP": "12.34.567.89", "PORT": 1324},  # Get from computer running the server
    S_MAJ_SYS: ["MC", "AH", "LC", "SP", "_S", "_U"],
    S_ALL_SYS: ["MC", "Ra",
                "AH", "Lh",
                "LC", "Sw",
                "SP", "Fs", "Ss", "Pr", "Th",
                "__", "_S", "_U", "ui"],
    S_STRUCTR: {"MC": ["Ra", "MC.__"],
                "AH": ["Lh"],
                "SP": ["Fs", "Ss", "Pr", "Th"],
                "LC": ["Sw"],
                "__": list(),
                "_S": list(),
                "_U": ["ui", ]
                },
    S_AUXILIS: {"MC": {"ra_sock": ["123.456.7.8", 8082]},  #
                "AH": {},
                "SP": {"door_com": "COM8",
                       "pcb_com": "COM9",
                       "omega_com": "COM6",
                       "valve_com": "COM12",
                       "lpx_com": "COM10",
                       "lift_com": "COM20",
                       "solar_omega": "COM15",
                       "solar_button": "COM21",
                       "solar_valve": None,
                       "anem_valve": "COM22"},
                "LC": {},
                "__": dict(),
                "_S": dict(),
                "_U": dict()
                },
    S_SYSTEM_: {S_NAME: "",
                S_FAULTS: list(),
                S_CHECKPOINTS: dict(),
                S_NETWORK: {N_CONNECTIVITY: V_OFFLINE, N_MODE: V_ONLINE, N_MEMBERS: list()}
                },
    S_DATABASES: {
        # name: [url, port]
        "Platform": ['http://12.34.567.89:12345', '67890'],  # Get from computer hosting platform database
        "Data": ["http://12.34.567.890:12345", '67890'],  # Get from computer hosting the data database
    },
    S_CAPACITY: {
        # Example stub for using the ResourceCapacityManager to avoid Busy-responses
        # would need to import ResourceCapacityManager as RCM from custom_classes.py
        # "liquid_handler_heater_shakers": RCM(
        #     capacity=6,
        #     context_on={
        #         'agent': "Lh",
        #         'operation': "start_stop_heater_shaker",
        #         'details': {
        #             'power': "on"
        #         }
        #     },
        #     context_off={
        #         'agent': "Lh",
        #         'operation': "start_stop_heater_shaker",
        #         'details': {
        #             'power': "off"
        #         }
        #     },
        # )
    }
}

# Timeout for waiting for a Status to unlock (seconds)
TIMEOUT = 5


def s_dir(parent_system: str) -> list:
    """ A system directory lookup

    :param parent_system: A system name ["MC", "AH", "LC", "SP", "_S", "_U"]
    :return: [parent_system, ] + MCN_CFG[S_STRUCTR].get(parent_system, list())
    """
    if parent_system not in MCN_CFG[S_MAJ_SYS]:
        raise ValueError(f"s_dir only defined for major systems, not '{parent_system}'")
    return [parent_system, ] + MCN_CFG[S_STRUCTR].get(parent_system, list())


def rs_dir(child_system: str) -> str:
    """ A reverse system directory lookup

    :param child_system: The (sub)system being searched
    :return: The Major system to which it belongs ("__" if not found)
    """
    if child_system == "_A":
        return "_A"
    if child_system in MCN_CFG[S_MAJ_SYS]:
        return child_system
    for maj_sys in MCN_CFG[S_STRUCTR].keys():
        if child_system in MCN_CFG[S_STRUCTR][maj_sys]:
            return maj_sys
    return "__"


# #### CLASS RETOBJ #### # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
class RetObj:
    """ A collection of completion, location, level, and information for an operation

    MCN Operations expect return objects, so the results of the operation can be properly reported and interpreted
    """
    def __init__(self, completion, *, location="Unknown", level=V_GOOD, data="No info provided"):
        """ Creates a Return Object

        USE NOTE: The data field for a fault should be a short, general description.  Use the logging functionality
        to provide details of the fault.

        :param completion: If the *operation* was completed successfully (False will generate a Fault)
        :param location: (If completion=False) Where the problem occurred [the system code, string]
        :param level: (If completion=False) Level of the fault [V_GOOD, V_BUSY, V_PROBLEM, V_FATAL]
        :param data: Information on the success/failure
        """
        self.completion = completion
        self.location = location
        self.level = level
        self.data = data

    @classmethod
    def complete(cls, data="NA"):
        """ Creates a return object for a successful operation

        :param data: Any information which may be needed
        :return: RetObj(True, data=data)
        """
        return RetObj(True, data=data)

    @classmethod
    def incomplete(cls, location, level, data):
        """ Creates a return object for a failed operation

        :param location: Where the problem occurred [the system code, string]
        :param level: Level of the fault [V_GOOD, V_BUSY, V_PROBLEM, V_FATAL]
        :param data: Pithy details on the fault (use Logging for more detail)
        :return: RetObj(False, location=location, level=level, data=data)
        """
        return RetObj(False, location=location, level=level, data=data)

    @classmethod
    def loads(cls, str_rep):
        """ Provides a Return object from a string

        :param str_rep: A json-readable string representation of the RetObj
        :return: RetObj via complete(d['data']) or incomplete(d['location'],
        d['level'], d['data']) [d = json.loads(str_rep)]
        """
        d = json.loads(str_rep)
        try:
            if d['completion']:
                return RetObj.complete(d['data'])
            else:
                return RetObj.incomplete(d['location'], d['level'], d['data'])
        except KeyError:
            system_log.exception(f"Return Object failed to be loaded from string: '{str_rep}'")
            return RetObj.incomplete("Unknown", V_PROBLEM, f"Received ill-formatted return object")

    def dumps(self):
        """ Provides a string representation of a return object

        :return: json.dumps({...})
        """
        if self.completion:
            return json.dumps({'completion': self.completion, 'data': self.data})
        else:
            return json.dumps({"completion": self.completion, "location": self.location,
                               "level": self.level, "data": self.data})

    def __repr__(self):
        return self.dumps()

    @classmethod
    def loadd(cls, dict_rep):
        """ Provides a Return Object loaded from a dictionary

        :param dict_rep: A dictionary representation of
        :return: RetObj via complete() or incomplete()
        """
        try:
            if dict_rep['completion']:
                return RetObj.complete(dict_rep['data'])
            else:
                return RetObj.incomplete(dict_rep['location'], dict_rep['level'], dict_rep['data'])
        except KeyError:
            system_log.exception(f"Return Object failed to be loaded from dict:\n{pformat(dict_rep)}")
            return RetObj.incomplete("Unknown", V_PROBLEM, f"Received ill-formatted return object")

    def dumpd(self):
        """ Provides a dictionary representation of the return object

        :return: {"completion": self.complete, "location": self.location, "level": self.level, "data": self.data}
        """
        return {"completion": self.complete, "location": self.location, "level": self.level, "data": self.data}

    def __dict__(self):
        return self.dumpd()

    def __str__(self):
        """ Human-readable string representation of RetObj for GUI

        :return: "({self.completion}[, {self.location}, {self.level}], {self.data})"
        """
        if self.completion:
            return f"({self.completion}, {self.data})"
        else:
            return f"({self.completion}, {self.location}, {self.level}, {self.data})"

    def __deepcopy__(self, memodict={}):  # noqa # This is the recommended syntax for deploying deepcopy
        """ Provides a dictionary

        Implemented for json

        :param memodict: Needed by deepcopy
        :return: A deep copy of the underlying dictionary
        """
        return RetObj(**deepcopy(self.dumpd()))


# #### CLASS FAULT #### # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
class Fault:
    """ MCN Fault Object

    Stores the location, level, info, and time of a fault on the platform
    """
    def __init__(self, location, level, data: str, timestamp=None, queue=None, **__):
        """ Creates a Fault object

        :param location: Location of the object [the system code, string]
        :param level: Level of the fault [V_GOOD, V_BUSY, V_PROBLEM, V_FATAL]
        :param data: Supplementary information [string]
        :param timestamp: Time of the fault [str: using TIME_FORMAT, datetime] (Default: use datetime.now())
        :param queue: A queue name with which the fault is associated
        """
        # Load location
        self.location = location
        if self.location not in MCN_CFG[S_ALL_SYS]:
            raise ValueError(f"Fault location must be a recognized system {MCN_CFG[S_ALL_SYS]} "
                             f"not '{self.location}'")

        # Load level
        self.level = level
        if self.level not in [V_GOOD, V_BUSY, V_PROBLEM, V_FATAL]:
            raise ValueError(f"Fault level must be {[V_GOOD, V_BUSY, V_PROBLEM, V_FATAL]} not '{self.level}'")

        # Load data
        self.data = data
        if not isinstance(self.data, str):
            raise ValueError(f"Fault data must be a string not '{data}' ({type(data)})")

        # Load timestamp
        self.timestamp = timestamp
        if timestamp is None:
            self.timestamp = datetime.now().strftime(TIME_FORMAT)
        elif isinstance(timestamp, datetime):
            self.timestamp = timestamp.strftime(TIME_FORMAT)
        elif isinstance(timestamp, str):
            self.timestamp = timestamp
        else:
            raise ValueError(f"A fault requires a timestamp which is a string or datetime object, not "
                             f"'{timestamp}' ({type(timestamp)})")

        self.queue = queue

    @classmethod
    def cast2fault(cls, representation):
        """ Converts a valid representation of a Fault into a real Fault Object

        Handles:
         - Fault
         - dict
         - list
         - tuple
         - str
         - :class:`RetObj`

        :param representation: A valid representation of a Fault
        :return: A Fault Object
        :raises ValueError, FaultValueError: When representation is not of a valid type or when the inputs forwarded to
          the constructor are invalid
        :raises JSONDecodeError: When a string representation cannot be parsed
        """
        if isinstance(representation, Fault):
            return representation
        if isinstance(representation, dict):
            return Fault.loadd(representation)
        if isinstance(representation, list) or isinstance(representation, tuple):
            return Fault(*representation)
        if isinstance(representation, str):
            return Fault.loads(representation)
        if isinstance(representation, RetObj):
            return Fault(location=representation.location, level=representation.level, data=representation.data)
        raise ValueError(f"Was not able to cast '{representation}' ({type(representation)}) to Fault")

    def __lt__(self, other):
        """ Compares Faults for sorting

        :param other: A Fault object
        :return: True if self is less than other, False otherwise
        """
        if not (self.location == other.location):
            return self.location < other.location
        # Location is equal
        if not (LEVEL_MAP.get(self.level, 5) == LEVEL_MAP.get(other.level, 5)):
            return LEVEL_MAP.get(self.level, 5) < LEVEL_MAP.get(other.level, 5)
        # Location and level are equalOptionMenu
        if not(self.timestamp < other.timestamp):
            return self.timestamp < other.timestamp
        # Location and level and timestamp are equal
        if not(self.__hash__() < other.__hash__()):
            return self.__hash__() < other.__hash__()
        # They are the same
        return False

    def __eq__(self, other):
        """ Compares Faults for sorting and for lookup

        Does not check for equivalence of associated queue names and the timestamps are allowed a buffer window

        :param other: A Fault object
        :return: True if self is equivalent than other, False otherwise
        """
        if not isinstance(other, Fault):
            return False
        loc = self.location == other.location
        lev = self.level == other.level
        dat = self.data == other.data
        my_time = datetime.strptime(self.timestamp, TIME_FORMAT)
        others_time = datetime.strptime(other.timestamp, TIME_FORMAT)
        tim = abs((my_time - others_time).total_seconds()) < 15
        # enough will be considered the same Fault
        # Neglected matching the name of the associated queue to prevent 2 queues from reporting the same fault
        #  it also makes clearing said faults easier.
        return loc and lev and dat and tim

    @classmethod
    def from_retobj(cls, retobj: RetObj, timestamp=None, queue=None):
        """ Creates a fault from a Return Object

        :param retobj: A (presumably failed) return object
        :param timestamp: Time of the fault [str: using TIME_FORMAT, datetime] (Default: use datetime.now())
        :param queue: Name of the queue associated with the fault
        :return: new Fault object
        """
        if not isinstance(retobj, RetObj):
            raise FaultValueError(f"Fault.from_retobj expected a RetObj not '{type(retobj)}' as first argument")
        try:
            return Fault(retobj.location, retobj.level, str(retobj.data), timestamp, queue)
        except AttributeError as ae:
            system_log.exception(f"Fault could not be generated from RetObj, RetObj ill-formatted: '{retobj}'")
            raise FaultValueError(*ae.args)

    def to_retobj(self):
        """ Casts a Fault to a Return Object

        :return: RetObj.incomplete(self.location, self.level, self.data)
        """
        return RetObj.incomplete(self.location, self.level, self.data)

    @classmethod
    def loads(cls, str_rep):
        """ Creates a Fault object from a string representation

        :param str_rep: A json-readable string representation of a Fault Object
        :return: Fault.loadd(json.loads(str_rep))
        """
        return Fault.loadd(json.loads(str_rep))

    def dumps(self):
        """ Creates a string representation of a Fault object

        :return: json.dumps(self.dumpd())
        """
        return json.dumps(self.dumpd())

    def __repr__(self):
        return self.dumps()

    @classmethod
    def loadd(cls, dict_rep):
        """ Creates a new Fault object from a dictionary

        :param dict_rep: A dictionary object with the correct keys [location, level, data[, timestamp]]
        :return: Fault(**dict_rep)
        """
        if not isinstance(dict_rep, dict):
            raise FaultValueError(f"Fault.loadd expected a dict not '{type(dict_rep)}' as only argument")
        try:
            return Fault(**dict_rep)
        except TypeError as te:
            system_log.exception(f"Fault could not be generated from dictionary, dictionary ill-formatted:\n"
                                 f"'{pformat(dict_rep)}'")
            raise FaultValueError(*te.args)

    def dumpd(self):
        """ Creates a dictionary representation of a Fault objecty

        :return: {"location": self.location, "level": self.level, "data": self.data, "timestamp": self.timestamp}
        """
        return {"location": self.location, "level": self.level, "data": self.data,
                "timestamp": self.timestamp, 'queue': self.queue}

    def __dict__(self):
        return self.dumpd()

    def __str__(self):
        """ Human-readable string representation of Fault (for GUI)

        :return: f"{self.location} - {self.level} - {self.timestamp} - [self.queue]{self.data}"
        """
        q_insert = f"[{self.queue}] " if self.queue else ''
        return f"{self.location} - {self.level} - {self.timestamp} - {q_insert}{self.data}"

    def __hash__(self):
        return hash((self.location, self.level, self.timestamp, self.data, self.queue))

    def __deepcopy__(self, memodict={}):  # noqa # This is the recommended syntax for deploying deepcopy
        """ Provides a dictionary

        Implemented for json

        :param memodict: Needed by deepcopy
        :return: A deep copy of the underlying dictionary
        """
        return Fault.loadd(deepcopy(self.dumpd()))


# #### CLASS Checkpoint #### # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
CDX = {'completion': 0, 'location': 1, 'level': 2, 'data': 3, 'queue': 4, 'step': 5, 'operation': 6, 'task_key': 7}


class Checkpoint(Event):
    """ Dedicated object for the monitoring of response-expecting operations of the platform

    Utilizes Event functionality to immediately trigger the handling of a return object
    A raised Checkpoint will not block scheduler but will persist such that it can be "revived" by the user to await
      a post-timeout response
    Raised Checkpoints shall have 3 'completion' values:
      (True) indicating a successful completion, (False) indicating a fault, and (None) indicating either timeout or
      the user raising the event (bypassing the Checkpoint)
    """
    def __init__(self, completion: Union[bool, None], location: str, level: Union[str, None], data: Any,
                 queue: str = None, step: Union[str, int] = None, operation: str = None, task_key: str = None,
                 *_, **__):
        """ Creates a Checkpoint object

        Used to await the completion of an operation across systems.

        Does not deploy validation on parameters

        :param completion: True if complete, False if Fault, or None incomplete/voided
        :param location: The name of the Agent
        :param level: (If completion is False) the level of the associated Fault
        :param data: Additional information
        :param queue: (Optional) The associated Queue document name
        :param step: (Optional) The associated Queue step number
        :param operation: (Optional) The associated operation name
        :param task_key: (Optional) The associated task key
        """
        super(Checkpoint, self).__init__()
        self.memory = [completion, location, level, data, queue, step, operation, task_key]
        self._backup = Event()

    def __del__(self):
        try:
            if not self.is_set():
                self.set()
        except:  # noqa
            pass
        try:
            if not self._backup.is_set():
                self._backup.set()
        except:  # noqa
            pass

    @classmethod
    def load_from_string(cls, string_repr):
        """ Converts a string (json) to a Checkpoint

        :param string_repr: A json-formatted string (representing a list/tuple (or dict)) of the constructor arguments
        :return: A Checkpoint object
        :raises JSONDecodeError: If string cannot be parsed
        :raises TypeError: If missing an argument for constructor or arguments cannot be unpacked
        """
        try:
            loaded_repr = json.loads(string_repr)
            if isinstance(loaded_repr, dict):
                return Checkpoint(**loaded_repr)
            else:
                return Checkpoint(*loaded_repr)
        except json.JSONDecodeError:
            return None

    def __getattr__(self, item):
        if isinstance(item, int):
            return self.memory[item]
        try:
            return self.memory[CDX[item]]
        except (IndexError, KeyError):
            return None

    def __setattr__(self, key, value):
        if key in CDX:
            self.memory[CDX[key]] = value
        else:
            super(Checkpoint, self).__setattr__(key, value)

    def unpack(self) -> list:
        """ Provides a copy (deepcopy) of the Checkpoint fields

        :return: A list representation of the Checkpoint
        """
        return deepcopy(self.memory)

    def update_checkpoint(self, **kwargs):
        """ Updates keyed values then triggers the Event

        :keyword completion:
        :keyword location:
        :keyword level:
        :keyword data:
        :keyword queue: (Optional)
        :keyword step: (Optional)
        :keyword operation: (Optional)
        :keyword task_key: (Optional)
        :return:
        """
        for k, v in kwargs.items():
            if k in CDX:
                self.memory[CDX[k]] = v
        self.set()
        self._backup.set()

    def __dict__(self):
        return {k: self.memory[v] for k, v in CDX.items()}

    def items(self):
        return self.__dict__().items()

    def wait_for_update(self, has_timeout=None,
                        get_time_est=lambda: 60, get_inw_times=lambda x: (max(x//4, 1), x), update_time_est=None):
        """ Blocking wait function for a checkpoint's update

        :param has_timeout: True if timeout should be checked, false if infinite loop
        :param get_time_est: A function which returns the initial estimate of the wait time
        :param get_inw_times: A function which returns the interval and wait times given est_time as an argument
        :param update_time_est: A function which updates the estimate of the wait time
        :return: (-1) The checkpoint is raised without update, (0) timeout, (1) complete, (2) incomplete
        """
        if update_time_est is None:
            update_time_est = get_time_est
        c_timer = datetime.now()
        est_time = get_time_est()
        while True:
            interval, wait_time = get_inw_times(est_time)
            if self.wait(interval):
                if self.completion is None:
                    return -1
                break
            else:
                if has_timeout and (datetime.now() - c_timer).total_seconds() > wait_time:
                    self.set()
                    return 0
            est_time = update_time_est()

        if self.completion:
            return 1
        else:
            return 2

    def wait_for_tardy(self, queue=None, step=None, operation=None):
        """ Used to wait 1 day for a checkpoint to be updated

        :param queue: (Optional) What the value of checkpoint.queue should be, if known
        :param step: (Optional) What the value of checkpoint.step should be, if known
        :param operation: (Optional) What the value of checkpoint.operation should be, if known
        :return: completion, queue, step, operation
        :raises TimeoutError: After waiting 1 day
        :raises ValueError: If the provided arg values can not be reconciled with those of the Checkpoint itself
        """
        def validate(internal, external):
            if (external is None) and (internal is None):
                raise ValueError("No checkpoint value can be found")
            elif external is None:
                return internal
            elif internal is None:
                return external
            elif internal != external:
                raise ValueError(f"Checkpoint values do not match: '{external}' (provided) vs '{internal}' (internal)")
            else:
                return internal  # they are equal

        validate(self.queue, queue)
        validate(self.step, step)
        validate(self.operation, operation)

        if not self._backup.wait(24*60*60):
            raise TimeoutError(f"Checkpoint {self.__str__()} is over 1 day tardy, abort silent monitoring")

        # NOTE: For _backup to set, update_checkpoint() must have been called, thus these may have changed
        return self.completion, self.queue, self.step, self.operation

    def __str__(self):
        return str(self.memory)

    def __repr__(self):
        return json.dumps(self.memory)

    def __eq__(self, other):
        """ Compared two Checkpoints

        If a task_key is provided for both, only checks task_key values.

        If a queue name, a step number, and operation name are provided for both, only checks said parameters. (The step
        numbers are cast to integers for comparison only)

        Otherwise, check location, level, and data for equivalence.

        :param other: A Checkpoint
        :return: True if equivalent, False otherwise
        """
        if not isinstance(other, Checkpoint):
            return False

        # completion, location, level, data, (queue, step, operation, task_key)
        task_key_not_none = (self.task_key is not None) and (other.task_key is not None)
        queue_not_none = (self.queue is not None) and (other.queue is not None)
        step_not_none = (self.step is not None) and (other.step is not None)
        operation_not_none = (self.operation is not None) and (other.operation is not None)

        if task_key_not_none:
            return self.task_key == other.task_key
        try:
            if queue_not_none and step_not_none and operation_not_none:
                queue_match = self.queue == other.queue
                step_match = int(self.step) == int(other.step)
                operation_match = self.operation == other.operation
                return queue_match and step_match and operation_match
        except (ValueError, TypeError):
            pass

        location_match = self.location == other.location
        level_match = self.level == other.level
        data_match = self.data == other.data
        return location_match and level_match and data_match

    def __deepcopy__(self, memodict={}):  # noqa # This is the recommended syntax for deploying deepcopy
        """ Provides a dictionary

        Implemented for json

        :param memodict: Needed by deepcopy
        :return: A deep copy of the underlying dictionary
        """
        return Checkpoint(**deepcopy(self.__dict__()))


# #### CLASS STATUS #### # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
class Status:
    """ MCN Status Object

    Repository for platform-relevant information as well as system-level information

    (See: MCN_CFG)
    """
    def __init__(self, name):
        """ Creates a Status object

        :param name: Name of the system being represented, should be a major system name
        """
        self.name = name
        self._memory: dict = deepcopy(MCN_CFG)
        self._memory[S_SYSTEM_][S_NAME] = self.name
        self._lock = Lock()

    def via_lock(func):
        """ Wrapper function for the lock/try-foo/finally-unlock workflow

        Implements a timeout of TIMEOUT for the lock

        :param func: Wrapped function
        :return: Return of Wrapped function
        :raise StatusAccessTimeout: On Timeout
        """
        @wraps(func)
        def wrapper_function(self, *args, **kwargs):
            if self._lock.acquire(timeout=TIMEOUT):
                time.sleep(0.001)
                try:
                    if func.__name__ == 'modify_checkpoint':
                        print("update_checkpoint")
                        print(self)
                        print(args)
                        print(kwargs)
                    return func(self, *args, **kwargs)
                finally:
                    self._lock.release()
                    time.sleep(0.001)
            else:
                raise StatusAccessTimeout("Timeout when fetching MCN Status Object")
        return wrapper_function

    def __enter__(self):
        """ Implemented to allow WITH statement use

        Implements a timeout of TIMEOUT for the lock

        :return: The protected object
        :raise StatusAccessTimeout: On Timeout
        """
        if self._lock.acquire(timeout=TIMEOUT):
            return self._memory
        else:
            raise StatusAccessTimeout("Timeout when fetching MCN Status Object")

    def __exit__(self, exc_type, exc_value, tb):
        """ Implemented to allow WITH statement use

        :param exc_type: Exception type
        :param exc_value: Exception value
        :param tb: Traceback object
        :return: None
        """
        if exc_type is not None:
            system_log.exception("MCN Status __exit__ has encountered an exception")
        if self._memory is None:
            raise ValueError(f"self._memory is None from an __exit__, caller destroyed object")
        self._lock.release()

    def __dict__(self):
        """ Provides the protected dictionary object at the core of a Status object

        :return: self._memory
        """
        return self._memory

    @via_lock
    def get(self, key: list, default: Any = Narg()):
        """ Generalized getter

        :param key: A *list* of keys (for single values use [key, ] syntax
        :param default: Default value to return if key not present
        :return: value of self._memory[...] or default
        """
        pointer = self._memory
        for k in key:
            try:
                pointer = pointer[k]
            except KeyError as ke:
                if not isinstance(default, Narg):
                    return default
                else:
                    raise ke
        return pointer

    def __getitem__(self, key):
        """ [] getter for Status

        Use: Status[key1, key2, ..., keyN] not Status[key1][key2][...][keyN]

        :param key: A tuple of keys generated from [] access syntax
        :return: self.get(list(key))
        """
        return self.get(list(key))

    @via_lock
    def set(self, key: list, value):
        """ Generalized setter

        Creates entry if one doesn't exist

        :param key: A *list* of keys (for single values use [key, ] syntax
        :param value: Value to assign to key
        :return: None
        """
        pointer = self._memory
        final_key = key.pop()
        for k in key:
            if hasattr(pointer, 'keys'):
                pointer[k] = pointer[k] if k in pointer.keys() else dict()
                pointer = pointer[k]
            else:
                pointer[k] = pointer[k] if k < len(pointer) else list()
                pointer = pointer[k]
        if final_key not in pointer:
            system_log.warning(f"MCN Status object given unexpected key '{final_key}'")
        pointer[final_key] = value

    def __setitem__(self, key, value):
        """ [] setter for Status

        Use: Status[key1, key2, ..., keyN] not Status[key1][key2][...][keyN]

        :param key: A tuple of keys generated from [] access syntax
        :param value: Value being assigned
        :return: None
        """
        return self.set(list(key), value)

    def keys(self):
        """ Recalls keys from self._memory

        :return: self._memory.keys(
        """
        return self._memory.keys()

    def get_short_copy(self) -> dict:
        """ Provides the system-level details of a system

        :return: self._memory[S_SYSTEM_]
        """
        return self._memory[S_SYSTEM_]

    @via_lock
    def update(self, other=None, **others):
        """ Updater for the self._memory field

        calls dict.update()

        :param other: A dictionary via dict object
        :param others: A dictionary via keyword arguments
        :return: None
        """
        if other is None:
            other = {}
        self._memory.update(**other, **others)

    @via_lock
    def short_update(self, other: dict):
        """ Updater for the system-level self._memory[S_SYSTEM_] field

        calls dict.update()

        :param other: A dictionary via dict object
        :return: None
        """
        others_faults: list = other[S_FAULTS]

        # clear old data from the child
        for fault in self._memory[S_SYSTEM_][S_FAULTS]:
            if fault.location in s_dir(other[S_NAME]):
                while fault in self._memory[S_SYSTEM_][S_FAULTS]:
                    self._memory[S_SYSTEM_][S_FAULTS].remove(fault)
        # Load in new data from the child
        for fault in others_faults:
            if fault not in self._memory[S_SYSTEM_][S_FAULTS]:
                self._memory[S_SYSTEM_][S_FAULTS].append(fault)

    # # # # # # # # Getters/Setters # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    def get_ip_and_port(self):
        """ Getter for socket information

        :return: (ip: str, port: int)
        """
        return self._memory[S_NETWORK]["IP"], self._memory[S_NETWORK]["PORT"]

    @via_lock
    def modify_network_settings(self, ip: str = None, port: int = None, *, pair: tuple = None):
        """ Changes the IP and PORT values

        :param ip: (str) the IP for a socket
        :param port: (int) the port number for a socket
        :param pair: (Usurps ip and port) a tuple of the ip and port
        :return: None
        """
        if pair:
            self._memory[S_NETWORK]["IP"] = pair[0]
            self._memory[S_NETWORK]["PORT"] = pair[1]
        else:
            self._memory[S_NETWORK]["IP"] = ip
            self._memory[S_NETWORK]["PORT"] = port

    @via_lock
    def network_is(self, conn, mode):
        """ Sets the network connectivity and mode

        :param conn: The connectivity (*is* connected), accepts V_ONLINE, V_OFFLINE, True, and False
        :param mode: The mode (*should* be connected), accepts V_ONLINE, V_OFFLINE, True, and False
        :return: None
        """
        if conn in [V_ONLINE, True]:
            conn = V_ONLINE
        elif conn in [V_OFFLINE, False]:
            conn = V_OFFLINE
        else:
            raise ValueError(f"Connection must be {V_OFFLINE}/False or {V_ONLINE}/True--not {conn}")
        if mode in [V_ONLINE, True]:
            mode = V_ONLINE
        elif mode in [V_OFFLINE, False]:
            mode = V_OFFLINE
        else:
            raise ValueError(f"Mode must be {V_OFFLINE}/False or {V_ONLINE}/True--not {mode}")
        self._memory[S_SYSTEM_][S_NETWORK][N_CONNECTIVITY] = conn
        self._memory[S_SYSTEM_][S_NETWORK][N_MODE] = mode

    @via_lock
    def network_conn_is(self, conn):
        """ Sets network connectivity

        Connectivity is a reporter of the current state. (For the prescriber of the connectivity state see Mode)

        :param conn: The current connection (True/V_ONLINE-->V_ONLINE) or (False/V_OFFLINE-->V_OFFLINE)
        :return: None
        """
        conn = V_ONLINE if conn is True else conn
        conn = V_OFFLINE if conn is False else conn
        if conn not in [V_ONLINE, V_OFFLINE]:
            raise ValueError(f"Network connectivity must be on/offline not '{conn}'")
        self._memory[S_SYSTEM_][S_NETWORK][N_CONNECTIVITY] = conn

    @via_lock
    def network_mode_is(self, mode):
        """ Sets network mode

        Mode is the prescriber of the network state.  (For the describer of the current network state see Connectivity)

        :param mode: The desired connect state (True/V_ONLINE-->V_ONLINE) or (False/V_OFFLINE-->V_OFFLINE)
        :return: None
        """
        mode = V_ONLINE if mode is True else mode
        mode = V_OFFLINE if mode is False else mode
        if mode not in [V_ONLINE, V_OFFLINE]:
            raise ValueError(f"Network mode must be on/offline not '{mode}'")
        self._memory[S_SYSTEM_][S_NETWORK][N_MODE] = mode

    def get_network_state(self):
        """ Fetches the connectivity and mode as a tuple

        :return: (self._memory[S_SYSTEM_][S_NETWORK][N_CONNECTIVITY], self._memory[S_SYSTEM_][S_NETWORK][N_MODE])
        """
        return self._memory[S_SYSTEM_][S_NETWORK][N_CONNECTIVITY], self._memory[S_SYSTEM_][S_NETWORK][N_MODE]

    @via_lock
    def update_network_members(self, member_list):
        """ Updates the members of the member list

        :param member_list: A list of members from the Server
        :return: None
        """
        member_list = [m for m in member_list if m in MCN_CFG[S_MAJ_SYS]]
        self._memory[S_SYSTEM_][S_NETWORK][N_MEMBERS] = member_list

    def get_subsystems(self):
        """ Provides the subsystems of a system

        :return: [ss for ss in self._memory[S_STRUCTR][self.name]]
        """
        return [ss for ss in self._memory[S_STRUCTR][self.name]]

    def is_faulty(self, subsystem=None):
        """ Reports if a (sub)system has active fatal faults

        :param subsystem: (defaults to the calling system) A subsystem to be investigated
        :return: True if any fatal faults are reported at subsystem, False otherwise
        """
        if subsystem is None:
            subsystem = self.name
        faults = self.get(DS_FAULTS, list())
        for f in faults:
            if (f.location == subsystem) and (f.level == V_FATAL):
                return True
        return False

    def get_faults(self, subsystem=None, up=False, *, _all=False) -> list:
        """ Provides faults reported for a (sub)system

        :param subsystem: (defaults to the calling system) A subsystem to be investigated
        :param up: False-Search this specific subsystem, True-Search the parent of this subsystem
        :param _all: True - return self.get(DS_FAULTS, list()), no matching (sub)system's location, supersedes others
        :return: [f for f in Faults if f.location == subsystem]
        """
        if _all:
            return self.get(DS_FAULTS, list())
        if subsystem is None:
            subsystem = self.name
        # _search = lambda x: rs_dir(x) if up else x
        def _search(x): rs_dir(x) if up else x
        faults = self.get(DS_FAULTS, list())
        return [Fault.cast2fault(f) for f in faults if _search(f.location) == subsystem]

    @via_lock
    def add_fault(self, fault: Union[Fault, RetObj, dict, str], queue=None):
        """ Adds a fault to the Status object

        ONLY CALL FROM SYNCHRONIZATION (otherwise other systems won't be aware of the fault)

        :param fault: A Fault (or object castable ot a Fault: RetObj, dict, str)
        :param queue: Allows the fault.queue field to be overwritten (used when creating from RetObj)
        :return: None
        """
        if not fault:
            system_log.info("Failed to add fault, was Falsey")
            return
        fault = Fault.cast2fault(fault)
        if queue:
            fault.queue = queue
        if fault not in self._memory[S_SYSTEM_][S_FAULTS]:
            self._memory[S_SYSTEM_][S_FAULTS].append(fault)

    @via_lock
    def remove_fault(self, fault: Union[Fault, RetObj, dict, str]):
        """ Adds a fault to the Status object

        ONLY CALL FROM SYNCHRONIZATION (otherwise other systems won't be aware of the fault)

        :param fault: A Fault (or object castable ot a Fault: RetObj, dict, str)
        :return: None
        """
        fault = Fault.cast2fault(fault)
        if not fault:
            system_log.info("Failed to remove fault, was Falsey")
            return
        while fault in self._memory[S_SYSTEM_][S_FAULTS]:
            self._memory[S_SYSTEM_][S_FAULTS].remove(fault)

    def is_busy(self, subsystem=None, bypass=None):
        """ Reports if a (sub)system has any active checkpoints

        :param subsystem: (defaults to the calling system) A subsystem to be investigated
        :param bypass: If subsystem == bypass, return False (ignores checkpoints)
        :return: True if any checkpoint's agent match the subsystem, False otherwise
        """
        if subsystem is None:
            subsystem = self.name
        if bypass and (subsystem == bypass):
            return False
        checkpoints = self.get(DS_CHECKPOINTS, dict())
        for _, ch in checkpoints.items():
            if ch.location == subsystem and not ch.is_set():
                return True

    def get_checkpoints(self, subsystem=None, *, _all=False) -> dict:
        """ Provides checkpoints working on a (sub)system

        :param subsystem: (defaults to the calling system) A subsystem to be investigated
        :param _all: Return all checkpoints (ignores subsystem)
        :return: {k:v for k,v in checkpoints.items() if v[1] == subsystem}
        """
        if _all:
            return self.get(DS_CHECKPOINTS, dict())
        if subsystem is None:
            subsystem = self.name
        checkpoints = self.get(DS_CHECKPOINTS, dict())
        return {k: v for k, v in checkpoints.items() if v[1] == subsystem}

    @via_lock
    def add_checkpoint(self, key: str, value: Checkpoint):
        """ Adds a checkpoint

        References :class:`Checkpoint`

        :param key: Message idempotency key
        :param value: A Checkpoint
        :return: None
        """
        if value.task_key is None:
            value.task_key = key
        self._memory[S_SYSTEM_][S_CHECKPOINTS].update({key: value})

    @via_lock
    def release_checkpoint(self, key):
        """ Removes a checkpoint

        :param key:
        :return:
        """
        checkpoints = self._memory[S_SYSTEM_][S_CHECKPOINTS]
        if hasattr(checkpoints.get(key, None), 'set'):
            checkpoints[key].set()
        try:
            del checkpoints[key]
        except KeyError:
            system_log.exception(f"Call to release checkpoint '{key}' moot, not present")

    def get_checkpoint(self, key) -> Checkpoint:
        """ Fetches a specific checkpoint

        References :class:`Checkpoint`

        :param key: Message idempotency key
        :return: A Checkpoint
        """
        return self._memory[S_SYSTEM_][S_CHECKPOINTS].get(key, None)

    @via_lock
    def modify_checkpoint(self, key, **kwargs):
        """ Edits the parameters of a Checkpoint by keyword

        Does not raise Exception if the Checkpoint does not exist

        References :class:`Checkpoint`

        :param key: A task key (used to look up the Checkpoint, but the task_key property is
          only modified by the task_key keyword in kwargs)
        :keyword completion: True if complete, False if Fault, or None incomplete/voided
        :keyword location: The name of the Agent
        :keyword level: (If completion is False) the level of the associated Fault
        :keyword data: Additional information
        :keyword queue: The associated Queue document name
        :keyword step: The associated Queue step number
        :keyword operation: The associated operation name
        :keyword task_key: The associated task key
        :return: None
        """
        if key not in self._memory[S_SYSTEM_][S_CHECKPOINTS]:
            system_log.info(f"Checkpoint update called on nonexistent checkpoint {{'{key}': {kwargs}}}")
            return
        sel_checkpoint: Checkpoint = self._memory[S_SYSTEM_][S_CHECKPOINTS][key]
        sel_checkpoint.update_checkpoint(**kwargs)

    # # # # # # # # Misc. # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    def __str__(self):
        """ Human-readable string representation of Status (for GUI)

        :return: pformat(self._memory[S_SYSTEM_])
        """
        return pformat(self._memory[S_SYSTEM_])

    def __repr__(self):
        return pformat(self._memory)

    # MUST be last thing in class! (Otherwise the compiler will freak out when making the @via_lock wrapper
    via_lock = staticmethod(via_lock)


# #### TESTING #### # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
if __name__ == '__main__':
    test = Status('MC')

    test.add_checkpoint("test_check_point_key",
                        Checkpoint(True, "_U", V_PROBLEM, "this is a test",
                                   "not_a_real_queue", "4", "test_check_point_key"))
    print("(1)", test.get(DS_CHECKPOINTS))

    temp = test.get_checkpoints(_all=True)
    print("(2)", temp)

    test_checkpoints_as_str = json.dumps(temp, default=lambda o: o.__dict__())
    print("(3)", test_checkpoints_as_str)

    test.release_checkpoint("test_check_point_key")
    print("(4)", test.get_checkpoints(_all=True))

    past_life = json.loads(test_checkpoints_as_str)
    print("(5)", past_life)
    for _k, _v in past_life.items():
        test.add_checkpoint(_k, Checkpoint(**_v))

    print("(6)", test.get_checkpoints(_all=True))
