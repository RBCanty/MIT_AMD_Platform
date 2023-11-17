""" task_management

Contains Coordinator

Coordinator coordinates tasks between DB and platform

Relevant Queue document fields
  * 'status': idle/in_progress/error/complete
  * 'dependency': parent /or/ [parents, ...]

The last operation of a queue (complete_queue) will change the status if necessary and move it to Historical

@author: Ben C
"""

from collections import namedtuple
import time
from copy import deepcopy
from datetime import datetime, timedelta
from pprint import pformat
from threading import Thread
from typing import List, Union

import custom_exceptions as cexc
import database_interface as dbi
import functionals
from ui_database_prompt import DatabasePrompt, tk
import mcn_status as mcs
import message as msgm
import operations as oprtn
import thread_safe_data as tsd
from constants import MSG_Q_IN, MSG_Q_OUT, CHILD_COM
from database_constants import *
from mcn_logging_manager import system_log
from robot_arm import LOCATION_CODEX
from task_logging_manager import tasking_log

# Necessary data used to track resource locations on platform such that collisions can be avoided
_LOCATION_CODEX = deepcopy(LOCATION_CODEX)
[_LOCATION_CODEX.pop(k) for k in list(_LOCATION_CODEX.keys()) if "liquid_handler" in k]

TRANSFER_OPERATIONS = ["move_wellplate", "Ss.stow", "Ss.request"]
LOCATION_CAPACITIES = {k: len(v) for k, v in _LOCATION_CODEX.items()}

# Some local "Data Types" for keeping track of scheduling information
Card = namedtuple("Card", ["q_name", "q_date", "q_status", "q_next", "q_agent"])
""" A quick reference for tasks containing properties of the queue and associated operation such as queue name, 
 creation date, queue status, the next operation to be run, said operation's agent """

Transfer = namedtuple("Transfer", ["container", "qn_nn", "destination", "start"])
""" A quick reference for transfer tasks containing properties of the container, queue name, the nickname of the 
 container, the destination, and the operation start time"""


class TimeEntry:
    """ Helper Class to store start and end times with validated types for comparison """
    def __init__(self, start: datetime, end: datetime = None):
        """ A quick reference for start-end time pairs for operations

        Times are automatically cast to datetime or NoneType

        :param start: The start time
        :param end: The end time
        """
        self.start = self.validator(start)
        self.end = self.validator(end)

    @staticmethod
    def validator(arg) -> Union[datetime, None]:
        """ Used to cast to datetime or None

        :param arg: A start or end time
        :return: A datetime object or None
        :raises ValueError: If arg could not be parsed
        """
        if isinstance(arg, datetime):
            return arg
        elif isinstance(arg, str):
            return datetime.strptime(arg, mcs.TIME_FORMAT)
        elif not arg:
            return None
        else:
            raise ValueError(f"Times must be a datetime, string of the form '{mcs.TIME_FORMAT}', or Falsy; not {arg}")

    @property
    def start(self) -> datetime:
        return self._start

    @start.setter
    def start(self, s):
        self._start = self.validator(s)

    @property
    def end(self) -> datetime:
        return self._end

    @end.setter
    def end(self, e):
        self._end = self.validator(e)

    def __str__(self):
        str_end = self._end.strftime(mcs.TIME_FORMAT) if self._end else self._end
        return f"TimeEntry(start={self._start.strftime(mcs.TIME_FORMAT)}, end={str_end})"


class Coordinator:
    """ Manages the workflow for the MCN

    my_coordinator = Coordinator(self.data_pipes, self.status)
    my_coordinator.run()
    """

    def __init__(self, pipes, status: mcs.Status, debug_mode=False):
        """ Creates a coordinator for the MCN

        :param pipes: Message queues and Child com
        :param status: Status object for MC
        :param debug_mode: (bool) For troubleshooting, prompts user before each action
        """
        self._pipes = pipes
        self.inbox = pipes[MSG_Q_IN]
        self.outbox = pipes[MSG_Q_OUT]
        self.child_com = pipes[CHILD_COM]
        self.status = status
        self.threads = tsd.ThreadSafeDataContainer(dict())
        self.queue_collection = dict()  # {queue_name: queue_doc}, {_: (name, date, next_step)}
        self.historical_collection = dict()  # {queue_name: queue_doc}
        self.register = list()
        self.debug_mode = debug_mode
        self.resource_ctrl = mcs.MCN_CFG[mcs.S_CAPACITY]
        self.penalty_box = dict()
        self.bypass = mcs.SCH_BYPASS
        self.polling_frequency = 7.5  # units: seconds
        self.mapping = {}
        self.connectivity = [True, 0]  # [are we connected, failed attempts]

    def run(self):
        """ Mainloop of Coordinator

        Spawns a processor_in thread and a janitor thread

        :return: None, upon a core thread dying
        """
        # startup
        tasking_log.info("Task manager has started")
        while not self.update(_raise=False):
            system_log.info("Coordinator startup delayed: Cannot connect to DB, request user Okay")
            user_input = functionals.quick_gui(title="User Intervention Requested",
                                               dialog="Coordinator cannot connect to Database\n"
                                                      "Press OK when resolved",
                                               ret_if_ok="resume")
            if user_input == 'resume':
                pass
            else:
                raise RuntimeError("User has exited prompt for Coordinator-DB communication error")

        # Activate processes
        t1 = Thread(target=self.processor, daemon=True)
        t2 = Thread(target=self.janitor, daemon=True)
        t1.start()
        t2.start()
        while True:
            if not t1.is_alive() or not t2.is_alive():
                msg_details = (item for item in [bool(not t1.is_alive()) * "processor",
                                                 bool(not t2.is_alive()) * "janitor"]
                               if item)
                system_log.warning(f"Coordinator died: check {', '.join(msg_details)}")
                tasking_log.info("Task manager has terminated")
                return
            time.sleep(5)

    def update(self, _raise=True):
        """ Pulls DB for most recent information

        :param _raise: If exceptions should be raised or the Return should be False
        :return: If _raise is True, no return; If _raise is False, return True on success and False on Exception
        :raises DatabaseRequestError: On failure to pull DB (If _raise)
        """
        # Pull DB
        try:
            q_collect, _, _ = dbi.query_collection('MC', 'queue')
        except cexc.DatabaseRequestError as dre:
            system_log.exception("Coordinator failed to get queue collection")
            tasking_log.info("Failed to pull queue collection from DB")
            if _raise:
                raise dre
            else:
                return False
        try:
            hist_collect, _, _ = dbi.query_collection('MC', 'historical')
        except cexc.DatabaseRequestError as dre:
            system_log.exception("Coordinator failed to get historical collection")
            tasking_log.info("Failed to pull historical collection from DB")
            if _raise:
                raise dre
            else:
                return False
        tasking_log.info("Queue update from DB complete")

        self.connectivity = [True, 0]

        # Create quick access lib for queue collection
        self.queue_collection = dict()
        self.queue_collection['_'] = list()
        for k, v in q_collect.items():
            q_name = v[Q_NAME]
            try:
                q_date = datetime.strptime(v[DBG_DATE_CREATED], DBQ_DATE_FORMAT)
            except ValueError:
                system_log.exception("Failed to parse queue creation date from DB collection")
                q_date = datetime.now()
            except TypeError:
                system_log.warning(f"Queue '{q_name}' missing a date created field, assuming 'now'")
                q_date = datetime.now()
            q_status = v[DBQ_STATUS]
            q_next = dbi.get_step_number(v)
            if q_next == '0':
                # This queue is actually complete
                continue
            q_agent = v[Q_OPERATIONS_LIST][q_next][QOP_AGENT]
            self.queue_collection[q_name] = v
            self.queue_collection['_'].append(Card(q_name, q_date, q_status, q_next, q_agent))

        # Create historical (contains more than just queues, so needs to be filtered)
        self.historical_collection = dict()
        for k, v in hist_collect.items():
            q_name = v.get(Q_NAME, None)
            if q_name is None:
                continue
            self.historical_collection[q_name] = v
        self.historical_collection['_'] = list(self.historical_collection.keys())

        return True

    def _check_dependency(self, parent):
        """ Checks to see if the parent is in historical and is complete

        :param parent: Name of a queue document
        :return: True - in historical AND is dbi.DBQS_DONE
        """
        if parent is None:
            return True
        elif parent == NO_DEPENDENCIES:
            return True
        elif parent in self.historical_collection['_']:
            tasking_log.info(f"{parent} found in historical")
            return self.historical_collection[parent][DBQ_STATUS] == DBQS_DONE
        else:
            tasking_log.info(f"{parent} not in historical, dependents must wait")
            return False

    def check_dependency(self, parents):
        """ Wrapper for _check_dependency

        Checks to see if all parents are in historical and are complete

        :param parents: A string or list of dependencies to check
        :return: True - if all parents are in historical and are dbi.DBQS_DONE
        (also True if parents == dbi.NO_DEPENDENCIES or None)
        """
        if parents is None:
            tasking_log.info(f"No dependencies detected")
            return True
        elif parents == NO_DEPENDENCIES:
            tasking_log.info(f"No dependencies detected")
            return True
        elif isinstance(parents, list):
            return all([self._check_dependency(parent) for parent in parents])
        elif isinstance(parents, str):
            return self._check_dependency(parents)
        else:
            raise ValueError(f"Dependency check expects a string or list not: '{parents}' ({type(parents)})")

    def check_resource(self, candidate_card: Card):
        """ Utilizes the Resource Capacity Manager to determine if a high-capacity agent can perform the task at hand

        :param candidate_card: (q_name, q_date, q_status, q_next, q_agent)
        :return: True - Good to run, False - Not good to run
        """
        try:
            queue_doc = self.queue_collection[candidate_card.q_name]
        except KeyError:
            msg = f"Queue '{candidate_card.q_name}' not loaded"
            system_log.exception(msg)
            tasking_log.info(msg)
            return False
        try:
            operations_list = queue_doc[Q_OPERATIONS_LIST]
        except KeyError:
            msg = f"Queue '{candidate_card.q_name}' missing operations list"
            system_log.exception(msg)
            tasking_log.info(msg)
            return False
        candidate_operation = operations_list.get(candidate_card.q_next, None)
        if candidate_operation is None:
            tasking_log.info(f"Candidate operation not found: {candidate_card}")
            return False
        for k, v in self.resource_ctrl.items():
            try:
                is_valid = v.manageable(candidate_operation, check_only=True)
            except ValueError:
                tasking_log.info(f"Context '{k}' poorly defined: {candidate_card}")
                # The context manager is bad, assume the operation cannot be run
                return False

            if (is_valid is None) or is_valid:
                # The context manager does not apply or operation is valid, either way, continue checking other
                # context managers
                continue
            else:
                # The context manager says this is bad, stop checking
                tasking_log.info(f"Context '{k}' failed for: {candidate_card}")
                return False
        # The context managers are not bad and, of the ones that applied, none of them took issue, this is fine
        return True

    def janitor(self):
        """ (Core) Cleans up finished tasks (step_processor calls)

        :return: None
        """
        try:
            while True:
                bad_keys = list()
                with self.threads as threads:
                    for key, value in threads.items():
                        if not value.is_alive():
                            bad_keys.append(key)
                    # delete dead thread
                    for key in bad_keys:
                        threads.pop(key, None)
                time.sleep(5)
        except Exception:  # noqa
            system_log.exception("Janitor encountered an unhandled exception")

    def processor(self):
        """ Looks for an operation to run, then spawns a thread to manage its execution

        Infinite loop until child_com is deleted

        :return: None
        """
        while True:
            # Pause if not in active mode
            while True:
                with self.child_com as cc:
                    try:
                        if cc['ui'][mcs.C_ACTIVE_MODE]:
                            break
                        else:
                            pass
                    except KeyError:
                        # cc is deleted, this is shutdown
                        return
                time.sleep(1)

            # cleaning out finished queues should no longer be necessary
            # will need to call update whenever a task finished--poor DB
            # # Have them call a common resource that forces the requests to be spaced out by 15 seconds?

            # Get next process
            try:
                next_process = self.get_next_process()  # returns a valid queue_name or raises an exception
            except cexc.DatabaseRequestError:
                # <<<
                _, n_attempts = self.connectivity
                self.connectivity = [False, n_attempts + 1]
                if n_attempts > 3:
                    new_database_settings = DatabasePrompt(
                        tk.Tk(),
                        defaults=(dbi.DATABASE_URL, dbi.DATABASE_PORT)
                    ).run()
                    dbi.update_database(*new_database_settings)
                    system_log.info(f"User has confirmed/updated the database network settings to:\n"
                                    f"\t (url, port)={(dbi.DATABASE_URL, dbi.DATABASE_PORT)}")
                    continue
                # >>> This is to halt operation until DB connectivity is restored by user
                #     If this breaks things, search for "self.connectivity" and "ui_database_prompt"
                #     and comment them out
                pass  # self.update failed, skip until we can sync with DB again
            except cexc.NoViableOptionsFound:
                pass  # Couldn't find a good next task
            except Exception:  # noqa
                system_log.exception("Scheduler encountered an Error (Setting 'next_process' to None)")
            else:
                system_log.debug(f"Selected queue {next_process}")
                with self.threads as threads:
                    threads[next_process] = Thread(target=self.tasque_manager, args=(next_process,), daemon=True)
                    threads[next_process].start()
            time.sleep(self.polling_frequency)  # Frequency with which we poll the database for queue stuff

    def get_next_process(self) -> str:
        """ Of the available processes (documents in the plate database queue collection), it returns the name
        of the process with the highest score subject to certain constraints

        :return: A name of the queue that has been selected
        """
        self.update(_raise=True)  # will raise DatabaseRequestError on failure

        # Check for overrides (empty or None --> [], str --> [str], list)
        with self.child_com as cc:
            cc_q_block: list = cc['ui'][mcs.C_QUEUE_BLOCK]

        # Valid idle queues are those which are:
        #   1) status is idle
        #   2) not in the blacklist
        #   3) have all dependencies fulfilled
        #   4) is not faulty
        #   5) is not busy
        #   6) (checking temporal conflicts is expensive, so we save this till later)
        # REF: (q_name, q_date, q_status, q_next, q_agent)
        idle_queue_reprs = [c for c in self.queue_collection['_'] if (c.q_status == DBQS_IDLE)]
        # print("Idle Queues Pre-Filter\n", [q_repr[0] for q_repr in idle_queue_reprs])
        idle_queue_reprs = [c for c in idle_queue_reprs if c.q_name not in cc_q_block]
        idle_queue_reprs = [c for c in idle_queue_reprs
                            if self.check_dependency(self.queue_collection[c.q_name].get(Q_DEPENDENCY, None))]
        # Conserving computation, silent_workers is a relative of idle_queues
        # Queues which are *waiting* (rather than truly idling) will have an idle status, but they will be both
        # status=IDLE and previous operation's IS_PAIRED = True
        silent_worker_reprs = [c for c in idle_queue_reprs if is_silent_worker(self.queue_collection[c.q_name])]
        # Back to computing idle_queues
        idle_queue_reprs = [c for c in idle_queue_reprs if not self.status.is_faulty(c.q_agent)]
        idle_queue_reprs = [c for c in idle_queue_reprs if not self.status.is_busy(c.q_agent, self.bypass)]

        idle_queue_names = [c.q_name for c in idle_queue_reprs]
        working_queue_reprs = [c for c in self.queue_collection['_'] if (c.q_status == DBQS_WORK)]

        if not idle_queue_names:
            tasking_log.info("No idle queues")
            raise cexc.NoViableOptionsFound

        # Hard mode, need to figure out what to run next
        # Build lists of:
        #   + Names (idle_queue_names)
        #   + Scores (S)
        #   + Doesn't conflict (C)
        # S*C (element-wise) will give a vector with scores (met criteria) or zeros (failed criteria)
        # max of this list is the best checklist to do next (first is fine if multiple maxes)
        # map index of max onto Names list to yield the return
        scores = [(datetime.now() - c.q_date).total_seconds()/(8.0*60*60) for c in idle_queue_reprs]
        #                                                      ^ Age in shifts
        self.clean_penalty_box()
        multipliers = [self.get_multiplier(c) if prev_step_was_not_paired(self.queue_collection[c.q_name])
                       else get_multiplier_for_is_paired(self.queue_collection[c.q_name])
                       for c in idle_queue_reprs]

        # print("Queues\n", idle_queue_names)
        # print("Queue scores\n", scores)
        # print("Queue multipliers\n", multipliers)

        no_conflict = [True] * len(idle_queue_reprs)
        # First Check for scheduling conflicts:
        for index, candidate in enumerate(idle_queue_reprs):
            try:
                if not self.check_resource(candidate):
                    raise StopIteration
                for worker in working_queue_reprs:
                    if self.detect_conflict(candidate, worker):
                        raise StopIteration
                for worker in silent_worker_reprs:
                    if self.detect_conflict(candidate, worker, silent_worker=True):
                        raise StopIteration
            except StopIteration:
                no_conflict[index] = False
        # Then Check for resource conflicts (more expensive, only do if necessary)
        if any(no_conflict):
            # Pull Location data from DB to initialize Collision search
            try:
                self._collision_register_update()
            except cexc.DatabaseRequestError as dre:
                system_log.warning("Failed to update location registers")
                tasking_log.warning("Failed to update location registers")
                raise dre
            # Perform search
            for index, candidate in enumerate(idle_queue_reprs):
                if not no_conflict[index]:
                    pass
                elif self.detect_collision(working_queue_reprs, candidate):
                    tasking_log.info(f"Resource collision found for {candidate} against workers")
                    no_conflict[index] = False
                elif self.detect_collision(silent_worker_reprs, candidate, silent_workers=True):
                    tasking_log.info(f"Resource collision found for {candidate} against silent workers")
                    no_conflict[index] = False

        filtered = list(map(lambda s, c, m: s * c * m, scores, no_conflict, multipliers))
        # tasking_log.debug("\n".join([str(i) for i in zip(idle_queue_names, filtered)]))
        if max(filtered) <= 0.0:
            tasking_log.info("All idle queues have a scheduling conflict")
            raise cexc.NoViableOptionsFound
        try:
            return idle_queue_names[filtered.index(max(filtered))]
        except ValueError:
            raise cexc.NoViableOptionsFound

    def detect_conflict(self, candidate_card: Card, worker_card: Card, silent_worker=False):
        """ Used to detect temporal conflicts

        :param candidate_card: (q_name, q_date, q_status, q_next, q_agent)
        :param worker_card: (q_name, q_date, q_status, q_next, q_agent)
        :param silent_worker: True if silent (idle due to waiting)/False otherwise
        :return: True - There is a conflict (do not schedule next)
        """
        if silent_worker:
            if candidate_card.q_name == worker_card.q_name:
                # Trying to put enforce paired wait here >>>
                # candidate = worker and worker is silent --> candidate is in its waiting period
                # Enforce the wait
                temp = dbi.get_step_number(self.queue_collection[candidate_card.q_name])
                prev_op_step = str(int(temp) - 1)
                prev_op = self.queue_collection[candidate_card.q_name][Q_OPERATIONS_LIST].get(prev_op_step, None)
                if prev_op is None:
                    print(f"How is a silent worker's previous step a NoneType? ({candidate_card.q_name})")
                    return False
                if dbi.get_pairedness(prev_op) == 'yes':
                    try:
                        prev_op_end_time = datetime.strptime(prev_op[QOP_END], mcs.TIME_FORMAT)
                    except TypeError:
                        # Mute the warning if the end time is None bc there is some delay updating the DB
                        try:
                            prev_op_start_time = datetime.strptime(prev_op[QOP_START], mcs.TIME_FORMAT)
                        except (ValueError, TypeError):
                            # If, while trying to see if we can mute the 'no end time' warning, we discover that
                            # we're missing a start time, then there is a problem
                            system_log.error(f"Start time missing/wrong by which to compare End time:\n\t{prev_op}")
                            raise RuntimeError("A pair task is missing requisite time data to function properly")
                        else:
                            time_elapsed = datetime.now() - prev_op_start_time
                            interval, wait_time = functionals.make_interval_and_wait_times(
                                dbi.get_reasonable_time(prev_op, 60, 15)
                            )
                            time_expected = timedelta(seconds=wait_time + interval)
                            if time_elapsed > time_expected:
                                system_log.warning(f"Previous end time not set:\n\t{prev_op}")
                            return True
                    except ValueError:
                        system_log.exception(f"Problem interpreting previous end time: {prev_op[QOP_END]}\n"
                                             f"\tRequesting USER to update database manually")
                        return True
                    else:
                        # prev_op_est_time = dbi.get_reasonable_time(prev_op)
                        _, prev_op_schedule_time = get_time_information(prev_op)
                        _now = datetime.now()
                        criterion = (_now - prev_op_end_time) < timedelta(seconds=prev_op_schedule_time)
                        # print(f"\t{_now} - {prev_op_end_time} "
                        #       f"= {_now - prev_op_end_time} "
                        #       f"< {timedelta(seconds=prev_op_schedule_time)} "
                        #       f"({criterion})")
                        # If less time than the time estimate has elapsed, it's not done yet,
                        # so step N conflicts with step N-1
                        if criterion:
                            return True
                        else:
                            return False
                else:
                    print(f"How is a silent worker's previous step not paired? ({candidate_card.q_name})")
                    return False
                # <<< end enforce wait code
            # print(f"Silent worker {worker_card[0]} is on step {worker_card[3]}, setting time reference back one")
            worker_card = worker_card._replace(q_next=str(int(worker_card.q_next) - 1))  # NamedTuple's _replace() method may have changed functionality between Python versions?
            # _replace() method makes a local a copy, editing the original throws everyone else off

        candidate_sequence = self._build_sequence(candidate_card)
        worker_sequence = self._build_sequence(worker_card)

        if self.debug_mode:
            w_start = datetime.strptime("11/10/2021 14:43:00", mcs.TIME_FORMAT)
            w_end = None
        else:
            w_op = self.queue_collection[worker_card.q_name][Q_OPERATIONS_LIST]
            worker_step = w_op[worker_card.q_next]
            try:
                w_start = datetime.strptime(worker_step[QOP_START], mcs.TIME_FORMAT)
            except TypeError:
                system_log.warning(f"Previous start time not set ({worker_card.q_name}):\n\t{worker_step}")
                return True
            except ValueError:
                system_log.exception(f"Problem with parsing start time of {worker_card.q_name}:\n\t{worker_step}")
                return True
            if silent_worker:
                try:
                    w_end = datetime.strptime(worker_step[QOP_END], mcs.TIME_FORMAT)
                except TypeError:
                    # Mute the warning if the end time is None bc there is some delay updating the DB
                    time_elapsed = datetime.now() - w_start
                    interval, wait_time = functionals.make_interval_and_wait_times(
                        dbi.get_reasonable_time(worker_step, 60, 15)
                    )
                    time_expected = timedelta(seconds=wait_time + interval)
                    if time_elapsed > time_expected:
                        system_log.warning(f"Previous end time not set:\n\t{worker_step}")
                    return True
                except ValueError:
                    system_log.exception(f"Problem with parsing end time of {worker_card.q_name}")
                    return True
            else:
                w_end = None

        c_time_table = self._build_timetable(datetime.now(), candidate_sequence)
        w_time_table = self._build_timetable(w_start, worker_sequence, w_end)

        if self._search_timetable(candidate_sequence, c_time_table, worker_sequence, w_time_table):
            tasking_log.info(f"Temporal conflict found between {candidate_card} "
                             f"and {worker_card} [silent worker = {silent_worker}]")
            return True
        return False

    def _build_sequence(self, initial_card: Card):
        """ Construct a sequence of paired operations from an initial step

        :param initial_card: The card of the initial step (q_name, q_date, q_status, q_next, q_agent)
        :return: A list of steps from the operations list of the associated queue
        """
        idx = str(int(initial_card.q_next))
        operations_list = self.queue_collection[initial_card.q_name][Q_OPERATIONS_LIST]
        initial_step = operations_list[idx]

        sequence = [initial_step, ]
        while True:
            idx = str(int(idx) + 1)
            if operations_list.get(idx, None) is None:
                break
            if dbi.get_pairedness(sequence[-1]) == 'yes':
                sequence.append(operations_list[idx])
            else:
                break
        return sequence

    def _collision_register_update(self):
        """ Updates registers for each location on the platform which requires resource tracking

        Attempts to pull the DB up to five times for each location to figure out how many resources occupy the specified
        locations.  Is the initialization procedure required by '_detect_collisions()'

        :return: None
        :raises DatabaseRequestError: on exhausting attempts
        """
        self.mapping = {k: list() for k in LOCATION_CAPACITIES}
        caught_exception = None

        # Initialize platform occupancy
        for location, _sub_locations in _LOCATION_CODEX.items():
            for _sub_location in _sub_locations:
                for _ in range(5):
                    try:
                        location_return, _, _ = dbi.query_location('MC', location, _sub_location)
                    except cexc.DatabaseRequestError as dre:
                        caught_exception = dre
                        time.sleep(0.1)
                        continue
                    else:
                        break
                else:
                    raise caught_exception
                if location_return['full_matches']:
                    _id = list(location_return['full_matches'].keys())[0]
                    plate_id = location_return['full_matches'][_id][DBG_CONTAINER_NAME]
                    self.mapping[location].append(plate_id)

    def detect_collision(self, worker_reprs: List[Card], candidate_card: Card, silent_workers=False):
        """ Looks for collisions of plates during transfers around the platform.

        Generates timetable for the candidate and a cumulative timetable for all workers, reduced to only transfer
        operations, and then merges them (based on start time) into a sequence of actions.  Given an initialization of
        the locations on platform, registers for each location are incremented and decremented when a transfer occurs.
        If a register goes over its maximum value, a collision is reported.

        If a queue operation is missing the requisite information to make these assessments, it was poorly formatted and
        a collision is reported regardless of the remaining time information.

        If a register were to become negative, a warning is printed to the logger (system level) and the register is
        reset to 0; however, a collision is not reported.

        :param worker_reprs: A list of worker cards [(q_name, q_date, q_status, q_next, q_agent), ...]
        :param candidate_card: The card for the candidate (q_name, q_date, q_status, q_next, q_agent)
        :param silent_workers: Bool if worker_reprs is of silent workers or not
        :return: True - Collision, False - No collisions
        """
        # Build for candidate (singleton)
        try:
            c_transfer_sequence = self._build_transfer_sequence(candidate_card)
        except (TypeError, AttributeError):
            return True
        c_transfer_sequence.sort(key=lambda x: x.start)

        # Compile all workers (multiple)
        w_transfer_sequence = list()
        for worker_card in worker_reprs:
            if silent_workers and (candidate_card.q_name == worker_card.q_name):
                continue
            if silent_workers:
                # Need a copy, editing the original throws everyone else off
                # A silent worker needs to consider the previous step
                worker_card = worker_card._replace(q_next=str(int(worker_card.q_next) - 1))
            try:
                w_transfer_sequence.extend(self._build_transfer_sequence(worker_card))
            except (TypeError, AttributeError):
                return True
        w_transfer_sequence.sort(key=lambda x: x.start)

        mapping = deepcopy(self.mapping)
        total_transfer_sequence = c_transfer_sequence + w_transfer_sequence
        total_transfer_sequence.sort(key=lambda x: x.start)
        for step in total_transfer_sequence:
            # step.container is the true name, it may not exist for all plates referenced in a queue, but it will be the
            #   value pulled from the DB during initialization since we can't look at a wellplate and know which queue
            #   it belongs to/what it's associated queue has nicknamed it.
            # step.qn_nn is f"{queue_name}_{plate_nickname}", it always exists for all plates referenced in a queue,
            #   but is not the name as plate_id.

            lookup = functionals.dictionary_recursive_search(mapping, step.container, sub_search=True)
            if lookup is None:
                lookup = functionals.dictionary_recursive_search(mapping, step.qn_nn, sub_search=True)
            if lookup is not None:
                # We have found the plate on the platform (if we haven't, it's because it's been created or is entering
                #   from a protected location and so no source decrement is necessary)
                m_key, m_index = lookup
                try:
                    del mapping[m_key][m_index]
                except (KeyError, IndexError):
                    warn_msg = f"Location reporting negative occupancy in:\n\t{total_transfer_sequence}"
                    system_log.warning(warn_msg)
                    tasking_log.warning(warn_msg)

            if step.destination in mapping:
                # We can find the destination on the platform (if we can't, it's because it's a protected location
                #   and no destination increment is necessary)
                mapping[step.destination].append(step.container if step.container else step.qn_nn)
                destination_occupancy = len(mapping[step.destination])
                destination_capacity = LOCATION_CAPACITIES[step.destination]
                if destination_occupancy > destination_capacity:
                    tasking_log.info(f"Location '{step.destination}' "
                                     f"over capacity ({destination_occupancy}/{destination_capacity}) "
                                     f"for '{step}' in :\n\t{total_transfer_sequence}")
                    return True
        return False

    def _build_transfer_sequence(self, card: Card) -> Union[None, List[Transfer]]:
        """ Converts a sequence of steps into a transfer sequence

        A transfer sequence is of the form: [(container_name, queue_underscore_nickname, destination, start_time)]
          * source is a database location [header, [details]]
          * destination is the value keyed by 'target_destination' keyword in the transfer
          * start_time is a datetime object

        :param card: A scheduling card (q_name, q_date, q_status, q_next, q_agent)
        :return: A transfer sequence, can be an empty list
        :raises AttributeError: If the transfer operation is missing a well-defined container or destination
        """
        transfer_sequence = list()
        """ List of step tuples: (container_name, queue_nickname, destination, start_time) """

        sequence = self._build_sequence(card)
        q_name = card.q_name
        time_table = self._build_timetable(datetime.now(), sequence)
        containers = self.queue_collection[q_name].get(Q_CONTAINERS, {})

        for step, time_entry in zip(sequence, time_table):
            try:
                this_operation = step[QOP_OPERATION]
                if this_operation not in TRANSFER_OPERATIONS:
                    continue
                if this_operation == "Ss.stow":
                    destination = None
                elif this_operation == "Ss.request":
                    destination = "lpx"
                else:
                    destination = step[QOP_DETAILS]['target_destination']
            except KeyError:
                msg = f"transfer operation missing a destination: {step}"
                system_log.exception(f"Scheduler found a {msg}")
                raise AttributeError(msg)
            try:
                nickname = step[QOP_CONTAINER]
                container_name = containers[nickname].get(DBG_CONTAINER_NAME, None)
            except KeyError:
                msg = f"transfer operation missing a container or container definition: {step}"
                system_log.exception(f"Scheduler found a {msg}")
                raise AttributeError(msg)
            transfer_sequence.append(Transfer(container_name, f"{q_name}_{nickname}", destination, time_entry.start))
        return transfer_sequence

    @staticmethod
    def _build_timetable(start: datetime, sequence, silent_worker=None) -> List[TimeEntry]:
        """ Creates a list of the form [TimeEntry(start: datetime, end: datetime), ...]

        :param start: The initial start time as a datetime object
        :param sequence: A list of steps
        :param silent_worker: The end time of the previous step (If the first step is a silent worker)
        :return: A list for generating the conflict matrix
        """
        time_table = [TimeEntry(start, silent_worker), ]
        for idx, step in enumerate(sequence):
            time_est, time_sch = get_time_information(step)  # (est, sch)
            if silent_worker and idx == 0:
                time_table.append(TimeEntry(time_table[idx].end + timedelta(seconds=time_sch)))
                continue
            time_table[idx].end = time_table[idx].start + timedelta(seconds=time_est)
            if time_sch is not None:
                time_table.append(TimeEntry(time_table[idx].end + timedelta(seconds=time_sch)))
        return time_table

    def _search_timetable(self, sequence_a, time_table_a: List[TimeEntry], sequence_b, time_table_b: List[TimeEntry]):
        """ Searches a pair of timetables for a conflict by constructing a conflict matrix

        :param sequence_a: A list of steps
        :param time_table_a: [TimeEntry(start: datetime, end: datetime), ...]
        :param sequence_b: A list of steps
        :param time_table_b: [TimeEntry(start: datetime, end: datetime), ...]
        :return: True - conflict detected, False - otherwise
        """
        # print([(i[0].strftime("%H:%M:%S"), i[1].strftime("%H:%M:%S")) for i in time_table_a])
        # print("vs")
        # print([(i[0].strftime("%H:%M:%S"), i[1].strftime("%H:%M:%S")) for i in time_table_b])
        # print("")

        for idx_a, step_a in enumerate(sequence_a):
            if self.bypass and (step_a[QOP_AGENT] == self.bypass):
                continue
            for idx_b, step_b in enumerate(sequence_b):
                if step_a[QOP_AGENT] == step_b[QOP_AGENT]:
                    # If a future agent in the chained set of operations is currently faulty we should not start
                    candidate_agent_is_faulty = self.status.is_faulty(subsystem=step_a[QOP_AGENT])
                    if candidate_agent_is_faulty:
                        return True

                    # if (B's step starts after A's step ends)
                    # or (A's step starts after B's step ends)
                    # then there is no conflict
                    # Otherwise, there is overlap
                    if not (
                            (time_table_a[idx_a].end < time_table_b[idx_b].start)
                            or
                            (time_table_b[idx_b].end < time_table_a[idx_a].start)
                    ):
                        # print(f"Collision {idx_a} @ {step_a[QOP_AGENT]} and {idx_b} @ {step_b[QOP_AGENT]} [ "
                        #       f"{tuple(i.strftime('%H:%M:%S') for i in time_table_a[idx_a])} & "
                        #       f"{tuple(i.strftime('%H:%M:%S') for i in time_table_b[idx_b])} ]")
                        return True
        return False

    def tasque_manager(self, next_process):
        """ Threaded processor for a task

        Runs in the background and checks for confirmation before exiting\n
        Exiting signals the task is complete and shifts self.working[key]-->self.waiting[key]

        :param next_process: Dictionary key for self.queue_collection
        :return: None
        """
        task = self.queue_collection.get(next_process, None)
        if task is None:
            return
        step_num = dbi.get_step_number(task)

        # move from waiting to working
        try:
            dbi.move(next_process, _from=DBQS_IDLE, _to=DBQS_WORK)
        except cexc.DatabaseRequestError:
            system_log.exception(f"Failed to change queue '{next_process}' from '{DBQS_IDLE}' to '{DBQS_WORK}'")
            return

        if self.debug_mode:
            print(task.get(Q_NAME, 'No Name'), step_num)
            return

        process: dict = task[Q_OPERATIONS_LIST][step_num]
        process.update({Q_NAME: next_process, DBQ_STEP_NUM: step_num})

        # Run
        try:
            self.step_processor(process)
        except cexc.UserVoidCheckpoint as uvc:
            q_move = {"_from": "Lookup", "_to": DBQS_IDLE}
            system_log.debug(f"User has voided checkpoint {uvc.get_checkpoint_key()} in {next_process}")
        except cexc.CoordinatorUserStop:
            # User said stop
            q_move = {"_from": "Lookup", "_to": DBQS_FAIL}
            system_log.exception(f"Step in {next_process} interrupted by error(s):")
        except cexc.TaskExecutionError:
            # Fatal errors and database exceptions
            q_move = {"_from": "Lookup", "_to": DBQS_FAIL}
            system_log.exception(f"Step in {next_process} interrupted by error(s):")
        except cexc.DatabaseGeneralException:
            # Connection failed
            q_move = {"_from": "Lookup", "_to": DBQS_FAIL}
            system_log.exception(f"Step in {next_process} interrupted by error(s):")
        except cexc.ConfirmationTimeout:
            q_move = {"_from": "Lookup", "_to": DBQS_FAIL}
            system_log.exception(f"Step in {next_process} timed out:")
        except Exception:  # noqa
            q_move = {"_from": "Lookup", "_to": DBQS_FAIL}
            system_log.exception(f"Coordinator.step_processor({next_process}) encountered an unhandled exception")
        else:
            # Things went normally
            q_move = {"_from": DBQS_WORK, "_to": DBQS_IDLE}
            system_log.debug(f"Step in {next_process} terminated")

        if process[QOP_OPERATION] != 'complete_queue':
            try:
                dbi.move(next_process, **q_move)
            except cexc.DatabaseGeneralException:
                system_log.exception(f"Failed to move '{next_process}' from '{q_move['_from']}' to '{q_move['_to']}'")

    def step_processor(self, task: dict):
        """ Executes an operation and confirms with the agent

        :param task: dictionary of operation with added q_name and step_num keywords
        :return: None
        :raises CoordinatorUserStop: If the user chooses to halt/abort operation.
        :raises TaskExecutionError: If user response not recognized, If a queue cannot be moved to Historical,
          If a task cannot be confirmed via DB, If the DB contradicts the function return in the event of error, If
          the DB does not corroborate a completion after a given number of attempts, If the function response was FATAL.
        :raises UserVoidCheckpoint: If a checkpoint is raised without update.
        :raises ConfirmationTimeout: If a checkpoint times out.
        :raises BadQueueFormatError(DatabaseGeneralException): If marktime fails.
        """
        with self.child_com as cc:
            if cc['ui'][mcs.C_SAFE_MODE]:
                user_resp = functionals.quick_gui(title="Safe Mode",
                                                  dialog=f"Run\n"
                                                         f"Queue: {task[Q_NAME]}\n"
                                                         f"Step: {task[DBQ_STEP_NUM]}: {task[QOP_OPERATION]}\n"
                                                         f"--OK to run--    --Halt to pause this queue--",
                                                  buttons={"Halt": functionals.void},
                                                  ret_if_ok="OK")
                if user_resp != "OK":
                    raise cexc.CoordinatorUserStop(f"User aborted operation of {task[Q_NAME]} "
                                                   f"at step {task[DBQ_STEP_NUM]}.")

        safety_check = self.safety(task)
        if safety_check is not None:
            user_resp = functionals.quick_gui(title="Safety Check Timeout",
                                              dialog=f"The following step failed due to timing out\n"
                                                     f"on the Ra-Lh safety check\n"
                                                     f"Queue: {task[Q_NAME]}\n"
                                                     f"Step: {task[DBQ_STEP_NUM]}: {task[QOP_OPERATION]}\n"
                                                     f"--OK to run--    --Abort to not--",
                                              buttons={"Abort": functionals.void},
                                              ret_if_ok="OK")
            if user_resp != "OK":
                raise cexc.CoordinatorUserStop(f"User aborted operation of {task[Q_NAME]} "
                                               f"at step {task[DBQ_STEP_NUM]}.")

        task_ref = f"({task.get(Q_NAME, 'Unknown')}, " \
                   f"#{task.get(DBQ_STEP_NUM, '_')}: " \
                   f"{task.get(QOP_OPERATION, 'Unspecified')})"

        # DB Pre-check
        try:
            if dbi.check_queue_step(task[Q_NAME], task[DBQ_STEP_NUM]) == 'yes':
                system_log.warning(f"DB has {task[Q_NAME]} #{task[DBQ_STEP_NUM]} "
                                   f"({task[QOP_OPERATION]}) marked as complete already")
                user_resp = functionals.quick_select(title="Troubleshooter",
                                                     dialog=f"Database says this task\n"
                                                            f"({task[Q_NAME]}: {task[DBQ_STEP_NUM]})\n"
                                                            f"[{task[QOP_OPERATION]}]"
                                                            f"has already been completed.\n"
                                                            f"How should the Coordinator proceed?",
                                                     options=["Run anyway", "Skip this step", "Halt"],
                                                     default="Skip this step")
                if user_resp == "Run anyway":
                    pass
                elif user_resp == "Skip this step":
                    system_log.debug(f"User elected to skip {task_ref}")
                    return
                elif user_resp == "Halt":
                    raise cexc.CoordinatorUserStop(f"User aborted operation of {task[Q_NAME]} "
                                                   f"at step {task[DBQ_STEP_NUM]}.")
                else:
                    raise cexc.TaskExecutionError("quick_select gave an unknown response")
            else:
                pass
        except cexc.DatabaseGeneralException:
            system_log.warning(f"Failed to pre-check {task[Q_NAME]} for "
                               f"completion of step #{task[DBQ_STEP_NUM]}",
                               exc_info=True)
            user_resp = functionals.quick_select(title="Troubleshooter",
                                                 dialog=f"Failed to pre-check\n"
                                                        f"{task_ref}\n"
                                                        f"How should the Coordinator proceed?",
                                                 options=["Run anyway", "Skip this step", "Halt"],
                                                 default="Run anyway")
            if user_resp == "Run anyway":
                pass
            elif user_resp == "Skip this step":
                system_log.debug(f"User elected to skip {task_ref}")
                return
            elif user_resp == "Halt":
                raise cexc.CoordinatorUserStop(f"User aborted operation of {task[Q_NAME]} "
                                               f"at step {task[DBQ_STEP_NUM]}.")
            else:
                raise cexc.TaskExecutionError("quick_select gave an unknown response")

        # RUN
        system_log.debug(f"Assigning {task[QOP_CONTAINER]} to {task[QOP_AGENT]}")
        task_key = msgm.Message.gen_idemp_key("DB")
        directive = msgm.Message.build_from_args(
            "CMD", "MC", task[QOP_AGENT], task_key, "__.Run_From_DB",
            oprtn.Operation({oprtn.FUNC: task[QOP_OPERATION],
                             oprtn.KWARGS: {Q_NAME: task[Q_NAME],
                                            oprtn.PID: task[QOP_CONTAINER],
                                            DBQ_STEP_NUM: task[DBQ_STEP_NUM]},
                             oprtn.AGENT: task[QOP_AGENT]
                             }).package(),
            msgm.P_CMD_RR)
        self.outbox.enqueue(directive)

        checkpoint = mcs.Checkpoint(completion=None,
                                    location=task[QOP_AGENT],
                                    level=None,
                                    data=task[QOP_OPERATION],
                                    queue=task[Q_NAME],
                                    step=task[DBQ_STEP_NUM],
                                    operation=task[QOP_OPERATION],
                                    task_key=task_key)
        self.status.add_checkpoint(task_key, checkpoint)

        checkpoint_update = checkpoint.wait_for_update(
            True,
            get_time_est=lambda: dbi.get_reasonable_time(task, default=60, minimum=15),
            get_inw_times=lambda x: functionals.make_interval_and_wait_times(x),
            update_time_est=lambda: max(dbi.get_time_est(checkpoint.queue, checkpoint.step, 60), 15)
        )
        system_log.info(f"{task_key} has been updated to {checkpoint.completion}")

        if checkpoint_update == -1:
            dbi.mark_time(task[Q_NAME], task[DBQ_STEP_NUM], task[QOP_OPERATION], stop=True, logger=system_log)
            raise cexc.UserVoidCheckpoint(checkpoint=task_key)
        elif checkpoint_update == 0:
            # The task timed out, but we may want to keep watching for the response
            Thread(target=dbi.silent_monitor_checkpoint,
                   args=(checkpoint,
                         self._pipes, self.status, task_key, task_ref,
                         task[Q_NAME], task[DBQ_STEP_NUM], task[QOP_OPERATION],
                         system_log),
                   daemon=True).start()
            raise cexc.ConfirmationTimeout(f"Task '{task_key}' timed out after "
                                           f"{max(dbi.get_time_est(checkpoint.queue, checkpoint.step, 60), 15)} "
                                           f"seconds")

        if checkpoint_update > 0:
            if task[QOP_OPERATION] != 'complete_queue':
                dbi.mark_time(task[Q_NAME], task[DBQ_STEP_NUM], task[QOP_OPERATION], stop=True, logger=system_log)

        if checkpoint_update == 1:
            system_log.info(f"{task_key} is completed (based on function returns), asking DB to confirm...")

            if task[QOP_OPERATION] == 'complete_queue':
                for _ in range(0, 4):
                    if task[Q_NAME] in dbi.get_historical_queue_names():
                        self.status.release_checkpoint(task_key)
                        return
                    else:
                        time.sleep(5)
                system_log.info(f'A complete_queue operation returned success '
                                f'but after exhausting the maximum number of DB access attempts, '
                                f'its queue document is still not in historical')
                raise cexc.TaskExecutionError(f"Task {task_ref} failed to move {task[Q_NAME]} to historical")

            # May need to wait for DB
            for _ in range(0, 4):
                try:
                    if dbi.check_queue_step(task[Q_NAME], task[DBQ_STEP_NUM]) == DB_YES:
                        # DB confirmed, exit
                        system_log.info(f"DB confirmed completion of {task_ref}!")
                        self.status.release_checkpoint(task_key)
                        return
                    else:
                        time.sleep(5)
                except cexc.DatabaseGeneralException:
                    # DB gave bad response code, stop trying
                    system_log.exception(f'DatabaseGeneralException in confirmation step')
                    raise cexc.TaskExecutionError(f"DB gave a bad return code")
            else:
                system_log.warning(f'Maximum number of attempts to to check {task[Q_NAME]} for '
                                   f'completion of step #{task[DBQ_STEP_NUM]} has been reached.')
                raise cexc.TaskExecutionError(f"Exhausted maximum number of DB access attempts"
                                              f" when trying to confirm step completion")
        else:
            try:
                if dbi.check_queue_step(task[Q_NAME], task[DBQ_STEP_NUM]) == 'yes':
                    raise cexc.TaskExecutionError(f"Function returns indicated {task_ref} "
                                                  f"had failed, but DB has this task marked as complete!")
            except cexc.DatabaseGeneralException:
                system_log.exception(f"Could not check {task_ref} against DB to confirm task failure")

            if checkpoint.level == mcs.V_BUSY:
                dbi.mark_time(task[Q_NAME], task[DBQ_STEP_NUM], task[QOP_OPERATION], reset=True, logger=system_log)
            if checkpoint.level in [mcs.V_BUSY, mcs.V_PROBLEM]:
                self.add_penalty(task[Q_NAME], task[QOP_OPERATION])
            # completion was false, boo, now handle the response
            functionals.handle_false_ret_obj(self._pipes, self.status, task_key, task_ref)
        return

    def add_penalty(self, q_name, q_operation):
        """ Adds a queue to the penalty box

        :param q_name: The queue being penalized
        :param q_operation: The operation which caused the penalty
        :return: None
        """
        system_log.debug(f"DEBUG for penalty box: {self.penalty_box}")
        self.penalty_box[q_name] = (datetime.now(), q_operation)

    def get_multiplier(self, task_card: Card):
        """ Used to get a multiplier for non-paired tasks.  It is intended to prevent the scheduler from soft-locking
        queues on busy/problem responses.

        Why don't we penalize based on agent rather than queue?  This allows us to not penalize an agent for a recovery
        based busy/problem response (e.g. don't penalize the Lh for a request for more tips)

        Scores
          * 1 - This queue-operation pair has not responded "busy/problem" in the last 5 minutes
          * 1e-5 - This queue-operation pair has responded "busy/problem" in the last 5 minutes
          * 0 - The operation improperly specified

        :param task_card: (q_name, q_date, q_status, q_next, q_agent)
        :return: A multiplier for the score
        """
        q_name = task_card.q_name
        q_next = task_card.q_next
        candidate_operation = self.queue_collection[q_name][Q_OPERATIONS_LIST].get(q_next, {}).get(QOP_OPERATION, None)
        if candidate_operation is None:
            tasking_log.info(f"Operation name missing for {task_card}")
            return 0
        prev = self.penalty_box.get(q_name, None)
        if prev is None:
            # There is no previous busy/problem response
            return 1
        elif prev[1] == candidate_operation:
            return 1e-5
        else:
            # The previous busy/problem response was for a different operation (the operation has changed-->there was a
            # recovery attempted, do not penalize)
            return 1

    def clean_penalty_box(self):
        """ Removes members of the penalty box who have aged out (5 minutes)

        :return: None
        """
        hit_list = []
        for q_name, v in self.penalty_box.items():
            if q_name in hit_list:
                continue
            if not v:
                hit_list.append(q_name)
            elif (datetime.now() - v[0]) > timedelta(minutes=5):
                hit_list.append(q_name)
        if hit_list:
            system_log.debug(f"DEBUG for hit_list: {hit_list}")
        for q_name in hit_list:
            self.penalty_box.pop(q_name, None)

    def safety(self, task, timeout=1):
        """ Issues a 10-second wait between Ra-Lh transfer/move commands within the last 10 seconds.

        :param task: self.waiting[next_process]
        :param timeout: Number of minutes for the safety check
        :return: Fault/None
        """
        timeout = 11.0/60.0 if timeout*60.0 < 10.0 else timeout
        timer = datetime.now()
        timer_str = timer.strftime(mcs.TIME_FORMAT)
        task_disp_name = f"[{task.get(Q_NAME, '<q_name?>')}: {timer_str}, {task[QOP_OPERATION]}]"

        # What situation are we in?
        if task[QOP_OPERATION] == "move_wellplate":
            # I am Ra and I need to check if Lh has a call to move to transfer hotel
            criterion = 'transfer_hotel'
        elif task[QOP_OPERATION] == "transfer_wellplate":
            # I am Lh (process of elimination) and I need to check if the Ra has a call to move to me
            criterion = 'liquid_handler'
        else:
            # No safety check needed
            return

        while True:
            # Loop through and cleanup
            for t, s in list(self.register):
                if datetime.now() - t > timedelta(seconds=10):
                    s_disp_name = f"[{s.get(Q_NAME, '<q_name?>')}: {t.strftime(mcs.TIME_FORMAT)}, {s[QOP_OPERATION]}]"
                    tasking_log.info(f"Removing {s_disp_name} from register")
                    self.register.remove([t, s])

            # With the clean register, check for problems
            for t, s in list(self.register):
                destination = s.get(QOP_DETAILS, {}).get('target_destination', '')
                if destination == criterion:
                    # Potential problem, wait
                    tasking_log.info(f'potential problem for {task_disp_name}, waiting 5 seconds')
                    time.sleep(5)
                    break
            else:  # No problems found (no break invoked)
                tasking_log.info(f"Adding {task_disp_name} to register")
                self.register.append([timer, task])
                return

            # Timeout for safety check
            # datetime.now() is computed whenever called, timer was instantiated outside this while loop
            if datetime.now() - timer > timedelta(minutes=timeout):
                new_fault = mcs.Fault("MC",
                                      mcs.V_PROBLEM,
                                      f"Timeout on scheduling:\n{pformat(task)}",
                                      queue=task.get(Q_NAME, None))
                system_log.error(f"Safety check has been locked for more than {timeout} minutes")
                tasking_log.warning(f"Safety check has been locked for more than {timeout} minutes on {task_disp_name}")
                return new_fault


def get_time_information(task, default_working=1, minimum_working=1):
    """ Retrieves time estimate and scheduling time information from a task.

    :param task: A step in an operations list
    :param default_working: Default duration for the completion of a task
    :param minimum_working: Minimum duration for the completion of any task
    :return: (Time estimate for the task, Time requirement for paired-scheduling)
    """
    working_time = max(task.get(QOP_TIME, default_working), minimum_working)
    is_paired = task.get(QOP_DETAILS, dict()).get(QDET_PAIRED, 'no')
    if is_paired == 'no' or is_paired is None:
        scheduled_time = None
    else:
        scheduled_time = task.get(QOP_DETAILS, dict()).get(QDET_TIME, 0)
    return working_time, scheduled_time


def get_prev_step(q_doc: dict) -> dict:
    """ Given a queue document, it returns the last completed step

    :param q_doc: A queue document
    :return: The last completed step or an empty dictionary
    """
    prev_step_index = str(int(dbi.get_step_number(q_doc)) - 1)
    operations_list = q_doc.get(dbi.Q_OPERATIONS_LIST, dict())
    return operations_list.get(prev_step_index, dict())


def prev_step_was_not_paired(q_doc: dict):
    """ Previous step was not paired

    :param q_doc: A queue document
    :return: if the is_paired flag on the previous step (see: get_prev_step()) is 'no'
    """
    return dbi.get_pairedness(get_prev_step(q_doc)) == 'no'


def get_multiplier_for_is_paired(queue_doc: dict) -> float:
    """ Gets a score multiplier from a queue document based on the current step being overdue via is_paired

    Assumes it is called if the document's previous step was paired

    :param queue_doc: A queue document
    :return: The total seconds of lateness for this queue_doc or 3
    """
    try:
        previous_step = get_prev_step(queue_doc)
        prev_end_time = datetime.strptime(previous_step[dbi.QOP_END], mcs.TIME_FORMAT)
        scheduled_time = timedelta(seconds=int(previous_step.get(QOP_DETAILS, dict()).get(QDET_TIME, 0)))
    except (ValueError, TypeError, KeyError):
        return 3
    else:
        tardiness = (datetime.now() - (prev_end_time + scheduled_time)).total_seconds()
        return max(tardiness + 3, 3)  # Would it be better to have this to use relative tardiness?


def is_silent_worker(queue_doc):
    """ Helper method, is effectively "not prev_step_was_not_paired(queue_document)"

    Alt name: was_prev_step_paired

    :param queue_doc: A queue document
    :return: not prev_step_was_not_paired(queue_doc)
    """
    return not prev_step_was_not_paired(queue_doc)


if __name__ == '__main__':
    from custom_classes import Nop
    from thread_safe_data import ThreadSafeDataContainer
    import mcn_queues

    test_pipes = {
        MSG_Q_IN: Nop(),
        MSG_Q_OUT: Nop(),
        CHILD_COM: ThreadSafeDataContainer(dict()),
    }
    with test_pipes[CHILD_COM] as _cc:
        _cc['ui'] = {mcs.C_ACTIVE_MODE: True, mcs.C_QUEUE_BLOCK: list(), mcs.C_SAFE_MODE: True,
                     mcs.C_EXPIRED: [False, ], mcs.C_TASK_QUEUE: mcn_queues.LogQueue(4)}

    test_status = mcs.Status('MC')

    test = Coordinator(test_pipes, test_status, True)

    # test.run()
