""" Podalirius

Library of analysis options for the MCN

Statistics Reported:
  * Utilization (Optimal Time, Actual Time, Max Reasonable Time, Estimated Downtime)
  * Number of Queues (Complete, Idle, In-Progress, Error, NR)
  * Number of Queues with at least one Fault (Problem, Fatal, Total)
  * Number of Steps (By Queue, By Location, and Total)
  * Number of Wells (By Queue, By Location, and Total)
  * Number of faults (Total, Problem, Fatal, NR) (By Queue, By Location, and Total)

Metric Definitions:
  * Utilization is a comparison of actual time to the best and worst case times, evaluates how well the
  scheduler (a greedy, blind horizon algorithm) performs
  * Performance is a success rate (1 - error rate) which can be evaluated at the queue or instrument level and can
  be calculated on the basis or steps or wells
  * Reliability is a queue-wise metric which evaluates the ultimate success of the queue's intent

(Podalirius, brother of Aceso and associated with diagnostics)

Author: Ben C
"""

import time
from copy import deepcopy
from datetime import timedelta, datetime
from os import path, walk
from pprint import pformat
from typing import Tuple, Union

import yaml

import custom_exceptions as cexc
import database_interface as dbi
import mcn_status as mcs
from custom_classes import pathfinder
from mcn_logging_manager import UP, DOWN
from metric_logging_manager import metric_log

try:
    import collections
    from ortools.sat.python import cp_model
except ModuleNotFoundError:
    print("ortools not installed--ignore if not MC")
    pass

try:
    path_to_logs = pathfinder(r"Logs", "mcn_logfile.log")
except:  # noqa
    metric_log.info("Failed to load the path of the MCN logs")
    path_to_logs = None

# Typing Hint
Timeframe = Tuple[Union[int, float, datetime, timedelta, None], ...]

CP_MODEL_DECODER = {
    0: "Unknown",
    1: "Invalid",
    2: "Feasible",
    3: "Infeasible",
    4: "Optimal"
}
""" Maps CP Model codes (0-4) with human readable forms """

VERB = {
    'n_q_with_fatal':      "(1) # queues with fatal errors",
    'n_q_with_problem':    "(1) # queues with problem errors",
    'n_q_with_fault':      "(1) # queues with faults",
    'n_steps_total':       "(3) Total # steps covered by all queues",
    'n_wells_total':       "(3) Total # wells covered by all queues",
    'n_complete_q':        "(2) # completed queues",
    'n_idle_q':            "(2) # idle queues",
    'n_in-progress_q':     "(2) # queues in progress",
    'n_error_q':           "(2) # queues in error",
    'n_not-reported_q':    "(2) # queues with unrecognized status",
    'n_fatal_faults':      "(4) Total # of fatal faults",
    'n_problem_faults':    "(4) Total # of problem faults",
    'n_total_faults':      "(4) Total # of faults",
    'optimal_time':        "(5) Optimal time [hr]",
    'max_reasonable_time': "(5) Maximum reasonable time [hr]",
    'actual_time':         "(5) Actual time taken [hr]"
}
""" Maps keyword argument names to human readable forms """

WIDTH = 255
LOG_TIME_FORMAT = "%Y-%m-%d %H:%M:%S"  # e.g. "2022-07-06 12:56:39"
BYPASS_NUM = max(v for _, v in mcs.SYS_ENUM.items())


def _fetch_data_from_file(filepath):
    """
    Invokes yaml.safe_load to pull a file into dict format

    Re-keys the collection from database _id to queue_name at the toplevel

    :param filepath: A path to a yaml file
    :return: {queue name: queue document}
    """
    with open(filepath, "r") as fh:
        doc = yaml.safe_load(fh)
    return {v.get(dbi.Q_NAME, "?"): v for _, v in doc.items()}


def pad(k, src):
    """
    Helper method for printing

    :param k: A string to print
    :param src: The collection the string belongs to (Iterable)
    :return: The required ' ' padding required for all elements of the collection to line up
    """
    longest = max(len(_key) for _key in src)
    delta = longest - len(k)
    return ' ' * delta


class MetricMonitor:
    """
    Calculates and reports platform metrics
    """
    def __init__(self, name='MC', trimming=1, calc_util=False):
        """
        Constructor for the MetricMonitor class

        :param name: A name for Database credentials
        :param trimming: Trim queue documents to remove incomplete steps: 0 trim nothing, 1 trim active, 2 trim both
        :param calc_util: Flag to calculate utilization (expensive)
        """
        self.name = name
        """ [reliability, availability, performance, performance (by system)] """

        self.trimming = trimming
        self.calc_util = calc_util
        self._time_mode_selector = {
            'first': self._get_first_timestamp,
            'last': self._get_last_timestamp,
            'created': self._get_creation_timestamp,
        }
        self.queue_collection = {}
        self.logs = path_to_logs

    # Time stuff # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    @staticmethod
    def _get_timestamps(queue_doc: dict):
        """
        Converts a queue document into a list of datetime objects pulled form start queue operation state/end times

        :param queue_doc: A queue document
        :return: A list of all start and end times (not ordered, not start-end paired)
        """
        times = list()
        for _, operation in queue_doc.get(dbi.Q_OPERATIONS_LIST, dict()).items():
            try:
                times.append(datetime.strptime(operation[dbi.QOP_START], mcs.TIME_FORMAT))
            except (ValueError, TypeError):
                pass
            try:
                times.append(datetime.strptime(operation[dbi.QOP_END], mcs.TIME_FORMAT))
            except (ValueError, TypeError):
                pass
        return times

    def _get_last_timestamp(self, queue_doc: dict):
        """
        Identifies the latest time stamp recorded in a queue document's operations list (or datetime.min)

        :param queue_doc: A queue document
        :return: The last timestamp as a datetime object (or datetime.min if no timestamps recorded)
        """
        try:
            return max(self._get_timestamps(queue_doc))
        except ValueError:
            return datetime.min  # min bc "now - this.return < cutoff" used as criteria for inclusion in calculation

    def _get_first_timestamp(self, queue_doc: dict):
        """
        Identifies the earliest time stamp recorded in a queue document's operations list (or datetime.min)

        :param queue_doc: A queue document
        :return: The first timestamp as a datetime object (or datetime.min if no timestamps recorded)
        """
        try:
            return min(self._get_timestamps(queue_doc))
        except ValueError:
            return datetime.min  # min bc "now - this.return < cutoff" used as criteria for inclusion in calculation

    def _get_bookends(self):
        """
        Looks through all queues loaded and pulls the earliest and latest timestamps from the operations list

        :return: (min, max)
        """
        starts = [self._get_first_timestamp(q_doc) for _, q_doc in self.queue_collection.items()]
        ends = [self._get_last_timestamp(q_doc) for _, q_doc in self.queue_collection.items()]
        all_times = starts + ends
        return min(all_times), max(all_times)

    @staticmethod
    def _get_creation_timestamp(queue_doc: dict):
        """
        Pulls the 'date_created' field from a queue document

        :param queue_doc: A queue document
        :return: A datetime representation of the date created or datetime.min if not present/properly formatted
        """
        try:
            return datetime.strptime(queue_doc[dbi.DBG_DATE_CREATED], dbi.DBQ_DATE_FORMAT)
        except (KeyError, ValueError, TypeError):
            return datetime.min  # min bc "now - this.return < cutoff" used as criteria for inclusion in calculation

    @staticmethod
    def _queue_to_interval_list(queue_doc: dict):
        """
        Converts the operation list of a queue document into a list of interval objects

        Interval: (step number, agent, duration, schedule time before next step, soft_tail_time)

        :param queue_doc: A queue document
        :return: A list of intervals sorted by step number
        """
        master = list()
        operations_list = queue_doc.get(dbi.Q_OPERATIONS_LIST, dict())
        # name_of_queue = queue_doc.get(dbi.Q_NAME, "nr")
        for step_num, operation in operations_list.items():
            if operation.get(dbi.QOP_OPERATION, None) == "soft_wait":
                # soft_waits are counted by appending tail times to the previous step
                continue
            agent = operation.get(dbi.QOP_AGENT, None)
            if agent is None:
                continue
            step_num = int(step_num)
            try:
                start = datetime.strptime(operation[dbi.QOP_START], mcs.TIME_FORMAT)
            except ValueError:
                continue
            except TypeError:
                continue
            try:
                end = datetime.strptime(operation[dbi.QOP_END], mcs.TIME_FORMAT)
            except ValueError:
                continue
            except TypeError:
                continue
            if dbi.get_pairedness(operation) == dbi.DB_YES:
                tail_time = operation.get(dbi.QOP_DETAILS, dict()).get(dbi.QDET_TIME, 0)
            else:
                tail_time = -1
            # # Paired actions have the trait that they are often of the form:
            # # 1: start thing (start: 0:00, end: 0:03)
            # # 2: stop thing  (start: 5:00, end: 5:02)
            # # So the true end time is 5:00 and not 0:03
            # if dbi.get_pairedness(operation) == dbi.DB_YES:
            #     next_step = operations_list.get(str(step_num + 1), None)
            #     if next_step:
            #         try:
            #             _next_start = datetime.strptime(next_step[dbi.QOP_START], mcs.TIME_FORMAT)
            #         except ValueError:
            #             _next_start = end
            #         end = _next_start
            duration = end - start
            if duration < timedelta(0):
                continue

            next_step = operations_list.get(str(step_num + 1), None)
            soft_tail_time = -1
            if next_step and isinstance(next_step, dict):
                if next_step.get(dbi.QOP_OPERATION, None) == "soft_wait":
                    soft_tail_time = next_step.get(dbi.QOP_DETAILS, dict()).get(dbi.QDET_TIME, -1)
                    # check the next step to ensure that this soft_wait isn't part of a larger holding pattern
                    next_next_step = operations_list.get(str(step_num + 2), None)
                    if next_next_step and isinstance(next_next_step, dict):
                        if next_next_step.get(dbi.QOP_OPERATION, None) == "detect_liquid_level":
                            hist = next_next_step.get(dbi.QOP_DETAILS, {}).get('volume_history', list())
                            n_reps = len(hist)
                            soft_tail_time *= n_reps  # (See note below*)
                        # elif ... (if others get added)
                        # *Note:  It would be more accurate to extract the times from volume_history and the start/end
                        #         times from detect_liquid_level to systematically reconstruct the steps of
                        #         iteratively soft_wait<-->detect_liquid_level.  However, this would require either
                        #         offsetting the step_num field or sub-indexing them (4a, 4b, 4c, etc.) which may
                        #         break things.  This is more accurate as it reflects the fact that Lh is used in those
                        #         iterations; however, it is used for ~33 seconds and soft_waits are typically >1 hour.
                        #         So this approximation is used. [future: the indexing must abide sorted()]
            # print(name_of_queue, ":", (step_num, agent, duration, tail_time, soft_tail_time))
            master.append((step_num, agent, duration, tail_time, soft_tail_time))
        master.sort(key=lambda x: x[0])

        return master

    def _queue_dependency_list(self, schedule_intervals):
        """ Used to find the dependencies of the jobs in a scheduling interval

        :param schedule_intervals: Enumerable of dict objects with queue names and queue dependencies keyed
        :return: {queue name: [dependency names, ...], ...}
        """
        dependency_lookup = {}
        for q_name, q_doc in self.queue_collection.items():
            dependencies = q_doc.get(dbi.Q_DEPENDENCY, None)
            if not dependencies:
                dependencies = []
            elif dependencies == dbi.NO_DEPENDENCIES:
                dependencies = []
            elif isinstance(dependencies, str):
                dependencies = [dependencies, ]
            if not isinstance(dependencies, list):
                print(f"Bad dependency list for {q_name}: {q_doc.get(dbi.Q_DEPENDENCY, None)}")
                continue
            dependency_lookup.update({q_name: dependencies})
        mapping = {k: i for i, k in enumerate(schedule_intervals)}
        result = {mapping.get(q_name, -1): [mapping.get(_d, -1) for _d in _dependencies]
                  for q_name, _dependencies in dependency_lookup.items()}
        result.pop(-1, None)
        return result

    def _extract_operational_time(self):
        """
        Examines log files for Online/Offline signals and compiles the amount of Offline time

        :return: Amount of offline time in hours or a string with error details
        """

        first_time, last_time = self._get_bookends()
        time_log = []
        rest_time = timedelta(seconds=0)
        log_collection = list(self.get_logging_options())
        if not log_collection:
            return "Empty log collection"
        for a_log_file in self.get_logging_options():
            with open(a_log_file, "r") as fh:
                for line in fh.readlines():
                    try:
                        timestamp = datetime.strptime(line[:19], LOG_TIME_FORMAT)
                    except ValueError:
                        return f"Error pulling up/down times from {a_log_file}"
                    if UP in line:
                        time_log.append((timestamp, True))
                    elif DOWN in line:
                        time_log.append((timestamp, False))

        time_log.sort(key=lambda x: x[0])
        if not time_log:
            # If no markers, assume all time is Uptime
            time_log.insert(0, (datetime.min, True))
            time_log.insert(-1, (datetime.max, False))
        else:
            time_log.insert(0, (datetime.min, not time_log[0][1]))
            time_log.insert(-1, (datetime.max, not time_log[-1][1]))

        # Safely clear all entries after last_time
        i = 0
        while True:
            try:
                curr = time_log[i]
            except IndexError:
                break
            if curr[0] > last_time:
                reg = curr[1]
                del time_log[i:]
                time_log.append((last_time, reg))
                break
            i += 1
        # Safely clear all entries before first_time
        i = -1
        while True:
            try:
                curr = time_log[i]
            except IndexError:
                break
            if curr[0] < first_time:
                reg = curr[1]
                del time_log[:i]
                time_log.append((first_time, reg))
                break
            i -= 1
        if not time_log:
            return "Logs contain no up/down time markers after cleaning"

        # Must alternate between Up and Down
        checker = None
        for item in time_log:
            if item[1] == checker:
                return f"Repeated {checker} marker in sorted logs"
            checker = item[1]

        # Pull DOWN intervals
        for i in range(len(time_log)-1):
            _c = time_log[i]
            if not _c[1]:
                rest_time += time_log[i+1][0] - _c[0]
        return -rest_time.total_seconds()/(60*60)

    # Queue Setters # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    def set_queue_to_database(self, include_active=False):
        """
        Populates the MetricMonitor's assembly of queue documents using the database

        :param include_active: True if incomplete queues should be included, and
        False if only completed (historical) queues should be included
        :return: A dictionary of queue documents keys by queue_name
        """
        hist_queues = dict()
        actv_queues = dict()

        # pull historical queues
        try:
            hist_collect, _, _ = dbi.query_collection('MC', dbi.COL_HISTORICAL)
        except cexc.DatabaseRequestError:
            metric_log.exception("MetricMonitor failed to get historical collection")
            hist_collect = {}
        for k, v in hist_collect.items():
            q_name = v.get(dbi.Q_NAME, None)
            if q_name is None:
                continue
            hist_queues[q_name] = v

        # pull active queues, if requested
        if include_active:
            try:
                queue_collect, _, _ = dbi.query_collection('MC', dbi.COL_QUEUE)
            except cexc.DatabaseRequestError:
                metric_log.exception("MetricMonitor failed to get queue collection")
                queue_collect = {}

            for k, v in queue_collect.items():
                q_name = v.get(dbi.Q_NAME, None)
                if q_name is None:
                    continue
                actv_queues[q_name] = v

        # trimming: 0 trim nothing, 1 trim active, 2 trim both
        hist = hist_queues if self.trimming < 2 else {k: self._trim_queue_doc(v) for k, v in hist_queues.items()}
        actv = actv_queues if self.trimming < 1 else {k: self._trim_queue_doc(v) for k, v in actv_queues.items()}
        self.queue_collection = {**actv, **hist}

    def set_queue_to_file(self, file_path):
        """
        Populates the MetricMonitor's assembly of queue documents using a file

        :param file_path: The path to a yaml file of queue documents (all assumed to be completed)
        :return: A dictionary of queue documents keys by queue_name
        """
        file_queues = _fetch_data_from_file(file_path)
        self.queue_collection = file_queues \
            if self.trimming < 2 \
            else {k: self._trim_queue_doc(v) for k, v in file_queues.items()}

    def set_log_path_to_default(self):
        """
        Have the MetricMonitor attempt to locate the log files for on/off-line data extraction

        :return:
        """
        self.logs = path.dirname(path_to_logs)

    def set_log_path_to_file(self, filename):
        """
        Have the MetricMonitor use the given file for on/off-line data extraction

        :param filename: A log file path
        :return:
        """
        self.logs = filename

    def set_log_path_to_folder(self, directory):
        """
        Have the MetricMonitor use the given folder for on/off-line data extraction

        :param directory: A path to a collection of log files
        :return:
        """
        self.logs = path.dirname(directory) if ".log" in directory else directory

    def get_logging_options(self):
        """
        :return: Yields paths to log files, files located according to the 'set_log_path_to_...' methods
        """
        if ".log" in self.logs:
            yield self.logs
            return
        for _root, _, _files in walk(self.logs):
            for name in _files:
                if (".log" in name) and ("mcn_logfile" in name):
                    yield path.join(_root, name)

    def ignore_user_errors(self):
        """ Scrubs fault record of faults with a location at "_U" (user)

        :return: None (edits dictionary in-place)
        """
        for _, q_doc in self.queue_collection.items():
            q_doc[dbi.Q_RECORD] = [fault
                                   for fault in q_doc.get(dbi.Q_RECORD, [])
                                   if fault.get('location', 'NR') != "_U"]

    @staticmethod
    def _trim_queue_doc(queue_doc: dict):
        """
        Removes steps which have not been completed from a queue document

        :param queue_doc: A queue document
        :return: A version of the queue document with incomplete steps removed
        """
        if not queue_doc:
            return queue_doc

        operations_list = queue_doc.get(dbi.Q_OPERATIONS_LIST, {})
        record = [int(ki) for ki, operation in operations_list.items()
                  if operation.get(dbi.QOP_COMPLETED, dbi.DB_NO) == dbi.DB_NO]

        if record:
            critical_step = min(record) - 1
        else:
            critical_step = max([int(k) for k in operations_list.keys()], default=99999)

        queue_doc[dbi.Q_OPERATIONS_LIST] = {k: v for k, v in operations_list.items() if int(k) <= critical_step}
        return queue_doc

    def down_select_queues_by_name(self, selected_queues: list, deselected_queues: list = None):
        """
        Reduces the internal library of queue documents to those specified

        :param selected_queues: A list of permissible queue document names
        :param deselected_queues: A list of non-permissible queue document names
        :return:
        """

        self.queue_collection = {k: self.queue_collection[k]
                                 for k in selected_queues
                                 if (k in self.queue_collection) and (k not in deselected_queues)}

    def down_select_queues_by_campaign(self, campaign_root: Union[list, str]):
        """
        Reduces the internal library of queue documents to those from a campaign

        :param campaign_root: Queues which contain this string are included
        :return:
        """
        if isinstance(campaign_root, str):
            campaign_root = [campaign_root, ]
        self.queue_collection = {k: v for k, v in self.queue_collection.items()
                                 if any(_root in k for _root in campaign_root)}

    def down_select_queues_by_time(self, timeframe: Timeframe = (7, None), time_mode='first'):
        """
        Reduces the internal library of queue documents to those which fall within a temporal window

        :param timeframe: (Oldest age, Youngest age), ages given in days, A NoneType age will remove that bound
        :param time_mode: The basis for determining an age: 'first' or 'last' (operation), or 'created' (queue)
        :return:
        """
        now = datetime.now()
        max_age, min_age = self.validate_ages(now, timeframe, time_mode)

        def determine_time_requirement(sel_queue) -> bool:
            """ Does this queue meet the time requirements?

            :param sel_queue: A queue document
            :return: Within the time window (True), else False
            """
            then = self._time_mode_selector[time_mode](sel_queue)
            if (max_age is not None) and (now - then > max_age):
                return False
            if (min_age is not None) and (now - then < min_age):
                return False
            return True

        self.queue_collection = {
            k: v for k, v in self.queue_collection if determine_time_requirement(v)
        }

    def validate_ages(self, now: datetime, timeframe: Timeframe = (7, None), time_mode='first'):
        """
        Validator to convert a variety of temporal inputs into the ages required by other methods

        :param now: A datetime object for what constitutes "now"
        :param timeframe: A tuple for (Oldest, Youngest): Int (age in days), Timedelta (age),
        Datetime (age determined by comparison with 'now'), None (no bound)
        :param time_mode: The basis for determining an age: 'first' or 'last' (operation), or 'created' (queue)
        :return:
        :raises ValueError: If the time_mode is not recognized, if timeframe mixes types (other than None),
        if an element of timeframe is not int, datetime, timedelta, or None
        """
        if time_mode not in self._time_mode_selector:
            raise ValueError(f"time_mode must be in {list(self._time_mode_selector.keys())} not '{time_mode}'")

        if any(isinstance(v, int) for v in timeframe) and any(isinstance(v, datetime) for v in timeframe):
            raise ValueError("timeframe must be of one type")
        if any(isinstance(v, int) for v in timeframe):
            max_age = None if timeframe[0] is None else timedelta(days=timeframe[0])
            min_age = None if timeframe[1] is None else timedelta(days=timeframe[1])
        elif any(isinstance(v, datetime) for v in timeframe):
            max_age = None if timeframe[0] is None else now - timeframe[0]
            min_age = None if timeframe[1] is None else now - timeframe[1]
        elif any(isinstance(v, timedelta) for v in timeframe):
            max_age, min_age = timeframe
        elif all(v is None for v in timeframe):
            max_age, min_age = timeframe
        else:
            raise ValueError("timeframe must be composed of int, datetime, timedelta, or None types")

        return max_age, min_age

    # Metrics # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    def _get_utilization(self):
        """
        Performs calculations to obtain utilization metrics

        Uses Google OrTools.CpModel to determine the optimal schedule for the queue operations.
          Operations are assigned to machines and constrained by their real duration, their scheduled 'tail' times,
          precedence ordering, and a 15-second buffer between consecutive operations.

        :return: optimal_time (s), max_reasonable_time (s), actual_time (s),
        ret_val (list of [solver, jobs_data, assigned_task_type, all_tasks, all_machines])
        :raises ValueError: If failed to make a makespan, solver status is not good,
        """
        queue_collection = {k: v for k, v in self.queue_collection.items() if v.get(dbi.Q_OPERATIONS_LIST, {})}
        if not queue_collection:
            return None
        max_reasonable_time, optimal_time, actual_time, ret_val = None, None, None, list()

        # Maximum (reasonable) makespan
        # k: (step_num, agent, duration, tailing_time, soft_tail_time) -> (int, str, timedelta, int, int)
        schedule_intervals = {k: self._queue_to_interval_list(v)
                              for k, v in queue_collection.items()}

        all_intervals = list()
        for k, intervals in schedule_intervals.items():
            for _step in intervals:
                tail = timedelta(seconds=(_step[3] if _step[3] > 0 else 0) + (_step[4] if _step[4] > 0 else 0))
                all_intervals.append(_step[2] + tail)
        max_reasonable_time = sum(all_intervals, timedelta(seconds=15 * len(all_intervals)))

        # Actual makespan
        starts = [self._get_first_timestamp(q_doc) for _, q_doc in queue_collection.items()]
        ends = [self._get_last_timestamp(q_doc) for _, q_doc in queue_collection.items()]
        try:
            actual_time = max(ends) - min(starts)
        except ValueError as ve:
            raise ValueError("Failed to get actual makespan") from ve

        # Theoretical minimum makespan
        # [[(agent #, duration, tail, soft), ...], ...]
        # grouped by queue_name then sorted by step number
        jobs_data = [
            [(mcs.SYS_ENUM.get(t[1], BYPASS_NUM), t[2].total_seconds(), t[3], t[4])
             for t in sorted(v, key=lambda x: x[0])]
            for _, v in schedule_intervals.items()
        ]
        machines_count = 1 + max(task[0] for job in jobs_data for task in job)
        all_machines = range(machines_count)
        horizon = sum(task[1] + abs(task[2]) + 15 for job in jobs_data for task in job)
        horizon = int(horizon + 0.5)
        model = cp_model.CpModel()
        task_type = collections.namedtuple('task_type', 'start end interval')
        assigned_task_type = collections.namedtuple('assigned_task_type', 'start job index duration')
        all_tasks = dict()
        machine_to_intervals = collections.defaultdict(list)

        for job_id, job in enumerate(jobs_data):
            for task_id, task, in enumerate(job):
                machine = task[0]
                duration = task[1]
                duration = int(duration + 0.5)
                suffix = '_%i_%i' % (job_id, task_id)
                start_var = model.NewIntVar(0, horizon, 'start' + suffix)
                end_var = model.NewIntVar(0, horizon, 'end' + suffix)
                interval_var = model.NewIntervalVar(start_var, duration, end_var, 'interval' + suffix)
                all_tasks[job_id, task_id] = task_type(start=start_var, end=end_var, interval=interval_var)
                machine_to_intervals[machine].append(interval_var)
        for machine in all_machines:
            model.AddNoOverlap(machine_to_intervals[machine])
        for job_id, job in enumerate(jobs_data):
            # job_id --> index # of queue name
            # job --> a sorted list of operations
            for task_id, task in enumerate(job[:-1]):
                # print(job_id, task)
                if task[3] >= 0:
                    # There was pairedness but also soft_waits
                    delay = int(max(task[2], 0) + task[3])
                    model.Add(all_tasks[job_id, task_id + 1].start >= all_tasks[job_id, task_id].end + delay)
                elif task[2] >= 0:  # and task[3] < 0
                    # If there was pairedness things going on, but no soft_waits, restrict the time
                    delay = int(task[2])
                    model.Add(all_tasks[job_id, task_id + 1].start >= all_tasks[job_id, task_id].end + delay)
                    model.Add(all_tasks[job_id, task_id + 1].start <= all_tasks[job_id, task_id].end + delay + 15)
                else:  # task[2] < 0 and task[3] < 0
                    # If no pairedness things or soft_wait things going on, just enforce precedence ordering
                    model.Add(all_tasks[job_id, task_id + 1].start >= all_tasks[job_id, task_id].end)
            # for task_id in range(len(job) - 1):
            #     model.Add(all_tasks[job_id, task_id + 1].start >= all_tasks[job_id, task_id].end)
        # Enforce the dependencies ordering
        dependencies = self._queue_dependency_list(schedule_intervals)
        for job_id, _ in enumerate(jobs_data):
            for dependency in dependencies[job_id]:
                if dependency < 0:  # Failed to lookup dependent (ignore dependency)
                    continue
                dependant = jobs_data[dependency]
                last = len(dependant) - 1
                model.Add(all_tasks[job_id, 0].start >= all_tasks[dependency, last].end)

        obj_var = model.NewIntVar(0, horizon, 'makespan')
        model.AddMaxEquality(obj_var, [all_tasks[job_id, len(job) - 1].end for job_id, job in enumerate(jobs_data)])
        model.Minimize(obj_var)
        solver = cp_model.CpSolver()
        status = solver.Solve(model)
        ret_val.append([solver, jobs_data, assigned_task_type, all_tasks, all_machines])
        if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
            optimal_time = solver.ObjectiveValue()
            optimal_time = timedelta(seconds=optimal_time)
        else:
            raise ValueError(f"Warning: status != optimal or feasible: "
                             f"{CP_MODEL_DECODER[status]} & {solver.ObjectiveValue()}")

        return optimal_time.total_seconds(), max_reasonable_time.total_seconds(), actual_time.total_seconds(), ret_val

    def get_utilization(self):
        """
        Public utilization calculation method, formats results

        :return: optimal time (hr), max reasonable time (hr), actual time (hr), down time (hr), and scheduling
        diagnostics; with times being -1 if NA
        :raises: From _get_utilization
        """
        optimal_time, max_reasonable_time, actual_time, ret_val = self._get_utilization()
        rest_time = self._extract_operational_time()

        sec_to_hour = 60*60
        optimal_time = -1 if optimal_time is None else round(optimal_time/sec_to_hour, 1)
        max_reasonable_time = -1 if max_reasonable_time is None else round(max_reasonable_time/sec_to_hour, 1)
        actual_time = -1 if actual_time is None else round(actual_time/sec_to_hour, 1)
        return optimal_time, max_reasonable_time, actual_time, rest_time, ret_val

    def get_utilization_report(self, source=None, _frmt=lambda x: x):
        """
        Used to generate a report of the utilization of the loose form

        "min: #, max: #, actual: # (#downtime) [#utilization_metric; #adjusted_utilization_metric]"

        * min - The minimum (optimal) time required to perform the operations
        * max - The maximum (reasonable) time required to perform the operations (no parallelization)
        * actual - How long the execution actually took
        * downtime - How long was DB control turned off
        * utilization_metric - (max - actual)/(max - min) in %s or "NA" if max ~= min or min/max/actual could not
        be calculated
        * adjusted_utilization_metric - (max - (actual - downtime))/(max - min) in %s or "NA" if max ~= min or
        min/max/actual could not be calculated

        :param source: A dictionary of time data to avoid recalculation on subsequent calls
        :param _frmt: Called on the report string before being returned
        :return: The report string
        """
        if source is None:
            try:
                optimal_time, max_reasonable_time, actual_time, rest_time, *_ = self.get_utilization()
            except ValueError as ve:
                return _frmt(f"Failed to calculate utilization:\n{repr(ve)}")
        else:
            optimal_time = source['optimal_time']
            max_reasonable_time = source['max_reasonable_time']
            actual_time = source['actual_time']
            rest_time = source.get('rest_time', "down time not reported")

        def safe_calc(_min, _max, val):
            if any(v <= 0 for v in [_min, _max, val]):
                return "NA"
            if 200*abs(_max - _min)/(_max + _min) < 5:
                return "NA"
            return f"{100*round((_max - val)/(_max - _min), 4)}%"

        util_metric = safe_calc(optimal_time, max_reasonable_time, actual_time)
        try:
            adj_util_metric = safe_calc(optimal_time, max_reasonable_time, actual_time - rest_time)
        except TypeError:
            adj_util_metric = "NA"

        try:
            print("*** Down time estimate: ", actual_time - (optimal_time + rest_time), "***")
        except Exception as e:
            print("*** Could not estimate down time:", repr(e), "***")

        return _frmt(f"min: {optimal_time}, max: {max_reasonable_time}; actual: {actual_time} ({rest_time})"
                     f"[{util_metric}; {adj_util_metric}]")

    def _queue_header(self):
        """
        Builds headers for the queues with their name, date created, and dates operational

        :return: {queue_name: padding + "[created ##/##/####]  (run: ##/##/#### - ##/##/####)"}
        """
        def grab(queue_doc):
            try:
                temp = self._get_timestamps(queue_doc)
                _min = min(temp).strftime("%m/%d")
                _max = max(temp).strftime("%m/%d")
            except ValueError:
                _min = "NR"
                _max = "NR"
            return self._get_creation_timestamp(queue_doc).strftime("%m/%d"), _min, _max
        return {k: pad(k, self.queue_collection) + "[created {}]  (run: {} - {})".format(*grab(v))
                for k, v in self.queue_collection.items()}

    @staticmethod
    def _get_metrics(queue_collection: dict):
        """
        Obtains queue and instrument level fault information

        Results are keyed by queue/instrument and contain data of the form:
          * Problem: # of problem level faults observed
          * Fatal: # of fatal level faults observed
          * NR: #, difference between total and (problem + fatal)
          * total: # of faults observed
          * steps: # of steps encountered
          * wells: # of "wells" encountered (wells uses the wellplate with the most wells per queue
          not the sum of all wellplates' wells in a queue)

        :param queue_collection: The collection of queues being processed
        :return: results by queue, results by instrument, a copy of the fault records {queue_name: list}
        """
        q_results = {}
        fault_record = {}
        i_results = {}
        default_entry = {mcs.V_PROBLEM: 0,
                         mcs.V_FATAL: 0,
                         'NR': 0,
                         'total': 0,
                         'steps': 0,
                         'wells': {}}

        for q_name in queue_collection.keys():
            active_queue = queue_collection[q_name]

            well_accumulator = [0, ]
            for _, plate in active_queue[dbi.Q_CONTAINERS].items():
                contents: dict = plate.get(dbi.DBG_CONTENTS, {})
                if (not contents) or isinstance(contents, str):
                    continue
                well_accumulator.append(len(contents.keys()))
            n_wells = max(well_accumulator)

            fault_record[q_name] = active_queue.get(dbi.Q_RECORD, list())
            n_fatal_faults = sum([1 for f in fault_record[q_name] if f.get('level', None) == mcs.V_FATAL])
            n_problem_faults = sum([1 for f in fault_record[q_name] if f.get('level', None) == mcs.V_PROBLEM])
            n_busy_faults = sum([1 for f in fault_record[q_name] if f.get('level', None) == mcs.V_BUSY])
            n_total_faults = len(fault_record[q_name])

            q_results[q_name] = {
                'steps': len(active_queue[dbi.Q_OPERATIONS_LIST].keys()),
                'wells': n_wells,
                mcs.V_PROBLEM: n_problem_faults,
                mcs.V_FATAL: n_fatal_faults,
                'NR': n_total_faults - (n_problem_faults + n_fatal_faults + n_busy_faults),
                'total': n_total_faults,
            }

            for _, step in active_queue.get(dbi.Q_OPERATIONS_LIST, {}).items():
                agent = step.get(dbi.QOP_AGENT, 'None')
                i_results.setdefault(agent, deepcopy(default_entry))
                i_results[agent]['steps'] += 1
                i_results[agent]['wells'].setdefault(q_name, 0)
                try:
                    container_nickname = step[dbi.QOP_CONTAINER]
                    wpc = len(active_queue[dbi.Q_CONTAINERS][container_nickname][dbi.DBG_CONTENTS].keys())
                except (KeyError, TypeError, AttributeError):
                    pass
                else:
                    i_results[agent]['wells'][q_name] = max(wpc, i_results[agent]['wells'][q_name])
                finally:
                    wpc = 0  # noqa (it's in a for-loop)
            for fault in fault_record[q_name]:
                location = fault.get('location', 'None')
                i_results.setdefault(location, deepcopy(default_entry))
                level = fault.get('level', 'NR')
                i_results[location][level] += 1
                i_results[location]['total'] += 1
                i_results[location]['wells'].setdefault(q_name, 0)
        for _, v in i_results.items():
            v['wells'] = sum([wpc for _, wpc in v['wells'].items()])
        return q_results, i_results, fault_record

    def get_metrics(self, calculate_utilization=False):
        """
        Compiles all metrics into a dictionary

        :param calculate_utilization: A flag for if utilization metrics should be calculated (expensive)
        :return: (the metrics document: dict, a copy of the fault record: list)
        """
        q_results, i_results, fault_record = self._get_metrics(self.queue_collection)

        q_summary = [v.get(dbi.Q_STATUS, "NR") for _, v in self.queue_collection.items()]
        metric_document = {
            'Queue': {
                k: {
                    mcs.V_FATAL: v[mcs.V_FATAL], mcs.V_PROBLEM: v[mcs.V_PROBLEM], 'NR': v['NR'], 'total': v['total'],
                    'steps': v['steps'], 'wells': v['wells']
                } for k, v in q_results.items()
            },
            'Instrument': {
                k: {
                    mcs.V_FATAL: v[mcs.V_FATAL], mcs.V_PROBLEM: v[mcs.V_PROBLEM], 'NR': v['NR'], 'total': v['total'],
                    'steps': v['steps'], 'wells': v['wells']
                } for k, v in i_results.items()
            },
            'Summary': {
                'n_q_with_fatal': len([1 for _, v in q_results.items() if v[mcs.V_FATAL]]),
                'n_q_with_problem': len([1 for _, v in q_results.items() if v[mcs.V_PROBLEM]]),
                'n_q_with_fault': len([1 for _, v in q_results.items() if v['total']]),
                'n_steps_total': sum([v['steps'] for _, v in q_results.items()]),
                'n_wells_total': sum([v['wells'] for _, v in q_results.items()]),
                'n_complete_q': q_summary.count(dbi.DBQS_DONE),
                'n_idle_q': q_summary.count(dbi.DBQS_IDLE),
                'n_in-progress_q': q_summary.count(dbi.DBQS_WORK),
                'n_error_q': q_summary.count(dbi.DBQS_FAIL),
                'n_not-reported_q': q_summary.count('NR'),
                'n_fatal_faults': sum([1
                                       for _, fault_list in fault_record.items()
                                       for fault in fault_list
                                       if fault.get('level', None) == mcs.V_FATAL]),
                'n_problem_faults': sum([1
                                         for _, fault_list in fault_record.items()
                                         for fault in fault_list
                                         if fault.get('level', None) == mcs.V_PROBLEM]),
                'n_total_faults': sum([len(fault_list) for _, fault_list in fault_record.items()])
            }
        }

        if calculate_utilization:
            optimal_time, max_reasonable_time, actual_time, rest_time, *_ = self.get_utilization()
            metric_document['Utilization'] = {
                'optimal_time': optimal_time,
                'max_reasonable_time': max_reasonable_time,
                'actual_time': actual_time,
                'rest_time': rest_time
            }

        return metric_document, fault_record

    def get_report(self, print_record='i', calculate_utilization=False, push=False):
        """
        Converts the metrics document into a printable report

        :param print_record: 'i' print an abbreviated fault record by instrument, 'q' print an abbreviated fault record
        by queue, else do not print fault record
        :param calculate_utilization: A flag for if utilization metrics should be calculated (expensive)
        :param push: A flag to update the metric log file
        :return: A "pretty" string reporting a summary of the results
        """
        metric_document, fault_record = self.get_metrics()
        if push:
            self.update_log(calculate_utilization, override=metric_document)

        section_break = f"+{'-' * (WIDTH - 1)}+"
        report = [f"#{'=' * (WIDTH - 1)}#", ]

        def add_to_report(title, item=None, _frmt=pformat):
            report.append(f"+ {title} {'-' * max(0, WIDTH - 3 - len(title))}+")
            if item:
                report.append(_frmt(item, width=WIDTH))
            report.append(section_break)

        def safe_divide(numerator, denominator):
            try:
                return round(1 - (numerator / denominator), 3)
            except ZeroDivisionError:
                return f"ZDE({numerator})"

        def fault_breakdown(v):
            return f"(fatal: {v[mcs.V_FATAL]:2d}, problem: {v[mcs.V_PROBLEM]:2d}, " \
                   f"nr: {v['NR']:2d}; total: {v['total']:2d})"

        def basis_breakdown(v):
            return f"(steps: {v['steps']:3d}, wells {v['wells']:4d}) " \
                   f"[{safe_divide(v[mcs.V_FATAL], v['steps'])}, {safe_divide(v['total'], v['steps'])}]"

        add_to_report(f"Performance Report ({datetime.now().strftime(mcs.TIME_FORMAT)})")

        if print_record == "i":
            _fault_record = dict()
            for q_name, fault_list in fault_record.items():
                for fault_record in fault_list:
                    loc = fault_record.get('location', 'None')
                    _fault_record.setdefault(loc, dict())
                    _fault_record[loc].setdefault('faults', list())
                    abbreviated_entry = {'data': fault_record.get('data', None),
                                         'level': fault_record.get('level', None),
                                         'queue': q_name,
                                         'step': fault_record.get('step', None),
                                         'timestamp': fault_record.get('timestamp', None)}
                    _fault_record[loc]['faults'].append(abbreviated_entry)
            add_to_report("Fault Record", _fault_record)
        elif print_record == 'q':
            _fault_record = dict()
            for _, fault_list in fault_record.items():
                for fault in fault_list:
                    fault.pop('file_locator', None)
            add_to_report("Fault Record", fault_record)

        add_to_report("Queues Analyzed", self._queue_header())

        add_to_report("Faults by Queue",
                      {k: f"{pad(k, metric_document['Queue'])}{fault_breakdown(v)} {basis_breakdown(v)}"
                       for k, v in metric_document['Queue'].items()
                       })
        add_to_report("Faults by Instrument",
                      {k: f"{pad(k, metric_document['Instrument'])}{fault_breakdown(v)} {basis_breakdown(v)}"
                       for k, v in metric_document['Instrument'].items()
                       })

        q_summary_1 = f"complete: {metric_document['Summary']['n_complete_q']}, " \
                      f"idle: {metric_document['Summary']['n_idle_q']}, " \
                      f"in-progress: {metric_document['Summary']['n_in-progress_q']}, " \
                      f"error: {metric_document['Summary']['n_error_q']}, " \
                      f"not-reported: {metric_document['Summary']['n_not-reported_q']}"
        q_summary_2 = f"fatal: {metric_document['Summary']['n_q_with_fatal']}, " \
                      f"problem: {metric_document['Summary']['n_q_with_problem']}; " \
                      f"either: {metric_document['Summary']['n_q_with_fault']}"
        q_summary_3 = f"fatal: {metric_document['Summary']['n_fatal_faults']}, " \
                      f"problem: {metric_document['Summary']['n_problem_faults']}; " \
                      f"total: {metric_document['Summary']['n_total_faults']}"
        q_summary_4 = safe_divide(metric_document['Summary']['n_fatal_faults'],
                                  metric_document['Summary']['n_steps_total'])
        q_summary_5 = safe_divide(metric_document['Summary']['n_total_faults'],
                                  metric_document['Summary']['n_steps_total'])

        if calculate_utilization:
            util_report = "\n| Utilization: " + self.get_utilization_report(
                source=metric_document.get('Utilization', None),
                _frmt=lambda x: pformat(x, width=WIDTH)
            )
        else:
            util_report = ""

        add_to_report("Summary",
                      f"| # Queues: {len(self.queue_collection)} ({q_summary_1})\n"
                      f"| # Queues with faults: ({q_summary_2})\n"
                      f"| # steps: {metric_document['Summary']['n_steps_total']}\n"
                      f"| # wells: {metric_document['Summary']['n_wells_total']}\n"
                      f"| # faults: ({q_summary_3})"
                      f" [{q_summary_4}, {q_summary_5}]"
                      f"{util_report}",
                      _frmt=lambda x, **__: x
                      )
        report = "\n".join(report)
        cleaners = [("\n '", "\n| "),  # Makes LHS into "| "
                    ("\n{'", "\n| "),  # Catches exception to finish LHS
                    ("}\n", "\n"),     # Removes terminal "}" to match "{" removal
                    ("'", ""),         # ' is not needed
                    ('"', ""),         # " is not needed
                    (",\n", "\n")]    # Terminal , not needed
        for c in cleaners:
            report = report.replace(*c)
        return report

    # Running # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    def run(self, interval_hour=8, terminator: int = None):
        """ Periodically calls the MetricMonitor to calculate and report reliability and availability

        :param interval_hour: How many hours to wait between updates
        :param terminator: None - run infinitum, positive int - runs that many times, negative int - run until signaled
        :return: None
        """
        _iter = 0
        if terminator <= 0:
            watcher = Watcher(self.name)
        else:
            watcher = None

        while True:
            self.update_log(self.calc_util)

            if terminator is None:
                time.sleep(interval_hour * 60 * 60)
            elif terminator > 0:
                _iter += 1
                if _iter > terminator:
                    break
                time.sleep(interval_hour * 60 * 60)
            elif watcher:
                watcher.wait_for_change()
            else:
                break

    def update_log(self, calculate_utilization, override=None):
        """
        Pushes a summary of the performance metrics to a log file

        :param calculate_utilization: A flag to indicate if utilization should be calcualted
        :param override: If not provided, self.get_metrics is called, otherwise override treated as the metrics document
        :return:
        """
        if override:
            metric_document = override
        else:
            metric_document, *_ = self.get_metrics(calculate_utilization)
        report = {
            'Queues': [k for k in metric_document["Queue"]],
            'Summary': {VERB[k]: v for k, v in metric_document["Summary"].items()}
        }
        report['Summary'].update({VERB[k]: v for k, v in metric_document.get("Utilization", {}).items()})

        metric_log.info("MCN Performance Metrics:\n" + pformat(report))


class Watcher:
    """
    Used to monitor the database for queue creation or moves to historical
    """
    def __init__(self, name='MC'):
        """

        :param name: For database credentials
        """
        self.name = name
        self.previous_q_col = None

    def detect_change_in_queues(self):
        """
        Inspects the current queue collection and compares to the previous queue collection
        :return: True - a change was detected, False - no change or database connection failed
        """
        try:
            doc, _, _ = dbi.query_collection(self.name, dbi.COL_QUEUE)
        except cexc.DatabaseRequestError:
            doc = None

        if doc is None:
            return False

        new_q_col = sorted([k for k in doc])
        change_detected = new_q_col != self.previous_q_col
        self.previous_q_col = new_q_col
        return change_detected

    def wait_for_change(self):
        """
        Checks every 60 seconds for a change in the queue collection
        :return:
        """
        while True:
            if self.detect_change_in_queues():
                return True
            time.sleep(60)


if __name__ == '__main__':
    sel_queues = "sulfur"

    mm = MetricMonitor(name="MC", trimming=2)
    mm.set_queue_to_database(include_active=True)
    mm.down_select_queues_by_campaign(sel_queues)
    mm.ignore_user_errors()
    print(mm.get_report(print_record='n', calculate_utilization=True, push=False))

    """
    '23':
        agent: MC.__
        completed: 'yes'
        container: ''
        details: {}
        end_time: 10/05/2022 17:21:42
        operation: soft_wait
        start_time: 10/04/2022 23:56:44
        time_est: 3600
    '24':
        agent: Lh
        operation: detect_liquid_level
        container: fraction_plate
        time_est: 75
        start_time: 10/05/2022 17:22:02
        end_time:   10/05/2022 17:22:35  # Only took 33 seconds...
        details:
            detection_well: H12
            volume_history:
                - - 10/04/2022 10:37:03
                  - 604.0
                - - 10/04/2022 16:05:25
                  - 597.0
                - - 10/04/2022 17:49:50
                  - 581.0
                - - 10/04/2022 18:50:52
                  - 567.0
                - - 10/04/2022 19:51:58
                  - 554.0
                - - 10/04/2022 20:53:03
                  - 537.0
                - - 10/04/2022 21:54:12
                  - 527.0
                - - 10/04/2022 22:55:22
                  - 511.0
                - - 10/04/2022 23:56:23
                  - 504.0
                - - 10/05/2022 14:29:50
                  - 348.0
                - - 10/05/2022 17:22:32  # Timestamps align with end_time
                  - -1.0
        completed: 'yes'
    """
