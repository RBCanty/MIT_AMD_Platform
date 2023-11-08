# Spark_Controller.py
# Created 12 Feb2022
# Creator: Ben Canty
# Python Version 3.7.2
# Needs to be Running As Admin in order to access the Spark's API

# Some of the code for controlling the plate reader (Spark) is more complicated than necessary: We had attempted to
#   allow the SP manage multiple jobs at once (e.g. let the Spark take a measurement for Job A while a plate for Job B
#   is baking in the photo-reactor if time would allow).

# import glob
# project_dirs = glob.glob(r'C:\Users\admin\PycharmProjects\AMD_Control_Platform\venv\*\\')
# [sys.path.append(project_dir) for project_dir in project_dirs if project_dir not in sys.path]


import os
import shutil
import time
import traceback
from collections.abc import Iterable
from copy import copy
from datetime import datetime
from pprint import pformat
from threading import Lock, Thread

import custom_exceptions as cexc
import database_interface as dbi
import mcn_status as mcs
import project_processing as pp
import spectrometer_project as sp
from constants import WELL_IDS
from custom_classes import safe_open, Nop, ParcelEvent
from database_constants import DBG_CONTENTS
from functionals import realtime_wait, list_subtract
from mcn_logging_manager import system_log
from spark_API import SparkAutomation
from ui_quick_dialog import tk, QuickUI
from ui_spark_state import SparkStateUI

# Paths and Directories
DATABASE_PATH = r".\Shell_SP\Spark_API\References"
_SPARK_LOG_DIRECTORY = "./"


class SparkController:
    """ SparkController (class: Operational Manager)

    Handles control over the Tecan Spark and the raw data files it produces
    """

    def __init__(self, **kwargs):
        self._lock = Lock()
        self.event_logger = sp.EventLogger()
        self.spark = SparkAutomation(self.event_logger, **kwargs)
        self.jobs = dict()
        """
        A dictionary, keyed by job_id, for running jobs
        
        References: :class:`sp.SparkMethod`, :class:`sp.Appellomancer`, :class:`sp.SparkProjectFileHandler`
        
        Core Keys:
          * state:                  'active' or 'queued'
          * plate:                  wellplate collection document
          * method_obj:             SparkMethod
          * appello:                Appellomancer(initialized to queue_name and assay_type)
          * data_handler:           SparkProjectFileHandler(initialized from Appellomancer)
        
        Additional Keys:
          * excluded_wells:         List
          * error_codes:            List of strings of error codes
          * task_obj:               Spark API's result object
          * method_mngr:            Spark API's method object
          * task_thread:            Thread object (deleted by Janus)
          * tracker:                ParcelEvent object (only created if needed)
          * cancel:                 bool
        """

        self.kill_switch = False
        self.shutdown_success = False
        # self.janitor = Thread(target=self.janus)
        # self.janitor.start()

        self.data_tamer = pp.DataTamer(db_path=DATABASE_PATH)
        system_log.info("SparkControl Initialized")

    @property
    def active_job_ids(self):
        """
        :return: A list of all jobs whose state is 'active'
        """
        return [k for k, v in self.jobs.items() if v['state'] == 'active']

    @property
    def queued_job_ids(self):
        """
        :return: A list of all jobs whose state is 'queued'
        """
        return [k for k, v in self.jobs.items() if v['state'] == 'queued']

    ####################################################################################################################

    def janus(self):
        """
        Watches the Threads of active jobs, once complete it saves the job to a log file then deletes the thread and job

        :return: None
        """
        while True:
            if self.kill_switch:
                active_jobs = self.active_job_ids
                if active_jobs:
                    self.log_job(active_jobs)
                break
            else:
                time.sleep(5)
            hit_list = []
            for job_id in self.active_job_ids:
                job_obj = self.jobs[job_id]
                task_thread: Thread = job_obj.get('task_thread', None)
                if task_thread is None:
                    pass
                if not task_thread.is_alive():
                    del task_thread
                    job_obj.pop('task_thread', None)
                    job_obj.pop('tracker', None)
                    if self.log_job(job_id):
                        self.event_logger.excise_job(job_id)
                        hit_list.append(job_id)
                    else:
                        job_obj['task_thread'] = Thread()  # w/o target it is already dead, Thread.start() not needed
            for hit in hit_list:
                self.jobs.pop(hit, None)

    def log_job(self, job_id):
        """ Writes information pertaining to a given job_id to a log file.

        :param job_id: A Key for self.jobs
        :return: True - success, False - FileNotFoundError or PermissionError
        """
        try:
            spark_log_file = f"Spark_job_{job_id}.log"
            try:
                log_directory = self.jobs[job_id]['appello'].get_project_directory()
            except:  # noqa
                log_directory = _SPARK_LOG_DIRECTORY
            with safe_open(os.path.join(log_directory, spark_log_file), 'a+') as fh:
                if isinstance(job_id, str):
                    fh.write(self.print_job(job_id, True))
                elif isinstance(job_id, Iterable):
                    fh.write("Spark Control Terminated with active jobs:")
                    for _job_id in job_id:
                        fh.write(f"\t{_job_id}\n")
                fh.write("\n")
        except (FileNotFoundError, PermissionError):
            return False
        else:
            return True

    ####################################################################################################################

    def define_method(self, queue_name, plate_doc, job_id, **kwargs):
        """
        Sets the method details for a job (self.job[job_id])

        :param queue_name: Name of queue (used to create Appellomancer)
        :param plate_doc: Used to select the base well ids
        :param job_id: Key for 'self.job' property
        :param kwargs: If the lookup of the project specified by the Appellomancer fails, these kwargs will be used to
        specify the method details; additional 'excluded_wells' key can be used to further exclude wells.
        :return: None
        """
        assay = kwargs.get('assay', "assay-not-specified")
        name_wizard = sp.Appellomancer(queue_name, assay)
        data_handler = sp.SparkProjectFileHandler(name_wizard.get_path_to_project())
        extant_details = data_handler.get_spark_method_details()
        if not extant_details:
            method_object = sp.SparkMethod(**kwargs)
            data_handler.set_spark_method_details(details=method_object)
        else:
            extant_details.update(kwargs)
            method_object = sp.SparkMethod(**extant_details)

        sel_plate_cont = plate_doc.get(DBG_CONTENTS, "Empty")
        if isinstance(sel_plate_cont, str):
            sel_plate_cont = dict()
        selected_wells = list(sel_plate_cont.keys())
        empty_wells = list_subtract(WELL_IDS, selected_wells)
        excluded_wells = kwargs.get('excluded_wells', [])
        x_wells = empty_wells + list_subtract(excluded_wells, empty_wells)

        self.jobs[job_id] = {
            'state': 'queued',
            'plate': plate_doc,
            'method_obj': method_object,
            'appello': name_wizard,
            'data_handler': data_handler,
            'manifest': list(),
            'excluded_wells': x_wells
        }
        return

    ####################################################################################################################

    def print_job(self, job_id: str, tabbed=False):
        """ Used to create a printable for a job

        The printable consists of:
          * header - job_id
          * method_details - measurement mode
          * data_output - Location of the project file
          * error_details - self.jobs[job_id]['error_codes'] or 'No error codes recorded'
          * loop_details - Loop properties from the method object
          * spectral_details - Spectral properties from the method object
          * temperature_details - Temperature properties from the method object
          * manifest - Extract of events logged for the given job or "<No entries found in manifest>"

        Which are then joined by the delimiter.

        :param job_id: A key for self.jobs
        :param tabbed: True - use newline-tab otherwise use newline for the delimiter
        :return: A string representation of a job's log
        """
        job_manifest_entries = self.event_logger.get_job(job_id)

        if job_id not in self.jobs:
            return f"{job_id} - Not Found"

        if tabbed:
            delimiter = "\n\t"
        else:
            delimiter = "\n"

        header = f"{job_id}"
        method_details = "method mode: " + str(self.jobs[job_id].get('method_obj', Nop()).mode)
        name_wizard = self.jobs[job_id].get('appello', Nop())
        data_output = "Project location: " + str(name_wizard.get_project_directory())
        loop_details = f"loop_properties: {self.jobs[job_id]['method_obj'].get_loop_properties()}"
        spectral_details = f"spectral_properties: {self.jobs[job_id]['method_obj'].get_spectral_properties()}"
        temperature_details = f"temperature_properties: {self.jobs[job_id]['method_obj'].get_temperature_properties()}"
        error_details = str(self.jobs[job_id].get('error_codes', 'No error codes recorded'))

        manifest = "<No entries found in manifest>" \
            if not job_manifest_entries \
            else delimiter.join(f"{v[0]} - {v[1]}" for v in job_manifest_entries)

        printable = [header, method_details, data_output, error_details,
                     loop_details, spectral_details, temperature_details,
                     manifest]

        return delimiter.join(printable)

    @staticmethod
    def atropos(str_source: str, str_target: str, size: int, out_type=None):
        """
        Cuts a string of size (size) from a larger string, located by the prefix (str_target).
        The string is then split by " " and the first element is taken.
        This excision has any quotations and spaces removed.

        :param str_source: A string being searched
        :param str_target: The prefix for the section desired
        :param size: The expected size of the desired string
        :param out_type: converts the string into the desired type
        :return: The string prefixed by the target within the source
        """
        try:
            start_index = str_source.index(str_target)
        except ValueError:
            return None
        start_index += len(str_target)
        excision = str_source[start_index:]
        excision = excision[:size]
        excision, *_ = excision.split(" ")
        excision = excision.translate({ord(c): None for c in " '\""})
        if out_type:
            return out_type(excision)
        else:
            return excision

    ####################################################################################################################

    def cancel(self, job_id):
        if job_id not in self.active_job_ids:
            return False, f"{job_id} not in active_jobs"
        self.event_logger.scan_cancel(self, job_id)
        self.jobs[job_id]['cancel'] = True
        if 'method_mngr' not in self.jobs[job_id]:
            return False, f"method_mngr not in {job_id}"
        try:
            self.spark.cancel(self.jobs[job_id]['method_mngr'])
        except Exception as mmce:
            print(repr(mmce))
            traceback.print_tb(mmce.__traceback__)
            return None, "See exception printout"
        time.sleep(2)
        job_obj = self.jobs[job_id]
        task_thread: Thread = job_obj.get('task_thread', None)
        if isinstance(task_thread, Thread) and (not task_thread.is_alive()):
            del task_thread
            job_obj.pop('task_thread', None)
            job_obj.pop('tracker', None)
            if self.log_job(job_id):
                self.event_logger.excise_job(job_id)
                self.jobs.pop(job_id, None)
        return True, "Cancelled"

    def is_cancelled(self, job_id):
        return self.jobs[job_id].get('cancel', False)

    def has_active_jobs(self):
        # TODO: Make non-blocking... (not sure why it _is_ blocking though)
        return bool(self.active_job_ids)

    @staticmethod
    def preheat(target=20, window=0.5):
        method_object = sp.SparkMethod(temperature=target, tol=window, off_when_done=False, wait_for_t=False)
        method_object.generate_method_file(override='preheat')
        name_wizard = sp.Appellomancer('na', 'preheat')
        new_job = {
            'state': 'queued',
            'appello': name_wizard,
            'method_obj': method_object,
            'manifest': copy(method_object.history),
            }
        return new_job

    def execute_job(self, job_id):
        """
        Spawns a thread to run a job.

        :param job_id: (str)
        :return: True/False (success), Str (message)
        """

        # Put the plate on the tray, and take it into the Spark
        if not self.spark.state['encumbered']:
            system_log.warning("Plate may be missing from tray!")
        self.spark.move_plate_in()

        # Move job from queued to active
        self.jobs[job_id]['state'] = 'active'

        # Start the job
        self.jobs[job_id]['task_thread'] = Thread(target=self.job_manager, args=(job_id,), daemon=True)
        self.jobs[job_id]['task_thread'].start()

        return True, f"Job Started [{job_id}]"

    def job_manager(self, job_id):
        """
        Threaded function used to manage the scan or loop of scans.

        :param job_id: (str)
        :return: None
        """
        active_job = self.jobs[job_id]
        loop_properties = active_job['method_obj'].get_loop_properties()
        q_name = active_job['appello'].queue_name
        n_iterations = loop_properties.get('n_iterations', 1)
        interval = loop_properties.get('interval', 0)
        is_photodeg = loop_properties.get('is_photo', False)
        photo_temp = active_job['method_obj'].get_photo_temperature()
        offset = loop_properties.get('offset', 0)
        interval = interval if interval else 0
        is_good = False
        self.event_logger.job_start(self, job_id)
        _start_iter = 0 + offset
        _end_iter = n_iterations + offset - 1

        try:
            if is_photodeg:
                self.validate_lamp(job_id, q_name, True)
            if photo_temp > 0:
                # self.spark.photoreactor.valves.open_valve()  # # Valve not installed, the air is just always on # #
                self.validate_heater(job_id, q_name, photo_temp)
        except RuntimeError:
            system_log.exception("Error in photodeg pre-run validation")

        for iteration in range(_start_iter, _end_iter + 1):
            if self.is_cancelled(job_id):
                system_log.info(f"User cancelled job '{job_id}' before iteration #{iteration}")
                break
            try:
                self.run_subjob(job_id, iteration, transient=(iteration != _end_iter))
            except cexc.BreakOuterLoop:
                system_log.exception(f"Job '{job_id}' terminated")
                break
            if not (iteration == _end_iter):
                # prepare for next iteration
                if is_photodeg:
                    system_log.info(f"Loading Photoreactor")
                    self.spark.transfer_to_photoreactor(job_id)

                timer = datetime.now()
                try:
                    system_log.info(f"Waiting {interval} minutes for iteration {iteration + 1}")
                    while realtime_wait(timer, minutes=interval):
                        if self.is_cancelled(job_id):
                            raise cexc.BreakOuterLoop(f"User cancelled job {job_id}")
                        time.sleep(5)
                except cexc.BreakOuterLoop:
                    system_log.info(f"User cancelled job '{job_id}' after iteration #{iteration}")
                    break

                try:
                    is_photodeg and self.spark.receive_from_photoreactor(job_id)
                except RuntimeError as re:
                    self.event_logger.job_fail(self, job_id, notes=re)
                    raise re
            else:
                # On final operation, do not prepare for next iteration
                pass
        else:
            is_good = True

        if not is_good:
            system_log.warning(f"Execution of job {job_id} was not good")
            self.event_logger.job_fail(self, job_id)

        closing_exceptions = [None, None, None]
        try:
            if is_photodeg:
                system_log.info("Turning Lamp Off")
                self.validate_lamp(job_id, q_name, False)
        except Exception as lamp_e:
            closing_exceptions[0] = lamp_e
        try:
            system_log.info("Turning Heater Off")
            self.validate_heater(job_id, q_name, 5)
        except Exception as heater_e:
            closing_exceptions[1] = heater_e
        # # Valve not installed, may not be installed, the air is just always on # #
        # try:
        #     system_log.info("Closing Valve")
        #     if closing_exceptions[1] is None:
        #         self.spark.photoreactor.valves.close_valve()
        #     else:
        #         raise RuntimeError("Cannot safely turn off valve if heater may be on")
        # except Exception as valve_e:
        #     closing_exceptions[2] = valve_e

        if any(closing_exceptions):
            if closing_exceptions[0]:
                system_log.error("FAILED TO TURN OFF LAMP")
            if closing_exceptions[1]:
                system_log.error("FAILED TO TURN OFF HEATER")
            if closing_exceptions[2]:
                system_log.error("FAILED TO TURN OFF VALVE")
            self.event_logger.job_fail(self, job_id, notes=closing_exceptions)
        else:
            self.event_logger.job_end(self, job_id)

    def run_subjob(self, job_id, iteration, live_method=None, attempt=0, transient: bool = None):
        """
        For each iteration or retry, this call function handles the individual tasks in a job

        Call without options for the initial run, options are for retries

        :param job_id: (str)
        :param iteration: (int) Which iteration in the loop is this task
        :param live_method: The working, reduced method file for a retry
        :param attempt: An index for tracking retries
        :param transient: If this scan is terminal or transient
        :return:
        :raises cexc.BreakOuterLoop:
        """
        if self.is_cancelled(job_id):
            system_log.info(f"User cancelled job '{job_id}' during iteration #{iteration}")
            raise cexc.BreakOuterLoop(f"User canceled job: {job_id}")

        # If initial use the job-specified values, otherwise this is a retry and use the inherited versions
        if live_method is None:
            x_wells = self.jobs[job_id].get('excluded_wells', None)
            live_method = self.jobs[job_id]['method_obj'].generate_method_file(excluded_wells=x_wells,
                                                                               transient=transient)

        system_log.info(f"Running {job_id}:{iteration} ({attempt})")
        try:
            self.event_logger.scan_start(self, job_id)
            results = self.spark.run_plate(self.jobs[job_id], live_method)
        except RuntimeError:
            raise cexc.BreakOuterLoop("Spark not operable")
        else:
            completed_wells = self._process_data_file(self.jobs[job_id], iteration, attempt)

        system_log.info(pformat(results))
        self.event_logger.log_scan_complete(self, pformat(results).replace("\n", "\n" + "\t" * 7), job_id)
        if results['post'][1] == "Finish":
            self.event_logger.scan_end(self, job_id)
            system_log.info(f"Finished {job_id}:{iteration}")
            return
        elif results['post'][1] == "Fail":
            self.event_logger.scan_partial(self, job_id)
            if results['post'][2]:
                system_log.warning(f"Spark failed a scan while partially operational, aborting {job_id}")
                raise cexc.BreakOuterLoop("Spark unstable")
            system_log.info(f"Retrying subset of {job_id}:{iteration}")
            if attempt < 24:
                time.sleep(5)
                live_method = sp.SparkMethod.deselect_wells(live_method, completed_wells)
                self.run_subjob(job_id, iteration, live_method, attempt + 1)
            else:
                system_log.error(f"More than {attempt} attempts made on {job_id}:{iteration}, skipping rest of plate")
        elif results['post'][1] == 'Cancel':
            raise cexc.BreakOuterLoop(f"User canceled job: {job_id}")
        else:
            self.event_logger.scan_fail(self, job_id)
            system_log.warning(f"Unhandled state for {job_id}:{iteration}: '{results['post']}'")

    def validate_lamp(self, job_id, q_name, lamp_on):
        """ Ensures the Lamp is in the desired state via call to controller or, failing that, a call on the user

        :param job_id:
        :param q_name:
        :param lamp_on:
        :return: None
        :raises: RuntimeError (if user fails to recover), DatabaseRequestError (if fails to post PROBLEM to DB)
        """
        lamp_exception = self.spark.photoreactor.set_lamp(lamp_on)
        if lamp_exception:
            try:
                raise lamp_exception
            except IOError:
                word = "on"*lamp_on + "off"*(not lamp_on)
                prompt = f"Failed to validate that the lamp is {word} (asking user to recover: {job_id})"
                dbi.record_fault(mcs.Fault(location="Pr",
                                           level=mcs.V_PROBLEM,
                                           data="User had to manually intervene with Solar Generator",
                                           queue=q_name))
                system_log.exception(prompt)
                root = tk.Tk()
                popup = QuickUI(root, title="Lamp Selector", dialog=prompt, buttons={}, ret_if_ok="ok")
                ret_val = popup.run()
                if ret_val == "ok":
                    pass
                else:
                    if lamp_on:
                        raise RuntimeError("User unable to recover lamp")
                    else:
                        system_log.warning("Lamp may still be on!")

    def validate_heater(self, job_id, q_name, photo_temp):
        """ Ensures the heater is in the desired state via call to controller or, failing that, a call on the user

        :param job_id:
        :param q_name:
        :param photo_temp:
        :return:
        """
        # If heat is required, make sure gas is flowing
        if photo_temp > 10:
            if not self.spark.photoreactor.is_gas_flowing():
                raise RuntimeError(f"Gas is not flowing, cannot start heater ({q_name})")

        # Set the setpoint (must be a positive number if pre-job, and will be 5 post-job)
        self.spark.photoreactor.set_temperature(photo_temp)

        # If the temperature is being set, wait for the temperature to be reached
        if photo_temp > 10:
            system_log.info("Waiting for photoreactor to reach temperature")
            status_code, history = self.spark.photoreactor.wait_for_temperature()
            if status_code == "Error":
                raise RuntimeError(f"Heater failed to reach setpoint temperature ({q_name})")
            system_log.info("Photoreactor temperature reached")

        # For logging
        if photo_temp > 5:
            # The temperature is being set
            self.jobs[job_id]['tracker'] = ParcelEvent()
            self.spark.photoreactor.start_temperature_tracker(job_id, self.jobs[job_id]['tracker'])
        else:
            # There is either no temperature control or the temperature control is being turned off
            if 'tracker' in self.jobs[job_id]:
                tracker = self.jobs[job_id]['tracker']
                tracker.set()
                # ^ This will cause the thread started by start_temperature_tracker to terminate
                # With a logger, this will copy the temperature history to the log

    @staticmethod
    def _process_data_file(job, iteration: int, attempt: int):
        """ Moves experimental data into the project (and project folder), removing it from the the original location

        :param job:
        :param iteration:
        :param attempt:
        :return: Completed wells (list)
        """
        src_output_file = job['appello'].get_path_to_source_data(None)
        if src_output_file is None:
            return []
        dst_output_file = job['appello'].generate_raw_data_dst_name(iteration, attempt)
        completed_wells = job['data_handler'].add_data(src_output_file, iteration, job['plate'])
        # ^ that is a SparkProjectFileHandler.add_data() not a dri.add_data()
        dst_output_file = job['appello'].generate_safe_name(dst_output_file, '.xlsx')
        shutil.move(src_output_file, dst_output_file)
        old_path, _ = os.path.split(src_output_file)
        old_path, _ = os.path.split(old_path)
        old_path, _ = os.path.split(old_path)
        try:
            shutil.rmtree(old_path)
        except PermissionError:
            pass
        return completed_wells

    def open_ui(self, pipes=None):
        """ Opens the spark control UI """
        prompt = SparkStateUI(tk.Tk(), self.spark, outbox=pipes)
        self.spark.set_instrument_state(**prompt.run())

########################################################################################################################
########################################################################################################################


if __name__ == '__main__':
    my_spark = SparkAutomation(event_logger=sp.EventLogger())
    print("Hello World")
    time.sleep(5)
    my_spark.move_plate_out(side="left")
    time.sleep(5)
    my_spark.move_plate_in()
    time.sleep(5)
    my_spark.shutdown()
    time.sleep(5)
    print("shutdown_success=", my_spark.shutdown_success)
    print("Goodbye World")
