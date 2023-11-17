""" Defines functionality of plate reader
@author: Ben C
"""

import re
import time
from copy import deepcopy
from datetime import datetime

import aceso
import ui_quick_dialog
from custom_exceptions import DatabaseRequestError
import database_interface as dbi
import functionals
import mcn_status as mcs
import project_processing as pp
import spectrometer_project as sp
from AMD_Spark_Control_v3 import SparkController
from custom_classes import Narg
from database_constants import *
from mcn_logging_manager import system_log
from constants import MSG_Q_OUT


@functionals.standard_exception
def run_platereader(*_, **kwargs) -> mcs.RetObj:
    """ Runs a plate reader method

    :param kwargs: Standard kwargs from MCN control
    :return: a RetObj
    """
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")

    if not kwargs.get('manual_mode', False):
        q_name, step_num, _, _ = functionals.unpack_db_kwargs(kwargs, d_pid="", d_a="")
        try:
            assay, plate_doc, method_kwargs = _lookup_spark_data_from_database(q_name, step_num)
        except (DatabaseRequestError, ValueError) as caught_exception:
            err_msg = f"Failed to retrieve requisite data from Database for {q_name}:{step_num}"
            system_log.exception(err_msg)
            return mcs.RetObj.incomplete('Pr', mcs.V_FATAL, err_msg + " " + str(caught_exception.args))
    else:
        q_name = kwargs['manual_mode'].get('q_name',
                                           f'manual_mode_{datetime.now().strftime(sp.TIME_FORMAT).replace(":", "-")}')
        step_num = 0
        assay, plate_doc, method_kwargs = _lookup_spark_data_from_kwargs(kwargs['manual_mode'])
    if 'offset' in kwargs:
        method_kwargs.setdefault('offset', kwargs['offset'])

    # Grab Spark
    with child_com as cc:
        spark_controller: SparkController = cc['Pr']['SparkControl']

    # Define Method for Spark
    try:
        job_name = f"{q_name}_{assay}_{step_num}"
        spark_controller.define_method(q_name, plate_doc, job_name, **method_kwargs)
    except ValueError:
        system_log.exception("Spark Method Loading error")
        return mcs.RetObj.incomplete('Pr', mcs.V_FATAL, f"Formatting Error for Method")

    # Run the Spark
    success, messages = spark_controller.execute_job(job_name)

    if success:  # method is started successfully (but end result is tbd)
        last_entry = [("na", f"Warning: Found no logs for {job_name}", job_name, "run_platereader"), ]
        for _ in range(72):
            time.sleep(5)
            job_logs = spark_controller.event_logger.get_job(job_name)  # (timestamp, message, job_id, caller,)
            if "Job End" in [x[1] for x in job_logs]:
                break
            if "Job Fail" in [x[1] for x in job_logs]:
                break
            if job_logs:
                last_entry = job_logs[-1]
        else:
            ui_quick_dialog.QuickUI(
                ui_quick_dialog.tk.Tk(),
                title="Spark method timeout",
                dialog="Press OK when Spark is done",
                buttons={}
            ).run()

        job_obj = spark_controller.jobs[job_name]
        task_thread = job_obj.get('task_thread', None)
        if hasattr(task_thread, 'is_alive') and (not task_thread.is_alive()):
            del task_thread
            job_obj.pop('task_thread', None)
        job_obj.pop('tracker', None)
        if spark_controller.log_job(job_name):
            job_logs = spark_controller.event_logger.excise_job(job_name)
            if job_logs:
                last_entry = job_logs[-1]
            spark_controller.jobs.pop(job_name, None)

        if (not last_entry) or ('Fail' in last_entry[1]):
            return mcs.RetObj.incomplete("Pr",
                                         mcs.V_FATAL,
                                         f"Spark job {job_name} failed with last known response: '{last_entry[1]}' "
                                         f"(See logs for more details)")
        if "Warning" in last_entry[1]:
            return mcs.RetObj.complete(f"Pr job {job_name} complete with Warnings: '{last_entry[1]}'")
        else:
            return mcs.RetObj.complete(f"Pr job {job_name} complete")
    else:
        system_log.info(f"Spark returned from run with: {messages} (See Pr log for more details)")
        return mcs.RetObj.incomplete("Pr", mcs.V_FATAL, messages)


@functionals.standard_exception
def check_platereader(*_, **kwargs) -> mcs.RetObj:
    """ Has the plate reader verify a signal strength for a plate and recommend corrections

    :param kwargs: Standard kwargs from MCN control
    :return: a RetObj
    """
    q_name, step_num, _, _ = functionals.unpack_db_kwargs(kwargs, d_pid="", d_a="")
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")

    try:
        main_doc, _, _ = dbi.query_document('SP', COL_QUEUE, 'queue_name', q_name)
    except DatabaseRequestError:
        exc_msg = "Pr DB request failed - Fetching Queue document"
        system_log.exception(exc_msg)
        return mcs.RetObj.incomplete("Pr", mcs.V_PROBLEM, exc_msg)
    queue_details = main_doc[Q_OPERATIONS_LIST][step_num][QOP_DETAILS]
    container = main_doc[Q_OPERATIONS_LIST][step_num][QOP_CONTAINER]
    assay = queue_details.get('assay', 'check')
    name_wizard = sp.Appellomancer(q_name, assay)
    project = sp.SparkProjectFileHandler(name_wizard.get_path_to_project())
    use_iteration = project.get_next_free_iteration()
    kwargs['offset'] = use_iteration

    scan_ret_obj = run_platereader(**kwargs)
    if not scan_ret_obj.completion:
        return scan_ret_obj

    with child_com as cc:
        spark_controller: SparkController = cc['Pr']['SparkControl']

    check_manager = pp.SignalCheck(project, spark_controller.data_tamer, **queue_details)
    check_manager.tamer.subtract_wellplate_background(project, on_iterations=[use_iteration, ])
    check_manager.tamer.subtract_background_by_solvent_volume(project, floor=0, on_iterations=[use_iteration, ])
    bad_wells = check_manager.check_bounds(None, use_iteration)
    target = check_manager.target

    job_name = f"{q_name}_{assay}_{step_num}"
    iter_num = int(queue_details.get('iter', -1))

    if not bad_wells:
        return mcs.RetObj.complete(f"Pr job {job_name} complete: no bad wells")
    if iter_num == 0:
        return mcs.RetObj.complete(f"Pr job {job_name} complete: max retry reached, bad wells (if any) ignored")

    track_changes = 0
    edited_doc = deepcopy(main_doc)
    try:
        edited_container = edited_doc[Q_CONTAINERS][container]
        temp_container_contents = edited_container.get(DBG_CONTENTS, "Empty")
        if isinstance(temp_container_contents, str):
            temp_container_contents = dict()

        # Grab highest iteration number:
        new_iter = 0
        for k, v in bad_wells.items():
            additional_cleanups = temp_container_contents.get(k, dict()).get('additional_cleanup', dict())
            existing_iters = [(('0',) + re.search(r"(\d+)(?!.*\d)", _k).groups())[-1]
                              for _k in additional_cleanups.keys()
                              if "spark_cleanup" in _k]
            candidate = max([int(n) for n in existing_iters], default=0)
            new_iter = max(candidate, new_iter)
        cleanup_name = f"spark_cleanup_{new_iter + 1}"

        for k, v in bad_wells.items():
            additional_cleanups = temp_container_contents.get(k, dict()).get('additional_cleanup', dict())
            if 'final' in additional_cleanups.keys():
                continue

            ratio = target / v
            if ratio > 1:
                cleanup_arg = {'concentration_factor': min(ratio, 50)}
            else:
                cleanup_arg = {'dilution_factor': ratio}

            if (DBG_CONTENTS not in edited_container) or ("Empty" == edited_container.get(DBG_CONTENTS, "Empty")):
                edited_container[DBG_CONTENTS] = dict()
            edited_container[DBG_CONTENTS].setdefault(k, dict())
            edited_container[DBG_CONTENTS][k].setdefault('additional_cleanup', dict())
            edited_container[DBG_CONTENTS][k]['additional_cleanup'].setdefault(cleanup_name, cleanup_arg)
            track_changes += 1
        if track_changes:
            system_log.info(f"Identified {track_changes} of {len(bad_wells)} for rectification")
            dbi.update_document('SP', COL_QUEUE, main_doc, edited_doc)
        else:
            return mcs.RetObj.complete(f"Pr job {job_name} complete: all bad wells were final")
    except DatabaseRequestError:
        exc_msg = "Pr DB request failed - Plate Cleanup Update"
        system_log.exception(exc_msg)
        return mcs.RetObj.incomplete("Pr", mcs.V_PROBLEM, exc_msg)

    _success, _msg, _exc = aceso.request_sample_adjustment(q_name, step_num, container, cleanup_name, iter_num, track_changes)
    if _success:
        system_log.info(_msg)
        return mcs.RetObj.complete(f"Pr job {job_name} complete: bad wells submitted for cleanup")
    else:
        if _exc:
            try:
                raise _exc
            except type(_exc):
                system_log.exception(_msg)
        else:
            system_log.info(_msg)
        return mcs.RetObj.incomplete("Pr", mcs.V_PROBLEM, _msg)


@functionals.standard_exception
def analyze_platereader(*_, **kwargs) -> mcs.RetObj:
    """ Process the results of a previous run

    :param kwargs: Standard kwargs from MCN control
    :return: a RetObj
    """
    q_name, step_num, _, _ = functionals.unpack_db_kwargs(kwargs, d_pid="", d_a="")
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    campaign_name = functionals.queue_name_to_campaign_name(q_name)

    with child_com as cc:
        spark_controller: SparkController = cc['Pr']['SparkControl']

    try:
        main_doc, _, _ = dbi.query_document('SP', 'queue', 'queue_name', q_name)
    except DatabaseRequestError:
        exc_msg = "Pr DB request failed - Fetching Queue document"
        system_log.exception(exc_msg)
        return mcs.RetObj.incomplete("Pr", mcs.V_PROBLEM, exc_msg)
    queue_details = main_doc[Q_OPERATIONS_LIST][step_num][QOP_DETAILS]
    assay = queue_details['assay']
    _assay = queue_details.get('from_assay', assay)
    name_wizard = sp.Appellomancer(q_name, _assay)
    project = sp.SparkProjectFileHandler(name_wizard.get_path_to_project())

    # Load the right manager
    system_log.info("Creating project manager")
    if assay == "peaks":
        project_manager = pp.Peaks(project, spark_controller.data_tamer)
    elif assay in ["oxideg", "photodeg"]:
        project_manager = pp.Kinetics(project, spark_controller.data_tamer)
    elif assay == "logd":
        project_manager = pp.LogD(project, spark_controller.data_tamer)
    elif assay in ["hyperpol", "logpi", "tox"]:
        return mcs.RetObj.incomplete("Pr", mcs.V_PROBLEM, f"Spark data processing not ready for assay '{assay}'")
    else:
        return mcs.RetObj.incomplete("Pr", mcs.V_PROBLEM, f"Spark data processing not defined for assay '{assay}'")

    # Subtract signal elements (only if not already done so)
    if not project.history_already_contains("subtract_wellplate_background"):
        system_log.info("Removing wellplate background")
        project_manager.tamer.subtract_wellplate_background(project)
    if not project.history_already_contains("subtract_background_by_solvent_volume"):
        system_log.info("Removing solvent background")
        project_manager.tamer.subtract_background_by_solvent_volume(project, floor=0)

    # Run the proper analysis
    system_log.info("Analyzing project")
    if assay == "peaks":
        project_manager.process_peaks(None, iteration=0, sort_by='wavelength', prominence=0.1)
        commit_kwargs = {'name': "SP", 'iteration': 0, 'peaks': True, 'spectra': True}
    elif assay in ["oxideg", "photodeg"]:
        if assay == "photodeg":
            project_manager.use_exposure_times(q_name, assay)
        else:
            project_manager.use_interval_times()
        project_manager.advanced_regression(assay, None, options=queue_details.get('regr_opt', None))
        commit_kwargs = {'name': "SP", 'assay_name': assay}
    elif assay == "logd":
        project_manager.process_log_d_by_integral(None)
        commit_kwargs = {'name': "SP"}
    else:  # should be unreachable
        return mcs.RetObj.incomplete("Pr", mcs.V_PROBLEM, f"Assay analysis fell through'{assay}'")

    # Commit data to project database
    system_log.info("Committing analysis results")
    _, failed_wells = project_manager.commit_to_database(**commit_kwargs, campaign=campaign_name)
    if failed_wells:
        system_log.warning(f"Database commit failed for: {failed_wells}")
        return mcs.RetObj.incomplete("Pr", mcs.V_PROBLEM, f"Database commit failed for {len(failed_wells)}")

    return mcs.RetObj.complete(f"Spark assay '{assay}' complete")


@functionals.standard_exception
def prepare_to_receive_and_send(*_, **kwargs) -> mcs.RetObj:
    """ Prepares the plate reader for a wellplate transfer

    :param kwargs: Standard kwargs from MCN control
    :return: a RetObj
    """
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    q_name, step_num, _, _ = functionals.unpack_db_kwargs(kwargs, d_qn="", d_sn="", d_pid="", d_a="Pr")
    with child_com as cc:
        spark_controller: SparkController = cc['Pr']['SparkControl']
    mode = kwargs.get(QOP_DETAILS, {}).get("direction", None)

    try:
        spark_controller.spark.reserve_spark()
        if not mode:
            return mcs.RetObj.incomplete("_U",
                                         mcs.V_PROBLEM,
                                         f"User did not specify 'direction' detail ({q_name}, {step_num})")
        elif mode.lower() == 'send':
            spark_controller.spark.prepare_send()
        elif mode.lower() == 'receive':
            spark_controller.spark.prepare_receive()
        else:
            return mcs.RetObj.incomplete("_U",
                                         mcs.V_PROBLEM,
                                         f"'{mode}' not recognized for 'direction' (expected 'send' or 'receive')\n"
                                         f"({q_name}, {step_num})")
    except (RuntimeError, RuntimeWarning):
        system_log.exception("")
        return mcs.RetObj.incomplete("Pr", mcs.V_FATAL, "Wellplate is mismanaged")
    finally:
        spark_controller.spark.release_spark()

    return mcs.RetObj.complete(f"Pr ready to {mode}")


@functionals.standard_exception
def return_to_initial_state(*_, **kwargs) -> mcs.RetObj:
    """ Returns the plate reader and its ancillary equipment to a default state

    :param kwargs: Standard kwargs from MCN control
    :return: a RetObj
    """
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    q_name, step_num, _, _ = functionals.unpack_db_kwargs(kwargs, d_qn="", d_sn="", d_pid="", d_a="Pr")
    with child_com as cc:
        spark_controller: SparkController = cc['Pr']['SparkControl']

    try:
        spark_controller.spark.reset_instrument()
    except (RuntimeError, RuntimeWarning):
        return mcs.RetObj.incomplete("Pr", mcs.V_FATAL, "Spark failed to return to initial state")

    return mcs.RetObj.complete("Spark returned to initial state")


@functionals.standard_exception
def preheat_platereader(*_, **kwargs) -> mcs.RetObj:
    """ Starts the platereader heating cycle (this can take a long time, so useful to call a few steps ahead)

    :param kwargs: Standard kwargs from MCN control
    :return: a RetObj
    """
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    q_name, step_num, _, _ = functionals.unpack_db_kwargs(kwargs, d_qn="", d_sn="", d_pid="", d_a="Pr")

    try:
        doc, _, _ = dbi.query_document("SP", COL_QUEUE, Q_NAME, q_name)
    except DatabaseRequestError:
        return mcs.RetObj.incomplete("Pr", mcs.V_PROBLEM, "Failed to access DB")
    temperature = doc[Q_OPERATIONS_LIST][step_num][QOP_DETAILS].get('temperature', 27)
    window = doc[Q_OPERATIONS_LIST][step_num][QOP_DETAILS].get('temperature', 7)

    with child_com as cc:
        spark_controller: SparkController = cc['Pr']['SparkControl']
        job_name = f"{q_name}_{step_num}_preheat_{datetime.now().strftime(DBQ_DATE_FORMAT).replace('/', '')}"
        spark_controller.jobs.update({job_name: spark_controller.preheat(temperature, window)})
        success, messages = spark_controller.execute_job(job_name)

    if success:
        return mcs.RetObj.complete(messages)
    else:
        return mcs.RetObj.incomplete('Pr', mcs.V_PROBLEM, messages)


@functionals.standard_exception
def is_spark_busy(*_, **kwargs) -> mcs.RetObj:
    """ Checks if the plate reader has active jobs

    :param kwargs: Standard kwargs from MCN control
    :return: a RetObj
    """
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    with child_com as cc:
        spark_controller: SparkController = cc['Pr']['SparkControl']
        return mcs.RetObj.complete(f"{spark_controller.has_active_jobs()}")


@functionals.standard_exception
def rebuild_plate_reader(*_, **kwargs) -> mcs.RetObj:
    """ Attempts to reconnect to plate reader

    :param kwargs: Standard kwargs from MCN control
    :return: a RetObj
    """
    _, internals, child_com = functionals.unpack_kwargs(kwargs, d_msg="")
    with child_com as cc:
        # Close/Clear old
        try:
            cc['Pr']['SparkControl'].spark.shutdown()
        except KeyError:
            pass
        else:
            cc['Pr'].pop('SparkControl', None)
        # Create new
        try:
            cc['Pr']['SparkControl'] = SparkController(
                **internals.get([mcs.S_AUXILIS, "SP"], {}))
        except:  # noqa
            system_log.exception("Failed to rebuild plate reader control")
            return mcs.RetObj.incomplete("Pr", mcs.V_FATAL, f"Failed to rebuild SparkController")
    return mcs.RetObj.complete("SparkController rebuilt")


@functionals.standard_exception
def cancel_job(*_, **kwargs) -> mcs.RetObj:
    """ Attempts to cancel a job

    Utility not checked as asynchronous operation never verified

    :param kwargs: Standard kwargs from MCN control
    :return: a RetObj
    """
    _, internals, child_com = functionals.unpack_kwargs(kwargs, d_msg="")
    job_name = kwargs.get("job_id", None)
    if job_name is None:
        q_name, step_num, _, _ = functionals.unpack_db_kwargs(kwargs, d_pid="", d_a="")
        try:
            doc, _, _ = dbi.query_document("SP", COL_QUEUE, Q_NAME, q_name)
        except DatabaseRequestError:
            return mcs.RetObj.incomplete("Pr", mcs.V_PROBLEM, "Failed to access DB")
        job_name = doc[Q_OPERATIONS_LIST][step_num][QOP_CONTAINER]
    with child_com as cc:
        success, messages = cc['Pr']['SparkControl'].cancel(job_name)
    if success:
        return mcs.RetObj.complete(f"{job_name} cancelled")
    elif success is None:
        return mcs.RetObj.incomplete("Pr", mcs.V_FATAL, messages)
    else:
        return mcs.RetObj.incomplete("Pr", mcs.V_PROBLEM, messages)


@functionals.standard_exception
def modify_spark_state(*_, **kwargs):
    """ Overrides internal state of the plate reader (and ancillary instruments)

    :param kwargs: Standard kwargs from MCN control
    :return: a RetObj
    """
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    with child_com as cc:
        spark_controller: SparkController = cc['Pr']['SparkControl']
    try:
        spark_controller.spark.set_instrument_state(
            tray_position=kwargs.get('tray_position', Narg()),  # Location of the tray: "In", "Right", or "Left" (str)
            thermostat=kwargs.get('thermostat', Narg()),  # If the thermostat is in use (number)
            encumbered=kwargs.get('encumbered', Narg()),  # If the plate is on the tray (bool)
            expose=kwargs.get('expose', Narg()),  # If the plate is on the forklift/in the photoreactor (bool)
        )
    except (ValueError, TypeError):
        system_log.exception("Bad arguments given to SparkAutomation.set_instrument_state()")
        return mcs.RetObj.incomplete("_U", mcs.V_PROBLEM, "Improperly defined SparkAutomation state variables")
    return mcs.RetObj.complete("Spark internal state updated")


@functionals.standard_exception
def spark_ui_spawn(*_, **kwargs):
    """ Spawns an instance of the GUI

    Change of state may still be buggy

    :param kwargs: Standard kwargs from MCN control
    :return: a RetObj
    """
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    with child_com as cc:
        spark_controller: SparkController = cc['Pr']['SparkControl']
    try:
        spark_controller.open_ui(kwargs.get(MSG_Q_OUT, None))
    except (TypeError, ValueError) as ce:
        return mcs.RetObj.complete(f"Spark state not set, bad value given: {repr(ce)}")
    return mcs.RetObj.complete("Spark GUI terminated")


@functionals.standard_exception
def lamp_on(*_, **kwargs):
    """ Turns solar generator on

    :param kwargs: Standard kwargs from MCN control
    :return: a RetObj
    """
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    with child_com as cc:
        spark_controller: SparkController = cc['Pr']['SparkControl']
    try:
        spark_controller.spark.photoreactor.lamp.on()
    except Exception as e:
        return mcs.RetObj.incomplete("Pr", mcs.V_FATAL, repr(e))
    return mcs.RetObj.complete("Lamp on")


@functionals.standard_exception
def lamp_off(*_, **kwargs):
    """ Turns solar generator off

    :param kwargs: Standard kwargs from MCN control
    :return: a RetObj
    """
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    with child_com as cc:
        spark_controller: SparkController = cc['Pr']['SparkControl']
    try:
        spark_controller.spark.photoreactor.lamp.off()
    except Exception as e:
        return mcs.RetObj.incomplete("Pr", mcs.V_FATAL, repr(e))
    return mcs.RetObj.complete("Lamp off")


@functionals.standard_exception
def raise_forklift(*_, **kwargs):
    """ Raises forklift in the photoreactor

    :param kwargs: Standard kwargs from MCN control
    :return: a RetObj
    """
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    with child_com as cc:
        spark_controller: SparkController = cc['Pr']['SparkControl']
    try:
        spark_controller.spark.photoreactor.lift_raise()
    except Exception as e:
        return mcs.RetObj.incomplete("Pr", mcs.V_FATAL, repr(e))
    return mcs.RetObj.complete("Forklift raised")


@functionals.standard_exception
def lower_forklift(*_, **kwargs):
    """ Lowers forklift in the photoreactor

    :param kwargs: Standard kwargs from MCN control
    :return: a RetObj
    """
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    with child_com as cc:
        spark_controller: SparkController = cc['Pr']['SparkControl']
    try:
        spark_controller.spark.photoreactor.lift_lower()
    except Exception as e:
        return mcs.RetObj.incomplete("Pr", mcs.V_FATAL, repr(e))
    return mcs.RetObj.complete("Forklift lowered")


@functionals.standard_exception
def set_photo_heater(*_, **kwargs):
    """ Sets heater setpoint for photoreactor

    :param kwargs: Standard kwargs from MCN control
    :keyword set: The setpoint (default: 5 <int>)
    :return: a RetObj
    """
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    setpoint = kwargs.get('set', 5)
    with child_com as cc:
        spark_controller: SparkController = cc['Pr']['SparkControl']
    try:
        spark_controller.spark.photoreactor.set_temperature(setpoint)
    except Exception as e:
        return mcs.RetObj.incomplete("Pr", mcs.V_FATAL, repr(e))
    return mcs.RetObj.complete(f"Heater set to '{setpoint}'")


@functionals.standard_exception
def move_tray(*_, **kwargs) -> mcs.RetObj:
    """ Moves the plate reader tray

    Can use kwarg "mode" or a detail "mode" for use with and without DB control
      - in  -> moves tray in
      - out -> moves tray out to default side (onlooker right, spark's left)
      - left -> moves tray out to onlooker's left side (spark's right)
      - right -> moves tray out to onlooker's right side (spark's left)

    :param kwargs: Standard kwargs from MCN control
    :return: a RetObj
    """
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    with child_com as cc:
        spark_controller: SparkController = cc['Pr']['SparkControl']
    mode = kwargs.get("mode", "").lower()
    if not mode:
        mode = kwargs.get(QOP_DETAILS, {}).get("mode", "in").lower()

    try:
        spark_controller.spark.reserve_spark()
        if mode == "in":
            spark_controller.spark.move_plate_in()
        elif mode == 'out':
            spark_controller.spark.move_plate_out()
        elif mode == 'right':
            spark_controller.spark.move_plate_out(side="right")
        elif mode == 'left':
            spark_controller.spark.move_plate_out(side="left")
        else:
            return mcs.RetObj.incomplete("_U", mcs.V_PROBLEM, f"Detail incorrectly specified (mode={mode})")
    except (RuntimeError, RuntimeWarning):
        system_log.exception("Error moving plate tray")
        return mcs.RetObj.incomplete("Pr", mcs.V_FATAL, "Error moving plate tray")
    finally:
        spark_controller.spark.release_spark()

    return mcs.RetObj.complete(f"Pr ready to {mode}")


@functionals.standard_exception
def read_anemometer(*_, **kwargs) -> mcs.RetObj:
    """ Checks value of gas flow

    :param kwargs: Standard kwargs from MCN control
    :return: a RetObj
    """
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    with child_com as cc:
        spark_controller: SparkController = cc['Pr']['SparkControl']

    try:
        spark_controller.spark.reserve_spark()
        anemometer_reading = spark_controller.spark.photoreactor.anemometer.read_wind_speed()
        anemometer_msg = f"Anemometer reads: {anemometer_reading}"
        system_log.info(anemometer_msg)
        return mcs.RetObj.complete(anemometer_msg)
    except (RuntimeError, RuntimeWarning):
        err_msg = "Error reading anemometer"
        system_log.exception(err_msg)
        return mcs.RetObj.incomplete("Pr", mcs.V_FATAL, err_msg)
    finally:
        spark_controller.spark.release_spark()


@functionals.standard_exception
def calibrate_anemometer(*_, **kwargs) -> mcs.RetObj:
    """ Establishes calibration for gas flow

    :param kwargs: Standard kwargs from MCN control
    :keyword n: The number of points to use in calibration
    :keyword override: The calibrated zero-point
    :return: a RetObj
    """
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    with child_com as cc:
        spark_controller: SparkController = cc['Pr']['SparkControl']
    n = kwargs.get('n', None)
    override = kwargs.get('override', None)
    czp_kwargs = dict()
    if n is not None:
        czp_kwargs['n'] = int(n)
    if override is not None:
        czp_kwargs['override'] = float(override)

    try:
        spark_controller.spark.reserve_spark()
        anemometer_reading = spark_controller.spark.photoreactor.anemometer.calibrate_zero_point(**czp_kwargs)
        anemometer_msg = f"Anemometer calibrated to: {anemometer_reading}"
        system_log.info(anemometer_msg)
        return mcs.RetObj.complete(anemometer_msg)
    except (RuntimeError, RuntimeWarning):
        err_msg = "Error calibrating anemometer"
        system_log.exception(err_msg)
        return mcs.RetObj.incomplete("Pr", mcs.V_FATAL, err_msg)
    finally:
        spark_controller.spark.release_spark()


def _lookup_spark_data_from_database(q_name, step_num):
    # Grab Queue document for details
    try:
        main_doc, _, _ = dbi.query_document('SP', 'queue', 'queue_name', q_name)
    except DatabaseRequestError:
        raise DatabaseRequestError("Pr DB request failed - Method type")

    method_kwargs = main_doc[Q_OPERATIONS_LIST][step_num][QOP_DETAILS]
    assay = method_kwargs.get('assay', 'assay-not-specified')

    # What is the plate we are scanning called?
    try:
        container_nickname = main_doc[Q_OPERATIONS_LIST][step_num][QOP_CONTAINER]
        plate_container_name = main_doc[Q_CONTAINERS][container_nickname][DBG_CONTAINER_NAME]
    except KeyError as ke:
        raise KeyError(f"Pr missing plate information: "
                       f"Queue '{q_name}', step #{step_num} "
                       f"[KeyError on '{ke.args[0]}']")

    # Grab Wellplate document for contexts
    try:
        plate_doc, _, _ = dbi.query_document('SP', COL_WELLPLATES, DBG_CONTAINER_NAME, plate_container_name)
    except DatabaseRequestError:
        raise DatabaseRequestError(f"DB missing plate information: "
                                   f"Wellplate '{plate_container_name}', "
                                   f"ref {container_nickname}")

    return assay, plate_doc, method_kwargs


def _lookup_spark_data_from_kwargs(manual_mode):

    method_kwargs = manual_mode.get('method_specs', {})
    plate_doc = manual_mode.get('plate_doc', None)

    if plate_doc is None:
        raise ValueError("plate_doc must be defined in manual mode")

    assay = method_kwargs.get('assay', 'assay-not-specified')

    return assay, plate_doc, method_kwargs
