# MAIN METHOD FOR LC
# imports the necessary modules for the Shimadzu wrapper
# (which controls the autosampler, HPLC-MS, and fraction collector)
from constants import *
import functionals
import time
import Shimadzu_API
from pprint import pformat
import mcn_status as mcs
from mcn_logging_manager import system_log
import database_interface as dbi
import data_repository_interface as dri
import custom_exceptions as cexc


# Define Mid-Level commands here
def get_functions():
    """
    Called to populate the command_list of the System Object
    Think of this like a *.h method from the age of C
    :return: A dictionary of keywords : function handles
    """
    function_list = {
        'Initialize': initialize,
        "Status": status,
        "Stop": lc_stop,
        "run_analytical_batch": run_analytical_batch,
        'run_semiprep_batch': run_semiprep_batch,
        'run_uplc_batch': run_uplc_batch,
    }
    return function_list

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


def initialize(*_, **kwargs):
    lcs = False
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    with child_com as cc:  # noqa
        try:
            lcs = bool(cc['shimadzu_wrapper'])
        except KeyError:
            pass

    system_log.info(f"Testing connection to Databases: "
                    f"(Platform={dbi.test_connection()}, Results={dri.test_connection()})\n\t"
                    f"If you need to change DB settings, "
                    f"use the '__.Update_Database_Settings' command with no (kw)args")

    if lcs:
        cc['is_initialized'] = True
        return mcs.RetObj.complete("Sw given data pipe")
    else:
        system_log.warning("Sw failed to initialize")
        return mcs.RetObj.incomplete('LC', mcs.V_FATAL, "Sw failed to initialize")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


@functionals.standard_exception
def status(*_, **kwargs):
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    with child_com as cc:  # noqa
        lc_object = cc['shimadzu_wrapper']
    stat, details = lc_object.status()
    return mcs.RetObj.complete([stat, details])


@functionals.standard_exception
def lc_stop(*_, **kwargs):
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    with child_com as cc:  # noqa
        lc_object = cc['shimadzu_wrapper']
    completion, info = lc_object.lc_stop()
    if completion:
        return mcs.RetObj.complete(info)
    else:
        system_log.info(f"LC failed to stop:\n"
                        f"{pformat(info)}")
        return mcs.RetObj.incomplete('LC', mcs.V_FATAL, "LC failed to stop")
    

@functionals.standard_exception
def run_analytical_batch(*_, **kwargs):
    """ Starts an analytical batch and waits for it to complete

    :param kwargs:
    :return:
    """
    return run_batch(mode="analytical", **kwargs)


@functionals.standard_exception
def run_uplc_batch(*_, **kwargs):
    """ Starts an analytical batch and waits for it to complete

    :param kwargs:
    :return:
    """
    return run_batch(mode="uplc", **kwargs)


@functionals.standard_exception
def run_semiprep_batch(*_, **kwargs):
    """ Starts a semiprep batch and waits for it to complete

    :param kwargs:
    :return:
    """
    return run_batch(mode="semiprep", **kwargs)


@functionals.standard_exception
def run_optimization_batch(*_, **kwargs):
    """ Starts an optimization batch and waits for it to complete

    :param kwargs:
    :return:
    """
    return run_batch(mode="optimization", **kwargs)


def run_batch(*_, mode: str, **kwargs):
    """ Starts an analytical batch and waits for it to complete

    :param mode: analytical or semiprep
    :param kwargs:
    :return:
    """
    document_queue_name, operation_number, _, _ = functionals.unpack_db_kwargs(kwargs, d_pid="", d_a="")
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")

    if mode == "analytical":
        required_column = 'analytical'
        user_prompt = "Please Install Analytical Column"
        default_prompt = required_column
        ret_obj_value = "LC analytical batch finished"
    elif mode == "semiprep":
        required_column = 'semiprep'
        user_prompt = "Please Install Semi-Prep Column"
        default_prompt = required_column
        ret_obj_value = "LC semiprep batch finished"
    elif mode == "uplc":
        required_column = 'uplc'
        user_prompt = "Please Install UPLC Column"
        default_prompt = required_column
        ret_obj_value = "LC UPLC iteration finished"
    else:
        raise ValueError(f"HPLC can only be run in analytical, optimization, semiprep, or UPLC mode, not {mode}")

    with child_com as cc:  # noqa
        lc_object = cc['shimadzu_wrapper']
    lc_object.lock.acquire()
    try:
        status_code, details = lc_object.status()
        if status_code != 'Ready':
            system_log.info(f"LC not ready for batch ({status_code}):\n"
                            f"{pformat(details)}")
            return mcs.RetObj.incomplete('LC', status_code, "LC not ready for batch")

        #while lc_object.lc_column != required_column:
        #    lc_object.lc_column = functionals.quick_gui(title='Column Check',
        #                                                dialog=user_prompt,
        #                                                has_entry=False, ret_if_ok=default_prompt)
        #    system_log.info(f"HPLC column set to '{lc_object.lc_column}'")

        if mode == "semiprep":
            batch_result, batch_msg = lc_object.run_semiprep_batch(document_queue_name, operation_number)
        elif mode == "analytical":
            batch_result, batch_msg = lc_object.run_analytic_batch(document_queue_name, operation_number)
        elif mode == "uplc":
            batch_result, batch_msg = lc_object.run_analytic_batch(document_queue_name, operation_number)
        elif mode == "optimization":
            batch_result, batch_msg = lc_object.run_optimization_batch(document_queue_name, operation_number)
        else:
            # should be unreachable
            batch_result = False
            batch_msg = f"mode unrecognized: '{mode}'"

        if batch_result:
            while lc_object.running == (False, 'Initialized'):
                time.sleep(5)
            time.sleep(30)
            while lc_object.running[0]:
                time.sleep(5)
            if lc_object.running[1] == 'Batch Ended':
                return mcs.RetObj.complete([ret_obj_value, batch_msg])
            else:
                system_log.warning(f"LC error state (LC not running but batch not ended):\n"
                                   f"{pformat(lc_object.running)}")
                # Future improvement, on PMax error, flush column with THF or ACN, return Problem first then busy
                return mcs.RetObj.incomplete('LC', mcs.V_FATAL, "LC not running but Batch not ended")
        else:
            system_log.info(f"LC Batch failed: {pformat(batch_msg)}")
            if batch_msg[0:2] == "P_":
                return mcs.RetObj.incomplete('LC', mcs.V_PROBLEM, "LC: " + batch_msg)
            elif batch_msg[0:2] == "F_":
                return mcs.RetObj.incomplete('LC', mcs.V_FATAL, "LC: " + batch_msg)
            elif batch_msg[0:2] == "B_":
                return mcs.RetObj.incomplete('LC', mcs.V_BUSY, "LC: " + batch_msg)

    finally:
        lc_object.lock.release()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


if __name__ == '__main__':
    from threading import Thread
    # from custom_classes import Nop
    import system

    lc_system = system.System('LC')
    try:
        with lc_system.data_pipes[CHILD_COM] as cc:
            cc['shimadzu_wrapper'] = Shimadzu_API.ShimadzuController(lc_system.data_pipes)
            safety_harness = cc['shimadzu_wrapper']
        c3 = Thread(target=lc_system.run)
        c3.start()
        c3.join()
    finally:
        try:
            safety_harness.disconnect()
        except NameError:
            print("System's Clear got to it first")
