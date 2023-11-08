""" The callable functionality of the Mater-Controler System
MAIN METHOD FOR THE MC
@author: Ben C
"""
import time

import robot_arm as ra
import mcn_status as mcs
import functionals
from pprint import pformat
from datetime import datetime

from database_constants import *
from mcn_logging_manager import system_log
import database_interface as dbi
import data_repository_interface as dri
from custom_exceptions import DatabaseRequestError


def get_functions():
    """
    Called to populate the command_list of the System Object
    Think of this like a *.h method from the age of C
    :return: A dictionary of keywords : function handles
    """
    function_list = {
        'Initialize': initialize,
        'soft_wait': soft_wait,
        'await_queue_step': await_queue_step,
        'move_wellplate': move_wellplate,
        'Disable_power': power_stop,
        'load_pfa_film': load_pfa_film,
        'unload_pfa_film': unload_pfa_film,
        'acquire_robotic_arm': acquire_robotic_arm,
        'release_robotic_arm': release_robotic_arm,
    }
    return function_list

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


@functionals.standard_exception
def initialize(*_, **kwargs):
    _, internals, child_com = functionals.unpack_kwargs(kwargs, d_msg="")
    ra_socket = internals.get([mcs.S_AUXILIS, "MC", "ra_sock"], None)
    with child_com as cc:
        cc['Ra'] = ra.RoboticArmController(child_com, ra_socket)
        cc['Ra'].initialize()
        cc['is_initialized'] = True

    system_log.info(f"Testing connection to Databases: "
                    f"(Platform={dbi.test_connection()}, Results={dri.test_connection()})\n\t"
                    f"If you need to change DB settings, "
                    f"use the '__.Update_Database_Settings' command with no (kw)args")

    return mcs.RetObj.complete("Robotic Arm given controller")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


@functionals.standard_exception
def soft_wait(*_, **kwargs):
    """ Generates a scheduled wait for a given period of time

    Consumes a given amount of time without consuming the resource.  The actual time between operations will be at least
    the time given.  If the next operation needs to be executed upon waiting the given time, use the schedule_time
    attribute of the scheduler on the previous operation.

    (X)(soft_wait)(Y) := At least soft_wait time required between X and Y.

    (X: is_paired w/ schedule_time)(Y) := (Almost) Exactly schedule_time between X and Y.

    :param _: Consumes args from calling function
    :param kwargs: Given kwargs from calling function
    :keyword time_est: The time to wait (number or numerical string)
    :return: A return object
    """
    q_name, step_num, _, _ = functionals.unpack_db_kwargs(kwargs, d_pid="", d_a="")
    wait_time = kwargs.get(QOP_TIME, None)
    try:
        doc, _, _ = dbi.query_document("MC", COL_QUEUE, Q_NAME, q_name)
    except DatabaseRequestError:
        pass
    else:
        wait_time = doc[Q_OPERATIONS_LIST][step_num].get(QOP_TIME, wait_time)
        # Potential thing to do would be look at previous step's END_TIME and adjust the wait time accordingly
        #   Maybe have this feature if needed and only when a detail requests the time adjustment to be made.
    if wait_time is None:
        return mcs.RetObj.incomplete("Pr", mcs.V_PROBLEM, "Failed to find wait time in DB")
    system_log.info(f"Initiating soft wait on {q_name} for {wait_time//60} minutes")
    time.sleep(float(wait_time))
    # try:
    #     dbi.update_mongo_queue('SP', q_name, step_num, DB_YES)
    # except DatabaseRequestError:
    #     system_log.exception(f"Failed to update DB for ({q_name}, {step_num})")
    #     return mcs.RetObj.incomplete('MC', mcs.V_PROBLEM, "Failed to update DB with completion")
    return mcs.RetObj.complete(f"Waited {wait_time} seconds; "
                               f"the operation following can now be run at the scheduler's leisure")


@functionals.standard_exception
def await_queue_step(*_, **kwargs):
    """
    A method to run from (agent: MC.__) that waits for a specified queue's step to complete.
    Give it a very long time_est field.

    :param _: Consumes args from calling function
    :param kwargs: Given kwargs from calling function
    :return: A return object
    """
    q_name, _, _, _ = functionals.unpack_db_kwargs(kwargs, d_sn="", d_pid="", d_a="")
    details = kwargs.get(QOP_DETAILS, {})
    try:
        await_queue = details['await_queue']
    except KeyError:
        return mcs.RetObj.incomplete('_U', mcs.V_FATAL,
                                     f"await_queue_step in '{q_name}' called without a queue to await")
    try:
        await_num = str(int(details['step_num']))
    except KeyError:
        return mcs.RetObj.incomplete('_U', mcs.V_FATAL,
                                     f"await_queue_step in '{q_name} called without a step number to await'")
    await_operation = details.get('operation', None)
    await_interval = details.get('interval', 5)  # minutes
    await_timeout = details.get('timeout', 5*24*60)  # minutes

    timer = datetime.now()
    while functionals.realtime_wait(timer, minutes=await_timeout):
        try:
            doc, _, _ = dbi.query_document("MC", COL_QUEUE, Q_NAME, await_queue)
        except DatabaseRequestError:
            time.sleep(1.0 + 60*await_interval/4)
            continue

        try:
            operations_list = doc[Q_OPERATIONS_LIST]
            if not operations_list:
                raise ValueError
        except (KeyError, ValueError):
            return mcs.RetObj.incomplete('_U', mcs.V_FATAL,
                                         f"await_queue_step in '{q_name}' missing a queue document ({await_queue}) "
                                         f"with a real operations list")

        if await_operation is not None:
            await_num = functionals.find_shifted_step_number(operations_list, await_num, await_operation, await_queue)

        awaited_step = operations_list[await_num]
        if awaited_step[QOP_COMPLETED] == DB_YES:
            return mcs.RetObj.complete(f"The awaited method is complete")

        time.sleep(60*await_interval)

    return mcs.RetObj.incomplete('MC', mcs.V_PROBLEM, "Timeout error on await_queue_step")


@functionals.standard_exception
def acquire_robotic_arm(*_, **kwargs):
    """ Requests the robotic arm to Lock

    :param kwargs: Given kwargs from calling function
    :keyword blocking: the Lock.acquire(blocking=...) keyword, default: True
    :keyword timeout: How long to wait, the Lock.acquire(timeout=...) keyword, default: -1
    :return: Complete or Incomplete RetObj given the success of the call
    """
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    lock_blocking = kwargs.get('blocking', True)
    lock_timeout = kwargs.get('timeout', -1)
    with child_com as cc:
        rac: ra.RoboticArmController = cc["Ra"]
    system_log.info(f"Requesting Ra Lock with args: blocking={lock_blocking} & timeout={lock_timeout}")
    if rac.acquire(lock_blocking, lock_timeout):
        system_log.info("Ra Lock acquired")
        return mcs.RetObj.complete("Robotic Arm Lock Acquired")
    else:
        system_log.info("Ra Lock not acquired")
        return mcs.RetObj.incomplete("Ra", mcs.V_PROBLEM, "Robotic Arm Lock NOT Acquired")


@functionals.standard_exception
def release_robotic_arm(*_, **kwargs):
    """ Requests the robotic arm to Unlock

    :param kwargs: Given kwargs from calling function
    :return: A Complete RetObj (always)
    """
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    with child_com as cc:
        rac: ra.RoboticArmController = cc["Ra"]
    try:
        rac.release()
    except RuntimeError:
        pass
    return mcs.RetObj.complete("Robotic Arm Lock Released")


@functionals.standard_exception
def move_wellplate(*_, **kwargs):
    """ Calls Robotic Arm (Ra) to move a plate

    By Specifying [Q_NAME, DBQ_STEP_NUM], it will run from a queue document

    By Specifying ['override'] with override={'container_name': str, 'destination': str}, it will run from override

    :param kwargs: Requires keys: [CHILD_COM] and (['override'] xor [Q_NAME, DBQ_STEP_NUM])
    :return: Return Object
    """
    override = kwargs.get('override', None)
    queue_id = kwargs.get(Q_NAME, None)
    operation_id = kwargs.get(DBQ_STEP_NUM, None)

    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")

    # This check is redundant as ra.move_plate() will raise a ValueError anyway
    db_mode = queue_id and operation_id
    if (override or db_mode) and not (override and db_mode):
        pass
    else:
        system_log.info(f"Improperly formatted kwargs for move_wellplate():\n{pformat(kwargs)}")
        return mcs.RetObj.incomplete("Ra", mcs.V_PROBLEM, f"Improperly formatted kwargs for move_wellplate")

    with child_com as cc:
        rac: ra.RoboticArmController = cc["Ra"]
    with rac as r:
        return r.move_plate(q_id=queue_id, q_step=operation_id, override=override)


@functionals.standard_exception
def load_pfa_film(*_, **kwargs):
    """ Loads a PFA film into the thermoreactor

    :param kwargs: Given kwargs from calling function
    :return: A Return Object
    """
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    with child_com as cc:
        rac: ra.RoboticArmController = cc["Ra"]
    with rac as r:
        return r.move_pfa_film(load_clean=True)


@functionals.standard_exception
def unload_pfa_film(*_, **kwargs):
    """ Unloads a PFA film from the thermoreactor

    :param kwargs: Given kwargs from calling function
    :return: A Return Object
    """
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    with child_com as cc:
        rac: ra.RoboticArmController = cc["Ra"]
    with rac as r:
        return r.move_pfa_film(load_clean=False)


@functionals.standard_exception
def power_stop(*args, **kwargs):
    """ Calls stop method for robotic arm

    :param args: Args passed to stop() method
    :param kwargs: Kwargs passed to stop() method
    :return: A Return Object
    """
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    with child_com as cc:
        rac: ra.RoboticArmController = cc["Ra"]
    with rac as r:
        return r.stop(*args, **kwargs)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


if __name__ == '__main__':
    from threading import Thread
    import system

    test_system1 = system.System('MC')
    c1 = Thread(target=test_system1.run)
    c1.start()
    c1.join()
