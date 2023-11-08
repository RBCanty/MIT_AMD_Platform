# MAIN METHOD FOR AH
# Note: Liquid handler known to crash when this is run from PyCharm

import glob
import queue
import sys
import threading
import time

project_dirs = glob.glob(r'C:\Users\admin\PycharmProjects\AMD_Control_Platform\venv\*\\')
[sys.path.append(project_dir) for project_dir in project_dirs if project_dir not in sys.path]

if r'C:\Users\admin\PycharmProjects\AMD_Control_Platform\venv\Shell_AH\Spark_API\\' not in sys.path:
    sys.path.append(r'C:\Users\admin\PycharmProjects\AMD_Control_Platform\venv\Shell_AH\Spark_API\\')

if r'C:\Users\admin\PycharmProjects\AMD_Control_Platform\venv\Ancillaries\Meta\\' not in sys.path:
    sys.path.append(r'C:\Users\admin\PycharmProjects\AMD_Control_Platform\venv\Ancillaries\Meta\\')

import synchronization as sync
import functionals
import mcn_status as mcs
from constants import *
from pprint import pprint, pformat
import database_interface as dbi
import data_repository_interface as dri
import operations as oprtn
from mcn_logging_manager import system_log
import aceso
import traceback
import numpy as np
from datetime import datetime
from copy import deepcopy


try:
    from Evoware_API.evoware_general_methods import status_lookup, initialize_evoware, check_evoware_status, check_evoware_systeminfo, status_lookup, execute_evoware_script
    from Evoware_API import evoware_api_support
    from Evoware_API.reaction_plate_preparation_function import reaction_plate_preparation
    from Evoware_API.transfer_wellplate_function import transfer_wellplate
    from Evoware_API.start_stop_heater_shaker_function import start_stop_heater_shaker
    from Evoware_API.filter_plate_function import filter_plate_sequence
    from Evoware_API.copy_wellplate_function import copy_wellplate
    from Evoware_API.liquid_liquid_extraction import liquid_liquid_extraction
    from Evoware_API.prepare_standard_wells import prepare_standard_wells
    from Evoware_API.prepare_characterization_plate import characterization_plate_preparation
    from Evoware_API.detect_liquid_level_function import detect_liquid_level
except ModuleNotFoundError as e:
    print("WARNING! Evoware_API could not be imported")
    print(traceback.format_exc())


# Define Mid-Level commands here
def get_functions():
    """
    Called to populate the command_list of the System Object
    Think of this like a *.h method from the age of C
    :return: A dictionary of keywords : function handles
    """

    # Based on Brent's code it seems like Lh just has one function since it manages the 'prepare_wellplate' vs
    # 'transfer_wellplate' vs etc stuff already.  However, the way running a method from the DB works, the
    # operation name stuff ('prepare_wellplate', etc.) is stored in kwargs[FUNC].  So, we need to do a thing
    # where all operation names for Lh ('prepare_wellplate', 'transfer_wellplate', etc) all map to the same function.

    function_list = {
        "Initialize": initialize,
        'prepare_wellplate': evoware_method_selector,  # See comment on line 27
        'transfer_wellplate': evoware_method_selector,        # ...
        'start_stop_heater_shaker': evoware_method_selector,  # ...
        'filter_wellplate': evoware_method_selector,          # ...
        'prepare_standard_wells': evoware_method_selector,
        'liquid_liquid_extraction': evoware_method_selector, 
        'copy_wellplate': evoware_method_selector,     # See comment on line 27
        'check_evoware_status': check_evoware_ready_status,
        'prepare_characterization_plate': evoware_method_selector,
        'detect_liquid_level': evoware_method_selector,
        'Lh_monitor': lh_monitor,
    }
    return function_list


def initialize(*_, **kwargs):
    _, internals, child_com = functionals.unpack_kwargs(kwargs, d_msg="")
    with child_com as cc:
        try:
            cc['Lh'] = queue.Queue()
            cc['Lh'].put(['Startup'])
            lhs = True
        except KeyError:
            lhs = False

    system_log.info(f"Testing connection to Databases: "
                    f"(Platform={dbi.test_connection()}, Results={dri.test_connection()})\n\t"
                    f"If you need to change DB settings, "
                    f"use the '__.Update_Database_Settings' command with no (kw)args")

    if lhs:
        with child_com as cc:
            cc['is_initialized'] = True
        return mcs.RetObj.complete("Lh and Pr given data pipes")
    else:
        exc_msg = f"Failed to give data pipes to children--(Lh:{lhs})"
        system_log.info(exc_msg)
        return mcs.RetObj.incomplete('AH', mcs.V_FATAL, exc_msg)


def lh_monitor(pipes):
    ret_obj = check_evoware_ready_status(**{CHILD_COM: pipes[CHILD_COM]})
    if ret_obj.level == mcs.V_FATAL:
        sync.synchronized_fault_add(pipes, mcs.Fault.from_retobj(ret_obj))
    elif ret_obj.level == mcs.V_BUSY:
        internals = pipes[INTERNAL]
        busy_confirm = internals.is_busy(subsystem='Lh')
        if not busy_confirm:
            system_log.warning("Lh reports it is busy but is not associated with a checkpoint")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


@functionals.standard_exception
def check_evoware_ready_status(*_, **kwargs) -> mcs.RetObj:
    """ Reads the Liquid Handler status queue

    :param kwargs: Used for kwargs[CHILD_COM]['Lh'] tp get liquid handler status queue
    :return: standard return object (idle and busy - True)(error - False, to create fault)
    """
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    with child_com as cc:
        evoware_status_queue = cc['Lh']
    if 'IDLE' in list(evoware_status_queue.queue):
        return mcs.RetObj.complete("tecan idle")
    elif 'BUSY' in list(evoware_status_queue.queue):
        system_log.info("Liquid handler reporting Busy")
        return mcs.RetObj.incomplete("Lh", mcs.V_BUSY, "tecan busy")
    elif 'ERROR' in list(evoware_status_queue.queue):
        system_log.warning("Liquid handler reporting error state")
        return mcs.RetObj.incomplete('Lh', mcs.V_FATAL, "tecan error msg")


@functionals.standard_exception
def evoware_method_selector(*_, **kwargs) -> mcs.RetObj:
    """ Starts a method then watches it, signals MCN when completed or faulted

    :param kwargs: Used for queue name, step number, function name, and Lh status queue
    :return: Standard return object
    """
    document_queue_name, operation_number, _, _ = functionals.unpack_db_kwargs(kwargs, d_pid="", d_a="")
    evoware_method_name = kwargs[oprtn.FUNC]
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    with child_com as cc:
        evoware_queue = cc['Lh']
    start_return = evoware_thread_spawn(document_queue_name, operation_number, evoware_method_name, evoware_queue)
    if start_return[0] != 'EVOWARE Started':
        if start_return[0] == 'EVOWARE already busy':
            system_log.info(f"Liquid handler prompted for ({document_queue_name}, #{operation_number}) but was busy")
            return mcs.RetObj.incomplete('Lh', mcs.V_BUSY, 'Tecan already busy')
        elif start_return[0] == 'EVOWARE ERROR':
            system_log.warning(f"Liquid handler reported error state when "
                               f"prompted for ({document_queue_name}, #{operation_number})")
            return mcs.RetObj.incomplete('Lh', mcs.V_FATAL, 'Tecan error')
        else:
            system_log.warning(f"Liquid handler reporting unkown state "
                               f"when prompted for ({document_queue_name}, #{operation_number}):\n"
                               f"{pformat(start_return)}")
            return mcs.RetObj.incomplete('Lh', mcs.V_FATAL, "Tecan unhandled error")
    while True:
        with child_com as cc:
            evoware_queue = cc['Lh'].queue[0]

        error_reporter = evoware_api_support.evo_error_watcher()
        error_return = error_reporter.check_for_errors()
        if error_return[0] == 'Error':
            system_log.warning(f"Liquid handler reporting error state:\n"
                               f"{pformat(error_return)}")
            return mcs.RetObj.incomplete('Lh', mcs.V_FATAL, "Tecan error return")
        if 'IDLE' in evoware_queue:
            return mcs.RetObj.complete("tecan success")
        elif 'ERROR' in evoware_queue:
            if 'No tips available' in evoware_queue:
                system_log.info("Liquid handler did not have enough tips, it has made a request to Ss for more tips")
                with child_com as cc:
                    evoware_queue = cc['Lh']
                evoware_queue.queue.popleft()
                evoware_queue.put(['IDLE', ''])
                return mcs.RetObj.incomplete('Lh', mcs.V_BUSY, "Insufficient tips")
            elif any('No suitable labware on the liquid handler' in queue_item for queue_item in evoware_queue):
                system_log.info("Liquid handler did not have enough labware, it has tried to make a request")
                return_string = deepcopy(evoware_queue[1])
                lpx_flag = return_string.split('>')[1].strip()
                with child_com as cc:
                    evoware_queue = cc['Lh']
                evoware_queue.queue.popleft()
                evoware_queue.put(['IDLE', ''])
                if lpx_flag == '1':
                    return mcs.RetObj.incomplete('Lh', mcs.V_BUSY, return_string)
                else:
                    return mcs.RetObj.incomplete("Lh", mcs.V_PROBLEM, return_string)
            elif 'Soft wait needed' in evoware_queue:
                system_log.info('Evoware needs a soft wait for an operation')
                with child_com as cc:
                    evoware_queue = cc['Lh']
                evoware_queue.queue.popleft()
                evoware_queue.put(['IDLE', ''])
                return mcs.RetObj.incomplete('Lh', mcs.V_BUSY, "Soft wait needed")
            elif 'Problem:' in evoware_queue[1]:
                return_string = deepcopy(evoware_queue[1])
                return mcs.RetObj.incomplete("Lh", mcs.V_PROBLEM, return_string)
            else:
                system_log.warning(f"Liquid handler reporting error state:\n"
                                   f"{pformat(evoware_queue)}")
                return mcs.RetObj.incomplete('Lh', mcs.V_FATAL, "tecan error msg")
        else:
            time.sleep(2)


def evoware_thread_func(name, document_queue_name, operation_number, evoware_method_name, evoware_queue):
    """ Parses message to determine which EVOWARE method to run, generates a suitable script to execute with
    the platform.

    :param name: Not used
    :param document_queue_name: Name of the queue document from thr DB
    :param operation_number: The operation number (DBQ_STEP_NUM)
    :param evoware_method_name: Name of the evoware method being run
    :param evoware_queue: Queue used to store EVOWARE status
    :return: on success it returns ['Success', <script_filepath>] and on error it returns ['Error', 'Error Message']
    """
    if evoware_method_name == 'prepare_wellplate':
        method_return = reaction_plate_preparation(document_queue_name, operation_number)
    elif evoware_method_name == 'transfer_wellplate':
        method_return = transfer_wellplate(document_queue_name, operation_number)
    elif evoware_method_name == 'start_stop_heater_shaker':
        method_return = start_stop_heater_shaker(document_queue_name, operation_number)
    elif evoware_method_name == 'filter_wellplate':
        method_return = filter_plate_sequence(document_queue_name, operation_number)
    elif evoware_method_name == 'copy_wellplate':
        method_return = copy_wellplate(document_queue_name, operation_number)
    elif evoware_method_name == 'prepare_standard_wells':
        method_return = prepare_standard_wells(document_queue_name, operation_number)
    elif evoware_method_name == 'liquid_liquid_extraction':
        method_return = liquid_liquid_extraction(document_queue_name, operation_number)
    elif evoware_method_name == 'prepare_characterization_plate':
        method_return = characterization_plate_preparation(document_queue_name, operation_number)
    elif evoware_method_name == 'detect_liquid_level':
        method_return = detect_liquid_level(document_queue_name, operation_number)
    else:
        evoware_queue.queue.popleft()
        evoware_queue.put(['ERROR', 'Problem: Method not implemented'])
        method_return = ['ERROR', 'Problem: Method not implemented']
    pprint(method_return)
    e_stop_flag = 0
    # If the message was successfully parsed then we can start up the EVOWARE API and execute the script
    
    if 'Nothing ' in method_return[1]:
        evoware_queue.queue.popleft()
        evoware_queue.put(['IDLE', ''])
    elif method_return[0] == 'Success':
        script_path = method_return[1]  # Absolute path
        # script_path = r'C:\Users\admin\PycharmProjects\AMD_Control_Platform\venv\Shell_AH\Evoware_API\process_scripts\testing_script.esc'
        evoware = evoware_api_support.EVO()
        login_return = evoware.logon()
        initialize_return = initialize_evoware(evoware)
        # print(initialize_return)
        if 'INITIALIZED' in initialize_return:
            evoware_queue.queue.popleft()
            evoware_queue.put(['INITIALIZED', ''])
            execute_return = execute_evoware_script(script_path, evoware)
            time.sleep(2)
            if execute_return != 'ERROR':
                evoware_queue.queue.popleft()
                evoware_queue.put(['BUSY', execute_return])
                while True:
                    current_evoware_status = check_evoware_status(evoware)
                    if 'IDLE' in current_evoware_status:
                        break
                    if 'EMERGENCYSTOP' in evoware_queue.queue[0]:
                        e_stop_flag = 1
                        evoware_queue.queue.popleft()
                        evoware_queue.put(['ERROR', 'Fatal: EMERGENCY STOP CALLED'])
                        break
                    time.sleep(1)
            else:
                evoware_queue.queue.popleft()
                evoware_queue.put(['ERROR', 'Fatal: Script Execution Error'])
        else:
            evoware_queue.queue.popleft()
            evoware_queue.put(['ERROR', 'Fatal: Initialization Error', initialize_return])
        logoff_return = evoware.logoff()
        if e_stop_flag == 1:
            return
        if evoware_method_name == 'detect_liquid_level':
            try:
                current_time = datetime.now()
                relative_filepath = r'.\Evoware_API\process_scripts\liquid_detection\liquid_detection_results.txt'
                with open(relative_filepath, 'r') as infile:
                    header = infile.readline()
                    value = infile.readline()
                    numerical_value = float(value.strip())
                queue_document, _, _ = dbi.query_document('AH', 'queue', 'queue_name', document_queue_name)
                initial_queue_document = deepcopy(queue_document)
                target_operation = queue_document['operations'][operation_number]
                insert_soft_wait = False
                if 'volume_history' not in target_operation['details'].keys():
                    target_operation['details']['volume_history'] = []
                    if numerical_value > 10:
                        insert_soft_wait = True
                target_operation['details']['volume_history'].append([current_time.strftime('%m/%d/%Y %H:%M:%S'),
                                                                      numerical_value])
                doc, err_details, resp_code = dbi.update_document('AH', 'queue', initial_queue_document, queue_document)
                if type(doc) is not str or doc != 'Success':
                    evoware_queue.queue.popleft()
                    evoware_queue.put(['ERROR', 'Problem: Issue updating liquid level detection'])
                if numerical_value < 10:
                    time_to_wait = 0
                elif numerical_value >= 10 and len(target_operation['details']['volume_history']) > 2:
                    volume_history = target_operation['details']['volume_history']
                    t_zero = datetime.strptime(volume_history[0][0], '%m/%d/%Y %H:%M:%S')
                    times = [(datetime.strptime(t_point[0], '%m/%d/%Y %H:%M:%S') - t_zero).total_seconds() for t_point in volume_history]
                    volumes = [t_point[1] for t_point in volume_history]
                    curve_fit = np.polyfit(times, np.log(volumes), 1)
                    target_volume = 25
                    forecast_time = (np.log(target_volume) - curve_fit[1]) / curve_fit[0]
                    time_to_wait = min([max([1800, forecast_time]), 3600 * 4])
                else:
                    time_to_wait = 3600

                if time_to_wait > 0 and insert_soft_wait:
                    new_steps = [{'operation': 'soft_wait',
                                  'container': '',
                                  'details': {},
                                  'completed': 'no',
                                  'time_est': int(time_to_wait),
                                  'agent': 'MC.__',
                                  'start_time': None,
                                  'end_time': None}]
                    dbi.insert_queue_steps(document_queue_name, int(operation_number)-1, new_steps)
                    evoware_queue.queue.popleft()
                    evoware_queue.put(['ERROR', 'Soft wait needed'])
                elif time_to_wait > 0 and not insert_soft_wait:
                    dbi.update_field('AH', 'operations.%s.completed' % str(int(operation_number)-1),
                                     document_queue_name, 'yes', 'no')
                    dbi.mark_time(document_queue_name,
                                  str(int(operation_number)-1),
                                  target_operation.get(dbi.QOP_OPERATION, None),
                                  reset=True,
                                  logger=system_log)
                    evoware_queue.queue.popleft()
                    evoware_queue.put(['ERROR', 'Soft wait needed'])
                else:
                    evoware_queue.queue.popleft()
                    evoware_queue.put(['IDLE', ''])
            except:
                evoware_queue.queue.popleft()
                evoware_queue.put(['ERROR', traceback.format_exc()])
        elif evoware_method_name == 'copy_wellplate':
            try:
            
                for return_statement in method_return[2]:
                    if 'Checking time for aliquot' in return_statement:
                        split_statement = return_statement.split('>')
                        wellplate_document, _, _ = dbi.query_document('AH', 'wellplates', 'container_name', split_statement[2])
                        original_document = deepcopy(wellplate_document)
                        
                        with open(r'C:\Users\admin\PycharmProjects\AMD_Control_Platform\venv\Shell_AH\Evoware_API\process_scripts\datetime_output\current_datetime.txt', 'r') as infile:
                            datetime_string = infile.readline()
                        wellplate_document['aliquot_time_taken'] = datetime_string.strip()
                        doc, err_details, resp_code = dbi.update_document('AH', 'wellplates', original_document, wellplate_document)
                        if type(doc) is not str or doc != 'Success':
                            evoware_queue.queue.popleft()
                            evoware_queue.put(['ERROR', 'Problem: Issue updating liquid level detection'])
                        break
                evoware_queue.queue.popleft()
                evoware_queue.put(['IDLE', ''])
            except:
                evoware_queue.queue.popleft()
                evoware_queue.put(['ERROR', traceback.format_exc()])
        else:
            evoware_queue.queue.popleft()
            evoware_queue.put(['IDLE', ''])
        time.sleep(2)
    elif method_return[1] == 'No tips available':
        aceso.tip_request(document_queue_name, operation_number)
        evoware_queue.queue.popleft()
        evoware_queue.put(['ERROR', method_return[1]])

    elif 'No suitable labware on the liquid handler' in method_return[1]:
        lpx_flag = method_return[1].split('>')[1].strip()
        print('lpx_flag', lpx_flag)
        container_nickname = method_return[1].split('>')[3].strip()
        print(container_nickname)
        if lpx_flag == '1':
            return_bool, nickname, extra_return = aceso.resource_request(document_queue_name, operation_number, 'liquid_handler', False, container_nickname)
            print(return_bool)
            print(nickname)
            print(extra_return)
        evoware_queue.queue.popleft()
        evoware_queue.put(['ERROR', method_return[1]])

    else:
        evoware_queue.queue.popleft()
        evoware_queue.put(['ERROR', method_return[1]])


def evoware_thread_spawn(document_queue_name, operation_number, evoware_method_name, evoware_status_queue):
    """ Checks the status queue to see if EVOWARE is ready and if so (is ready/IDLE) we can spawn a thread
    to generate a schedule to send into the remote mode of EVOWARE

    :param document_queue_name: Name of the queue document from thr DB
    :param operation_number: The operation number (DBQ_STEP_NUM)
    :param evoware_method_name: Name of the evoware method being run
    :param evoware_status_queue: Queue used to store EVOWARE status
    :return: [string code for current status, if error or busy evoware_status_queue.queue]
    """
    evoware_thread_name = 'evoware_thread'
    if 'ERROR' in list(evoware_status_queue.queue):
        return ['EVOWARE ERROR', evoware_status_queue.queue]
    if 'BUSY' in list(evoware_status_queue.queue):
        return ['EVOWARE already busy', evoware_status_queue.queue]
    evoware_thread = threading.Thread(target=evoware_thread_func, name=evoware_thread_name,
                                      args=(evoware_thread_name, document_queue_name, operation_number,
                                            evoware_method_name, evoware_status_queue),
                                      daemon=True)
    evoware_status_queue.queue.popleft()
    evoware_status_queue.put(['Starting EVOWARE Thread', ''])
    evoware_thread.start()
    return ['EVOWARE Started']


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


if __name__ == '__main__':
    import system

    test_system2 = system.System('AH')
    c2 = threading.Thread(target=test_system2.run)
    c2.start()
    c2.join()
