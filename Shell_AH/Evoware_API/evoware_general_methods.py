# -*- coding: utf-8 -*-
"""
Created on Mon, July 26th 2021

Short methods that are generally useful for EVOWARE execution, these are mostly convenience functions based on the
DLL methods/functions that are accessible and detailed in evoware_api_support

@author: Brent Koscher
"""

import time


def status_lookup(status_code):
    status_table = {
        1:{'8': 'INITIALIZING', '2': 'LOADING', '1': 'NOINTRUMENTS', '4': 'NOTINITIALIZED', '0': 'UNKNOWN'},
        2:{'1': 'INITIALIZED', '4': 'SHUTDOWN', '2': 'SHUTTINGDOWN', '8': 'UNLOADING'},
        3:{'4': 'PAUSED', '2': 'PAUSEREQUESTED', '8': 'RESOURCEMISSING', '1': 'RUNNING'},
        4:{'1': 'DEADLOCK', '2': 'EXECUTIONERROR', '4': 'TIMEVIOLATION'},
        5:{'4': 'ABORTED', '2': 'BUSY', '1': 'IDLE', '8': 'STOPPED'},
        6:{'2': 'ERROR', '8': 'LOGON_ERROR', '1': 'PIPETTING', '4': 'SIMULATION'},
        7:{'1': 'CONNECTION_ERROR'}}
    status_list = []
    for i in range(1, len(status_code)):
        if status_code[-i] == 'x':
            break
        elif status_code[-i] == '0':
            continue
        else:
            status_list.append(status_table[i][status_code[-i]])
    return status_list


def initialize_evoware(evoware):
    evoware_status = check_evoware_status(evoware)
    print(evoware_status)
    time.sleep(2)
    if 'NOTINITIALIZED' in evoware_status:
        # First we need to load a worktable into EVOWARE
        empty_script_path = r'C:\Users\admin\PycharmProjects\AMD_Control_Platform\venv\Shell_AH\Evoware_API\process_scripts\empty_script.esc'
        evoware.prepare_script(empty_script_path)
        evoware.initialize()
    evoware_status = check_evoware_status(evoware)
    time.sleep(1)
    if 'INITIALIZED' in evoware_status:
        return ['INITIALIZED', '']
    else:
        return ['ERROR', evoware_status]


def check_evoware_status(evoware):
    evoware_status_code = evoware.get_status()
    evoware_status = status_lookup(str(hex(evoware_status_code)))
    return evoware_status


def check_evoware_systeminfo(evoware):
    systeminfo = evoware.get_system_info()
    if systeminfo[1] == False:
        mode = 'REAL'
    else:
        mode = 'SIMULATION'
    return mode


def execute_evoware_script(script_path, evoware):
    evoware_status = check_evoware_status(evoware)
    if 'INITIALIZED' in evoware_status:# and 'IDLE' in evoware_status:
        print('Ready to execute schedule!')
        script_id = evoware.start_script(script_path)
        time.sleep(5)
        return script_id
    else:
        print('Not ready to execute schedule!')
        return 'ERROR'