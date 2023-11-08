# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 13:47:06 2021

Supporting functions for Evoware API for MCN module

@author: Brent Koscher
"""

from win32com.client import DispatchWithEvents
import pythoncom
import time
import sys
import re

if r'C:\ProgramData\Anaconda3\Lib\site-packages\win32\\' not in sys.path:
    sys.path.append(r'C:\ProgramData\Anaconda3\Lib\site-packages\win32\\')

import win32gui


class LoginInformation:
    username = "Admin"
    password = "admin3"


class EVO:
    def __init__(self):
        pythoncom.CoInitialize()
        self.system = DispatchWithEvents('EVOAPI.System', EVOSystemEvents)

    def cancel_process(self, processID):
        # self.system.CancelProcess(int ProcessID)
        return 'Not implemented'

    def execute_script_command(self, command):
        # self.system.ExecuteScriptCommand(string bstrCommand)
        return 'Not implemented'

    def get_device(self):
        # self.system.GetDeviceObject(object Index)
        return 'Not implemented'

    def get_device_count(self):
        number_of_devices = self.system.GetDeviceCount()
        return number_of_devices

    def get_lc_count(self):
        lc_count = self.system.GetLCCount()
        return lc_count

    def get_lc_info(self, lc_index):
        # self.system.GetLCInfo(int Index, out string lcName, out bool isDefault, out bool isCustomized)
        lc_info = self.system.GetLCInfo(lc_index)
        return lc_info

    def get_number_of_racks(self, script_id):
        number_of_racks = self.system.GetNumberOfRacks(script_id)
        return number_of_racks

    def get_process_status(self, process_id):
        process_status = self.system.GetProcessStatus(process_id)
        process_message = EVOConstants.PROCESS_STATUS[process_status]
        return [process_status, process_message]

    def get_process_status_ex(self, process_id):
        # process_status = self.system.GetProcessStatusEx(process_id)
        # process_message = EVOConstants.PROCESS_STATUS[process_status]
        # return [process_status, process_message]
        return 'Not implemented'

    def get_process_variable(self, process_id, variable_name, variable_scope):
        # process_variable = self.system.GetProcessVarable(int processID, string VariableName, SC_VARIABLEScope Scope)
        return 'Not implemented'

    def get_rack(self, script_id, index):
        # rack_details = self.system.GetRack(int index, out string Name, out string Label, out int Location, out int Grid, out int Site, out string CarrierName)
        return 'Not implemented'

    def get_script_status(self, script_id):
        returned_value = self.system.GetScriptStatus(script_id)
        return returned_value

    def get_script_status_ex(self, script_id):
        returned_value = self.system.GetScriptStatus(script_id)
        return returned_value

    def get_script_variable(self, script_id, variable_name, variable_scope):
        # process_variable = self.system.GetScriptVarable(int ScriptID, string VariableName)
        return 'Not implemented'

    def get_status(self):
        returned_value = self.system.GetStatus()
        return returned_value

    def get_sub_lc_count(self, lc_index):
        # sub_lc_count = self.system.GetSubLCCount(int index)
        return 'Not implemented'

    def get_sub_lc_info(self, lc_index, sub_lc_index):
        # sub_lc_info = self.system.GetSubLCInfo(int lcindex, int subLCIndex, out SC_TipType tipType, out bool isAllVolumes, out double minVolume, out double maxVolume)
        return 'Not implemented'

    def get_system_info(self):
        returned_value = self.system.GetSystemInfo()
        return returned_value

    def get_window_handle(self):
        # handles = self.system.GetWindowHandles(out int MainWinHandle, out int ScriptWinHandle, out int WorktableWinHandle, out int LogWinHandle)
        return 'Not implemented'

    def hide_gui(self, hide_val):
        # self.system.HideGUI(hide_val)
        return 'Not implemented'

    def initialize(self):
        self.system.Initialize()
        returned_value = self.system.GetStatus()
        return returned_value

    def logoff(self):
        self.system.Logoff()
        return 'Remote control turned off'

    def logon(self):
        # Logon(string UserName, string Password, int Plus, int Simulation)
        # We run EvoWare standard so Plus = 0
        self.system.Logon(LoginInformation.username, LoginInformation.password, 0, 0)
        while self.system.GetStatus() == EVOConstants.no_instruments_status:
            time.sleep(1)
        while self.system.GetStatus() == EVOConstants.loading_status:
            time.sleep(1)
        return 'Remote control turned on'

    def pause(self):
        self.system.Pause()
        return 'System Paused'

    def prepare_process(self, process_name):
        # self.system.PrepareProcess(process_name)
        return 'Not implemented'

    def prepare_script(self, script_name):
        # self.system.PrepareScript(string scriptName)
        return 'Not implemented'

    def read_liquid_classes(self):
        # self.system.ReadLiquidClasses()
        return 'Not implemented'

    def reset_stored_adh_info(self):
        # self.system.ResetStoredADHInfo()
        return 'Not implemented'

    def resume(self):
        self.system.Resume()
        return 'System Resumed'

    def set_door_locks(self, close_val):
        self.system.SetDoorLocks(close_val)
        return close_val

    def set_lamp(self, lamp_status):
        # LAMP_OFF=0
        # LAMP_GREEN=1
        # LAMP_GREENFLASHING=2
        # LAMP_REDFLASHING=3
        # self.system.SetLamp(SC_LampStatus Status)
        return 'Not implemented'

    def set_process_variable(self, process_id, variable_name, value):
        # self.system.SetProcessVariable(int ProcessID, string VariableName, SC_Variable Scope, object Value)
        return 'Not implemented'

    def set_rack(self, script_id, index, location, barcode):
        # self.system.SetRack(int ScriptID, int Index, int Location, string Barcode)
        return 'Not implemented'

    def set_remote_mode(self, enable_val):
        # self.system.SetRemoteMode(int Enable)
        return 'Not implemented'

    def set_script_variable(self, script_id, variable_name, value):
        # self.system.SetProcessVariable(int ScriptID, string VariableName, object Value)
        return 'Not implemented'

    def shutdown(self):
        self.logoff()
        self.system.Shutdown()
        return 'Instrument Shutdown'

    def start_adh(self):
        # self.system.StartADH()
        return 'Not implemented'

    def start_adh_load_process(self):
        # self.system.StartADHLoadProcess(SC_EmergencyLevel Emergency)
        # 0 = Normal, 1 = Emergency, 2 = High Emergency
        return 'Not implemented'

    def start_adh_unload_process(self):
        # self.system.StartADHUnloadProcess()
        return 'Not implemented'

    def start_script(self, scriptName, start_line=0, end_line=0):
        # Added the prepare script to the function for convience but the DLL function
        # this is based off of only starts the script, it does not load the script
        self.scriptID = self.system.PrepareScript(scriptName)
        self.system.StartScript(self.scriptID, start_line, end_line)
        return (self.scriptID)

    def stop(self):
        self.system.Stop()
        return 'System Stopped'

    def write_to_trace_file(self):
        # self.system.WriteToTraceFile(string bstrTrace)
        return 'Not implemented'

    def error_event(self):
        returned_value = self.system.ErrorEvent()
        # This requires an ErrorEventHandler
        return returned_value


# Event handling from evoware
class EVOSystemEvents:
    # Currently theses are messages that are being sent after they happened or have
    # been responded to, but as currently being used they do read the message
    # when it arrives for the user to make
    def OnStatusChanged(self, status):
        def debug(message):
            print(message)

        statusString = ""
        for key, value in EVOConstants.STATUS.items():
            if status & key == key:
                statusString += value

    def OnErrorEvent(*args):
        # Need to look into the McsSvr DLLs where these error events are handled
        # just the system DLL does not give access to the correct functions...
        print(args)


# Constants specific to EVOWARE
class EVOConstants:
    normal = 0
    emergency = 1
    high_emergency = 2

    idle_process = 0
    busy_process = 1
    finsihed_process = 2
    error_process = 3
    stopped_process = 4

    unknown_script = 0
    idle_script = 1
    busy_script = 2
    aborted_script = 3
    stopped_script = 4
    pipetting_script = 5
    paused_script = 6
    error_script = 7
    simulation_script = 8
    status_error_script = 9

    unknown_status = 0
    no_instruments_status = 1
    loading_status = 2
    not_initialized_status = 4
    initializing_status = 8
    initialized_status = 16
    shutting_down_status = 32
    shutdown_status = 64
    unloading_status = 128
    running_status = 256
    pause_request_status = 512
    paused_status = 1024
    resource_missing_status = 2048
    deadlock_status = 4096
    execution_error_status = 8192
    time_violation_status = 16384
    idle_status = 65536
    busy_status = 131072
    aborted_status = 262144
    stopped_status = 524288
    pipetting_status = 1048576
    error_status = 2097152
    simulation_status = 4194304
    logon_error_status = 8388608
    connection_error_status = 16777216

    PROCESS_STATUS = {0: "IDLE_PROCESS",
                      1: "BUSY_PROCESS",
                      2: "FINISHED_PROCESS",
                      3: "ERROR_PROCESS",
                      4: "STOPPED_PROCESS"}

    STATUS = {262144: "ABORTED",
              131072: "BUSY",
              16777216: "CONNECTION_ERROR",
              4096: "DEADLOCK",
              2097152: "ERROR",
              8192: "EXECUTION_ERROR",
              65536: "IDLE",
              16: "INITIALIZED",
              8: "INITIALIZING",
              2: "LOADING",
              8388608: "LOGON_ERROR",
              1: "NOINTRUMENTS",
              4: "NOTINITIALIZED",
              1024: "PAUSED",
              512: "PAUSE_REQUESTED",
              1048576: "PIPETTING",
              2048: "RESOURCE_MISSING",
              256: "RUNNING",
              64: "SHUTDOWN",
              32: "SHUTTING_DOWN",
              4194304: "SIMULATION",
              524288: "STOPPED",
              16384: "TIMEVIOLATION",
              0: "UNKNOWN",
              128: "UNLOADING"}


# This class monitors EvoWare for the creation of new gui windows, this does
# allow some windows like timers and user prompts to not trigger this code
class evo_error_watcher:
    def __init__(self):
        self._handle = None
        self._evoware_handle = None
        self._children_windows = {}
        self.evoware_children = []

    def _evoware_enum_callback(self, hwnd, wildcard):
        if re.match(wildcard, str(win32gui.GetWindowText(hwnd))) != None:
            self._evoware_handle = str(hwnd)

    def find_evoware_window(self, wildcard='Freedom EVOware'):
        self._handle = None
        win32gui.EnumWindows(self._evoware_enum_callback, wildcard)

    def _window_enum_callback(self, hwnd, wildcard):
        if str(win32gui.GetParent(hwnd)) not in self._children_windows.keys():
            self._children_windows[str(win32gui.GetParent(hwnd))] = []
        self._children_windows[str(win32gui.GetParent(hwnd))].append([str(win32gui.GetWindowText(hwnd)),
                                                                      hwnd, str(win32gui.GetClassName(hwnd)),
                                                                      str(win32gui.GetParent(hwnd))])

    def find_evoware_child_windows(self, wildcard):
        self._handle = None
        self.evoware_children = []
        win32gui.EnumWindows(self._window_enum_callback, wildcard)
        for window in self._children_windows[str(self._evoware_handle)]:
            if window not in self.evoware_children:
                self.evoware_children.append(window)
        for loop in range(0, 5):
            for evoware_child in self.evoware_children:
                try:
                    for window in self._children_windows[str(evoware_child[1])]:
                        if window not in self.evoware_children:
                            self.evoware_children.append(window)
                except Exception as e:
                    continue

    def check_for_errors(self):
        evo_acceptable_codes = ['Waiting', 'Default IME', 'MSCTFIME UI', 'EVOware remote mode', 'Runtime Controller']
        self.find_evoware_window()
        if self._evoware_handle is None:
            return ['Error', 'Evoware handle not found']
        self.find_evoware_child_windows('*')
        for window in self.evoware_children:
            if window[0] == '':
                continue
            if any(code in window[0] for code in evo_acceptable_codes):
                continue
            elif 'Initializing' in window[0]:
                continue
            else:
                return ['Error', 'Evoware encountered an error: "%s"' % window[0]]
        return ['Success', 'No errors encountered']
