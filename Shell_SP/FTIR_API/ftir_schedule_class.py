# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 09:51:02 2021

@author: kfj_AMD_FTIR
"""

import yaml
import time
import datetime
import queue
import traceback
import os
import re
from dde_client import DDEClient
from dde_client import *
from custom_classes import Narg


def directory_check(path_to_check):
    try:
        if not os.path.exists(path_to_check):
            os.makedirs(path_to_check)
            print('Made a new directory: ' + path_to_check)
            return ['Success', 'Made a new directory!']
        else:
            return ['Success', 'Directory already existed!']
    except:
        return ['Error', traceback.format_exc()]


def format_file_path(datapath):
    now = datetime.datetime.now()
    today = now.strftime("%Y%m%d")
    filepath_out = datapath+today+'\\'
    return filepath_out


def format_file_name(xpm_file, file_path, scan_id):
    now = datetime.datetime.now()
    today = now.strftime("%Y%m%d%H%M%S")
    xpm_prefix = xpm_file[0:3]
    suffix = 0
    while suffix<999:
        filename_out = scan_id + "_" + xpm_prefix + "_" + today + '.' + str(suffix)
        if os.path.exists(file_path + "\\" + filename_out) == False:
            break
        suffix=suffix+1
    return re.sub(r'[^\w\.\-\_]+', '', filename_out)


class ftir_controller:
    def __init__(self, test_mode=Narg()):
        config_filepath = r'C:\Users\kfj_AMD_FTIR\PycharmProjects\AMD_Control_Platform\venv\Shell_SP\FTIR_API\config.yaml'
        try:
            with open(config_filepath,'r') as yaml_file:
                config = yaml.load(yaml_file, Loader=yaml.FullLoader)
        except:
            raise Exception('Configuration file was not found, check directory...')
        self.default_paths = config['paths/files']
        self.startup_values = config['startup']
        if not isinstance(test_mode, Narg):
            self.startup_values['code_test_mode'] = test_mode
        self.hts_values = config['HTS_Default']
        self.my_status = queue.Queue()
        self.my_status.put({'status':'idle', 'tray':'closed'})
        try:
            test_opus_return = self.test_opus_connection()
        except Exception as e:
            raise Exception('Connection to OPUS not made: %s' % e)
        if test_opus_return[0] != 'Success':
            raise Exception('Problem with OPUS: %s' % test_opus_return[1])
    
    def current_log_folder(self):
        current_date = datetime.datetime.now()
        log_path = self.default_paths['python_log_path']
        log_path_date = log_path + current_date.strftime("%Y%m%d")
        log_path_return = directory_check(log_path_date)
        if log_path_return[0] != 'Success':
            return ['Error', log_path_return[1]]
        self.current_log_path = log_path_date
        return ['Success', self.current_log_path]
           
    def test_opus_connection(self, taskname = 'test_dde', test_mode_input = 'default'):
        # Start by checking the connection to OPUS
        # It is quick to do and should return the OPUS Version Date and Instrument
        if test_mode_input == 'default':
            test_mode = self.startup_values['code_test_mode']
        else:
            test_mode = test_mode_input
        if test_mode == "0":
            link = None
            link = DDEClient("OPUS","System")
            output1 = link.request("GET_VERSION")
            output2 = link.request("GET_BENCH")
            output = "OPUS Ver. " + str(output1.decode().strip("\n")) + \
                        " Inst. " + str(output2.decode().strip("\n"))
            DDEClient.__del__(link)
            link = None
            if str(output1.decode().strip("\n")) == '20190310':
                return ['Success', 'OPUS reached']
            else:
                return ['Error', 'Incorrect OPUS Version reached: %s' % str(output1.decode().strip("\n"))]
            return output
        else:
            return ['Testing', 'OPUS not reached, test mode: %s' % str(test_mode)]
    
    def eject_tray_function(self, test_mode_input = 'default'):
        if test_mode_input == 'default':
            test_mode = self.startup_values['code_test_mode']
        else:
            test_mode = test_mode_input
        inside_home = self.hts_values['inside_home_location']
        outside_home =self.hts_values['outside_location']
        if test_mode == "0":
            test_opus_return = self.test_opus_connection()
            if test_opus_return[0] != 'Success':
                return ['Error', 'Lost communication with OPUS: %s' % test_opus_return[1]]
            link = None
            link = DDEClient("OPUS","System")
            # First move to inside home position to prepare to eject the tray
            inside = link.request("COMMAND_LINE Stage Control(0,{NOX=" + str(inside_home[0]) + ",NOY=" + str(inside_home[1]) + "});")
            if 'OK' not in inside.decode():
                return ['Error', 'Eject command not properly executed: %s' % str(inside)]
            time.sleep(1)
            # Now move to outside position to wait for a new sample plate
            outside = link.request("COMMAND_LINE Stage Control(0,{NOX=" + str(outside_home[0]) + ",NOY=" + str(outside_home[1]) + "});")
            if 'OK' not in outside.decode():
                return ['Error', 'Eject command not properly executed: %s' % str(outside)]
            time.sleep(2)
            DDEClient.__del__(link)
            return ['Success', 'Tray was ejected correctly based on OPUS returns']
        else:
            return ['Testing', 'OPUS not reached, test mode: %s' % str(test_mode)]
    
    def inject_tray_function(self, test_mode_input = 'default'):
        if test_mode_input == 'default':
            test_mode = self.startup_values['code_test_mode']
        else:
            test_mode = test_mode_input
        inside_home = self.hts_values['inside_home_location']
        if test_mode == "0":
            test_opus_return = self.test_opus_connection()
            if test_opus_return[0] != 'Success':
                return ['Error', 'Lost communication with OPUS: %s' % test_opus_return[1]]
            link = None
            link = DDEClient("OPUS","System")
            # Return the tray to the inside home position
            inside = link.request("COMMAND_LINE Stage Control(0,{NOX=" + str(inside_home[0]) + ",NOY=" + str(inside_home[1]) + "});")
            if 'OK' not in inside.decode():
                return ['Error', 'Inject command not properly executed: %s' % str(inside)]
            time.sleep(2)
            DDEClient.__del__(link)
            return ['Success', 'Tray was injected correctly based on OPUS returns']
        else:
            return ['Testing', 'OPUS not reached, test mode: %s' % str(test_mode)]
    
    def move_tray_function(self, position, test_mode_input = 'default'):
        if test_mode_input == 'default':
            test_mode = self.startup_values['code_test_mode']
        else:
            test_mode = test_mode_input
        if test_mode == "0":
            test_opus_return = self.test_opus_connection()
            if test_opus_return[0] != 'Success':
                return ['Error', 'Lost communication with OPUS: %s' % test_opus_return[1]]
            link = None
            link = DDEClient("OPUS","System")
            movement = link.request("COMMAND_LINE Stage Control(0,{NOX=" + str(position[0]) + ",NOY=" + str(position[1]) + "});")
            if 'OK' not in movement.decode():
                return ['Error', 'Inject command not properly executed: %s' % str(movement)]
            time.sleep(2)
            DDEClient.__del__(link)
            return ['Success', 'Tray was moved correctly based on OPUS returns']
        else:
            return ['Testing', 'OPUS not reached, test mode: %s' % str(test_mode)]
    
    def run_opus_background_measure(self, xpmname, macro_name):
        test_mode = self.startup_values['code_test_mode']
        if test_mode == "0":
            end_of_file_character = chr(13) + chr(10)
            macroname = self.default_paths['macrofilespath'] + macro_name
            xpmpath = self.default_paths['xpmpath']
            param_num = " 2"
            macro_request = "RUN_MACRO " + macroname + param_num
            with open(self.current_log_path + r'\execution.log', 'a') as execution_log_file:
                execution_log_file.write('\t'+macro_request+'\n')
                execution_log_file.write('\t\t'+'Experiment: '+xpmpath+xpmname+'\n')
                execution_log_file.flush()

            test_opus_return = self.test_opus_connection()
            if test_opus_return[0] != 'Success':
                return ['Error', 'Lost communication with OPUS: %s' % test_opus_return[1]]
            link = None
            id1 = None
            link = DDEClient("OPUS","System")
            request_mode = link.request("REQUEST_MODE")
            link.request(macro_request)
            link.request(xpmname + end_of_file_character)
            id1 = link.request(xpmpath + end_of_file_character)
            macro_id = str(id1[4:-1].decode())
            with open(self.current_log_path + r'\execution.log', 'a') as execution_log_file:
                execution_log_file.write('\t\tRequest sent with macro ID: ' + macro_id + '\n')
            time.sleep(5)
            try:
                while link.request("MACRO_RESULTS "+macro_id).decode().split('\n')[2] == '0':
                    time.sleep(1)
                with open(self.current_log_path + r'\execution.log', 'a') as execution_log_file:
                    execution_log_file.write('\t\tCompleted execution of Macro: ' + macroname + '\n')
            except:
                with open(self.current_log_path + r'\execution.log', 'a') as execution_log_file:
                    execution_log_file.write('\t\tOPUS polling error, I will wait a little more and continue??' + '\n')
                time.sleep(5)
            DDEClient.__del__(link)
            link = None
            id1 = None
            return ['Success', 'Based on OPUS returns, background measurement completed']
        else:
            return ['Testing', 'OPUS not reached, test mode: %s' % str(test_mode)]

    def run_opus_sample_measurement(self, xpmname, outputname, outputpath, macro_name):
        test_mode = self.startup_values['code_test_mode']
        if test_mode == "0":
            # Here we need to define things that the OPUS Macros are expecting
            end_of_file_character = chr(13) + chr(10)
            macroname = self.default_paths['macrofilespath'] + macro_name
            xpmpath = self.default_paths['xpmpath']
            param_num = " 4"
            macro_request = "RUN_MACRO " + macroname + param_num
            with open(self.current_log_path + r'\execution.log', 'a') as execution_log_file:
                execution_log_file.write('\t'+macro_request+'\n')
                execution_log_file.write('\t\t'+'Experiment: '+xpmpath+xpmname+'\n')
                execution_log_file.write('\t\t'+'Output: '+outputpath+outputname+'\n')
                execution_log_file.flush()
            
            test_opus_return = self.test_opus_connection()
            if test_opus_return[0] != 'Success':
                return ['Error', 'Lost communication with OPUS: %s' % test_opus_return[1]]
            link = None
            id1 = None
            link = DDEClient("OPUS","System")
            request_mode = link.request("REQUEST_MODE")
            link.request(macro_request)
            link.request(outputpath + end_of_file_character)
            link.request(outputname + end_of_file_character)
            link.request(xpmname + end_of_file_character)
            id1 = link.request(xpmpath + end_of_file_character)
            macro_id = str(id1[4:-1].decode())
            with open(self.current_log_path + r'\execution.log', 'a') as execution_log_file:
                execution_log_file.write('\t\tRequest sent with macro ID: ' + macro_id + '\n')
            time.sleep(5)
            try:
                while link.request("MACRO_RESULTS "+macro_id).decode().split('\n')[2] == '0':
                    time.sleep(1)
                with open(self.current_log_path + r'\execution.log', 'a') as execution_log_file:
                    execution_log_file.write('\t\tCompleted execution of Macro: ' + macroname + '\n')
            except:
                with open(self.current_log_path + r'\execution.log', 'a') as execution_log_file:
                    execution_log_file.write('\t\tOPUS polling error, I will wait a little more and continue??' + '\n')
                time.sleep(5)
            DDEClient.__del__(link)
            link = None
            id1 = None
            return ['Success', 'Based on OPUS returns, sample measurement completed']
        else:
            return ['Testing', 'OPUS not reached, test mode: %s' % str(test_mode)]
    
    def convert_to_xy_data(self, outputfile, outputpath, xyfile):
        test_mode = self.startup_values['code_test_mode']
        if test_mode == "0":
            # Here we need to define things that the OPUS Macros are expecting
            end_of_file_character = chr(13) + chr(10)
            macroname = self.default_paths['macrofilespath'] + 'CONTOXY.mtx'
            param_num = " 3"
            macro_request = "RUN_MACRO " + macroname + param_num
            with open(self.current_log_path + r'\execution.log', 'a') as execution_log_file:
                execution_log_file.write('\t\t'+macro_request+'\n')
                execution_log_file.flush()
                
            test_opus_return = self.test_opus_connection()
            if test_opus_return[0] != 'Success':
                return ['Error', 'Lost communication with OPUS: %s' % test_opus_return[1]]
            link = None
            id1 = None
            link = DDEClient("OPUS","System")
            request_mode = link.request("REQUEST_MODE")
            link.request(macro_request)
            link.request(outputpath + end_of_file_character)
            link.request(outputfile + end_of_file_character)
            id1 = link.request(xyfile + end_of_file_character)
            macro_id = str(id1[4:-1].decode())
            with open(self.current_log_path + r'\execution.log', 'a') as execution_log_file:
                execution_log_file.write('\t\tRequest sent with macro ID: ' + macro_id + '\n')
            time.sleep(1)
            try:
                while link.request("MACRO_RESULTS "+macro_id).decode().split('\n')[2]=='0':
                    time.sleep(1)
                with open(self.current_log_path + r'\execution.log', 'a') as execution_log_file:
                    execution_log_file.write('\t\tCompleted execution of Macro: ' + macroname + '\n')
            except:
                with open(self.current_log_path + r'\execution.log', 'a') as execution_log_file:
                    execution_log_file.write('\t\tOPUS polling error, I will wait a little more and continue??' + '\n')
                time.sleep(5)
            DDEClient.__del__(link)
            link = None
            id1 = None
            return ['Success', 'Based on OPUS returns, convert to xy completed']
        else:
            return ['Testing', 'OPUS not reached, test mode: %s' % str(test_mode)]
    
    def prepare_send_recieve(self):
        print(list(self.my_status.queue))
        ftir_status = self.my_status.queue[0]
        if ftir_status['status'] == 'busy':
            return ['Error', 'FTIR busy']
        elif ftir_status['status'] == 'error':
            return ['Error', 'FTIR error']
        if ftir_status['tray'] == 'closed':
            return_statements = self.eject_tray_function()
            if return_statements[0] != 'Success':
                return ['Error', return_statements[1]]
            ftir_status['tray'] = 'open'
        elif ftir_status['tray'] == 'open':
            pass
        else:
            return ['Error', 'Unknown state present in tray-status']
        self.my_status.get()
        self.my_status.put(ftir_status)
        print(list(self.my_status.queue))
        
    def return_to_initial_state(self):
        print(list(self.my_status.queue))
        ftir_status = self.my_status.queue[0]
        if ftir_status['status'] == 'busy':
            return ['Error', 'FTIR busy']
        elif ftir_status['status'] == 'error':
            return ['Error', 'FTIR error']
        if ftir_status['tray'] == 'open':
            return_statements = self.inject_tray_function()
            if return_statements[0] != 'Success':
                return ['Error', return_statements[1]]
            ftir_status['tray'] = 'closed'
        elif ftir_status['tray'] == 'closed':
            pass
        else:
            return ['Error', 'Unknown state present in tray-status']
        self.my_status.get()
        self.my_status.put(ftir_status)
        print(list(self.my_status.queue))
    
    def execute_ir_experiment(self, incoming_dictionary, xpm = 'HTS_test_method_v2.XPM'):
        
        print(list(self.my_status.queue))
        ftir_status = self.my_status.queue[0]
        if ftir_status['status'] == 'busy':
            return ['Error', 'FTIR busy']
        elif ftir_status['status'] == 'error':
            return ['Error', 'FTIR error']
        else:
            ftir_status['status'] = 'busy'
        if ftir_status['tray'] == 'open':
            return_statements = self.inject_tray_function()
            if return_statements[0] != 'Success':
                return ['Error', return_statements[1]]
            ftir_status['tray'] = 'closed'
        elif ftir_status['tray'] == 'closed':
            pass
        else:
            return ['Error', 'Unknown state present in tray-status']
        self.my_status.get()
        self.my_status.put(ftir_status)
        print(list(self.my_status.queue))
        
        #if the tray is already out inject the tray to get ready for the measurement
        
        # Initialize the current log paths for execution
        self.current_log_folder()
        start_time = datetime.datetime.now()
        print(start_time)
        # Now open the file that contains the positions of the wells
        well_position_path = self.default_paths['wellpositionpath']
        try:
            with open(well_position_path, 'r') as yaml_file:
                well_positions = yaml.load(yaml_file, Loader=yaml.FullLoader)
        except FileNotFoundError:
            return ['Error', 'Well-position file was not found, check directory...']
        
        # Now get the output folder ready for OPUS to dump data into
        # OPUS does not make it's own folder directories
        output_path = format_file_path(self.default_paths['datapath'])
        output_return = directory_check(output_path)
        if output_return[0] != 'Success':
            return ['Error', output_return[1]]
        
        with open(self.current_log_path + r'\execution.log', 'a') as execution_log_file:
            execution_log_file.write('=========================================================\n')
            execution_log_file.write('Started schedule at ' + start_time.strftime("%Y%m%d-%H:%M:%S") + '\n')
        # Start executing the schedule from the incoming dictionary
        data_files_generated = []
        
        # The first step is to check the connection to OPUS
        test_opus_return = self.test_opus_connection()
        if test_opus_return[0] != 'Success':
            return ['Error', 'Lost communication with OPUS: %s' % test_opus_return[1]]
        
        # Then we check the incoming dictionary for blanks to run
        if 'blanks' in incoming_dictionary.keys():
            blank_wells = incoming_dictionary['blanks']
            for well in blank_wells.keys():
                move_return = self.move_tray_function(well_positions[well])
                if move_return[0] != 'Success':
                    return ['Error', move_return[1]]
                # Now we execute the actual experiment
                back_xpm = xpm
                back_macro = 'RUNBACK.MTX'
                backmeasure_results = self.run_opus_background_measure(back_xpm, back_macro)
                if backmeasure_results[0] != 'Success':
                    return ['Error', backmeasure_results[1]]
                print('\t'+'Background measurement executed: well ' + well)
        
        if 'samples' in incoming_dictionary.keys():
            sample_wells = incoming_dictionary['samples']
            for well in sample_wells.keys():
                move_return = self.move_tray_function(well_positions[well])
                if move_return[0] != 'Success':
                    return ['Error', move_return[1]]
                
                output_file_path = output_path
                output_file_name = format_file_name(sample_wells[well]['sample_name'], output_file_path,
                                                    '%s_%s_%s' %(well, sample_wells[well]['support'], 
                                                                 sample_wells[well]['matrix']))
                sample_xpm = xpm
                sample_macro = 'RUNMEAS.MTX'
                samplemeasure_results = self.run_opus_sample_measurement(sample_xpm, output_file_name, output_file_path, sample_macro)
                if samplemeasure_results[0] != 'Success':
                    return ['Error', samplemeasure_results[1]]
                print('\t'+'Sample measurement executed: %s %s' % (well, sample_wells[well]['sample_name']))
                output_xy_name = output_file_name.replace('.','_') + '.txt'
                data_files_generated.append(output_xy_name)
                convert_to_xy = self.convert_to_xy_data(output_file_name, output_file_path, output_xy_name)
                if convert_to_xy[0] != 'Success':
                    return ['Error', convert_to_xy[1]]
                with open(self.current_log_path + r'\execution.log', 'a') as execution_log_file:
                    execution_log_file.write('\t\t'+'Datafile Output: ' + output_file_path + output_xy_name + '\n')
                print('\t'+'Sample measurement converted to X-Y data: %s %s' % (well, sample_wells[well]['sample_name']))
                
        end_time = datetime.datetime.now()
        with open(self.current_log_path + r'\execution.log', 'a') as execution_log_file:
            execution_log_file.write('Completed the schedule at ' + end_time.strftime("%Y%m%d-%H:%M:%S") + '\n')
            execution_log_file.write('=========================================================\n')
        print('Completed the schedule at %s' % end_time.strftime("%Y%m%d-%H:%M:%S"))
        
        print(list(self.my_status.queue))
        ftir_status = self.my_status.queue[0]
        if ftir_status['status'] == 'error':
            return ['Error', 'FTIR error']
        elif ftir_status['status'] == 'idle':
            pass
        elif ftir_status['status'] == 'busy':
            ftir_status['status'] = 'idle'
        else:
            return ['Error', 'Unknown status present']
        self.my_status.get()
        self.my_status.put(ftir_status)
        print(list(self.my_status.queue))





testing_input = {'blanks' : {'B2': {'sample_name':'BLANK'}},
                 'samples' : {'B4': {'sample_name':'MaxSample1_5uL', 'matrix':'air', 'support':'Si'}, 
                              'B6': {'sample_name':'MaxSample1_2uL', 'matrix':'air', 'support':'Si'}}}

test = ftir_controller()
print(list(test.my_status.queue))
print(test.test_opus_connection())

print(test.prepare_send_recieve())
input("press enter")
print(test.current_log_folder())
print(test.execute_ir_experiment(testing_input))
print(test.prepare_send_recieve())
input("press enter")
print(test.return_to_initial_state())





