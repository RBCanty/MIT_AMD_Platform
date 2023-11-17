"""
Main Level of Code to interface between Python and OPUS
Version 0.0.1
Updated on 12/7/2019
@author: Brent
Based on ASAP Implementation and Bruker Reference Files
"""

import os
import datetime
import sys
import yaml
import pandas as pd
import time
from dde_client import *
from extra_function_library import *
import json
from pprint import pprint


# TODO Actually write the thing.... you can do it! :)
def prepare_send_and_recieve():
    # if the tray is open and empty everything is great, otherwise spit the tray out or return an error
    pass


def amd_ftir(schedule_object):
    def run_opus_sample_measurement(xpmname, xpmpath, outputname, outputpath, macro_name, macro_path):
        if startup_values['code_test_mode'] == "0":
            endoffile=chr(13)+chr(10)
            macroname = macro_path + macro_name
            param_num = " 4"
            macro_request = "RUN_MACRO " + macroname + param_num
            execution_log_file.write('\t\t'+macro_request+'\n')
            execution_log_file.write('\t\t\t'+'Experiment: '+xpmpath+xpmname+'\n')
            execution_log_file.write('\t\t\t'+'Output: '+outputpath+outputname+'\n')
            execution_log_file.flush()
            link = None
            id1 = None
            link = DDEClient("OPUS","System")
            request_mode = link.request("REQUEST_MODE")
            #print(macro_request)
            link.request(macro_request)
            link.request(outputpath+endoffile)
            link.request(outputname+endoffile)
            link.request(xpmname+endoffile)
            id1 = link.request(xpmpath+endoffile)
            macro_id = str(id1[4:-1].decode())
            #print('Request sent with macro ID: ' + macro_id)
            time.sleep(5)
            try:
                while link.request("MACRO_RESULTS "+macro_id).decode().split('\n')[2]=='0':
                    time.sleep(1)
                print('\t'+'Completed execution of Macro: ' + macroname)
            except:
                print("OPUS polling error, I will wait a little more and continue??")
                time.sleep(5)
            DDEClient.__del__(link)
            link = None
            id1 = None
            return 1
        else:
            print('OPUS not reached, test mode: ' + str(startup_values['code_test_mode']))
            return -1
    
    def run_opus_background_measure(xpmname, xpmpath, macro_name, macro_path):
        if startup_values['code_test_mode'] == "0":
            endoffile=chr(13)+chr(10)
            macroname = macro_path + macro_name
            param_num = " 2"
            macro_request = "RUN_MACRO " + macroname + param_num
            execution_log_file.write('\t\t'+macro_request+'\n')
            execution_log_file.write('\t\t\t'+'Experiment: '+xpmpath+xpmname+'\n')
            execution_log_file.flush()
            #print(macro_request)
            link = None
            id1 = None
            link = DDEClient("OPUS","System")
            request_mode = link.request("REQUEST_MODE")
            link.request(macro_request)
            link.request(xpmname+endoffile)
            id1 = link.request(xpmpath+endoffile)
            macro_id = str(id1[4:-1].decode())
            #print('Request sent with macro ID: ' + macro_id)
            time.sleep(5)
            try:
                while link.request("MACRO_RESULTS "+macro_id).decode().split('\n')[2]=='0':
                    time.sleep(1)
                print('\t'+'Completed execution of Macro: ' + macroname)
            except:
                print("OPUS polling error, I will wait a little more and continue??")
                time.sleep(5)
            DDEClient.__del__(link)
            link = None
            id1 = None
            return 1
        else:
            print('OPUS not reached, test mode: ' + str(startup_values['code_test_mode']))
            return -1
    
    def convert_to_xy_data(outputfile, outputpath, xyfile):
        if startup_values['code_test_mode'] == "0":
            endoffile=chr(13)+chr(10)
            macroname = default_paths['macrofilespath'] + 'CONTOXY.mtx'
            param_num = " 3"
            macro_request = "RUN_MACRO " + macroname + param_num
            execution_log_file.write('\t\t\t'+macro_request+'\n')
            execution_log_file.flush()
            #print(macro_request)
            link = None
            id1 = None
            link = DDEClient("OPUS","System")
            request_mode = link.request("REQUEST_MODE")
            link.request(macro_request)
            link.request(outputpath+endoffile)
            link.request(outputfile+endoffile)
            id1 = link.request(xyfile+endoffile)
            macro_id = str(id1[4:-1].decode())
            #print('Request sent with macro ID: ' + macro_id)
            time.sleep(1)
            try:
                while link.request("MACRO_RESULTS "+macro_id).decode().split('\n')[2]=='0':
                    time.sleep(1)
                print('\t'+'Completed execution of Macro: ' + macroname)
            except:
                print("OPUS polling error, I will wait a little more and continue??")
                time.sleep(5)
            DDEClient.__del__(link)
            link = None
            id1 = None
        return 1
    
    def schedule_builder(string_schedule):
        schedule_name = string_schedule['name']
        experiment = 'HTS_test_method.XPM'
        back_macro = 'RUNBACK.MTX'
        sample_macro = 'RUNMEAS.MTX'
        schedule = []
        schedule.append('\t'.join(['test_dde','None','None','None','None']))
        schedule.append('\t'.join(['eject_tray','None','None','None','None']))
        schedule.append('\t'.join(['inject_tray','None','None','None','None']))
        if len(string_schedule['background']) != 0:
            for element in string_schedule['background']:
                if element[1] == 'Same Plate':
                    schedule.append('\t'.join(['move_tray',element[0],'None','None','None']))
                    schedule.append('\t'.join(['back_measure',element[0],back_macro,experiment,'None']))
                elif element[1] == 'Diff Plate':
                    schedule.append('\t'.join(['eject_tray','None','None','None','None']))
                    schedule.append('\t'.join(['inject_tray','None','None','None','None']))
                    schedule.append('\t'.join(['move_tray',element[0],'None','None','None']))
                    schedule.append('\t'.join(['back_measure',element[0],back_macro,experiment,'None']))
        if len(string_schedule['sample']) != 0:
            for element in string_schedule['sample']:
                schedule.append('\t'.join(['move_tray',element[0],'None','None','None']))
                schedule.append('\t'.join(['sample_measure',element[0],sample_macro,experiment,element[1]]))
        schedule.append('\t'.join(['eject_tray','None','None','None','None']))
        schedule.append('\t'.join(['inject_tray','None','None','None','None']))
        built_schedule_name = 'FTIR_Schedule_' + schedule_name + '.txt'
        schedule_path = default_paths['schedulespath'] + 'FTIR_Schedule_' + schedule_name + '.txt'
        with open(schedule_path, 'w') as schedule_file:
            for element in schedule:
                schedule_file.write(element + '\n')
        return built_schedule_name
    
    """
    Right now the schedule is accepted as a sys.stdin.read() argument from the terminal
    Example of sequence to start the execution: echo $variable | python '...\ftir_schedule_execute.py'
    Or if not variable: echo 'STRING' | python '...\ftir_schedule_execute.py'
    """
    
    # Start by reading in the configuration file from the directory
    # Gives us some default paths and variables to start with
    os.chdir('E:\\AMD_FTIR_v2\\')
    
    #schedule_to_run = 'testschedule.txt'
    config_filename = os.getcwd() + '\\config.yaml'
    try:
        with open(config_filename,'r') as yaml_file:
            config = yaml.load(yaml_file, Loader=yaml.FullLoader)
    except FileNotFoundError:
        print('Configuration file was not found, check directory...')
        sys.exit(-1)
    default_paths = config['paths/files']
    startup_values = config['startup']
    hts_values = config['HTS_Default']
    
    if len(sys.argv) != 1:
        print('Give me the schedule json dictionary to execute')
    else:
        incoming = schedule_object
        schedule_to_run = schedule_builder(incoming)
    
    # Prepare the python logs for the measurement (Execution and Error logs)
    # Will check for a directory and make one if needed
    current_date = datetime.datetime.now()
    log_path = default_paths['python_log_path']
    log_path_date = log_path + current_date.strftime("%Y%m%d")
    log_path_return = directory_check(log_path_date)
    if log_path_return == 1:
        print('New log folder made!')
    execution_log_path = log_path_date + '\\execution.log'
    error_log_path = log_path_date + '\\error.log'
    #error_log_file = open(error_log_path,'a')
    #error_log_file.write('=========================================================\n')
    #error_log_file.write('Started schedule: ' + schedule_path + ' at ' + current_date.strftime("%Y%m%d-%H:%M:%S") + '\n')
    
    
    # Open the schedule that needs to be executed
    # The code just prints a function and ends if there is not a schedule
    schedule_path = default_paths['schedulespath'] + schedule_to_run
    try:
        schedule = pd.read_csv(schedule_path, delimiter='\t', header=None)
    except FileNotFoundError:
        print(schedule_path)
        print('Schedule file was not found, check directory...')
        #error_log_file.write('\t'+'Schedule file was not found, check directory...')
    
    # Open the well-plate locations dictionary
    # The code prints a problem and ends if there is not a well position configuration
    well_position_path = default_paths['wellpositionpath']
    try:
        with open(well_position_path, 'r') as yaml_file:
            well_positions = yaml.load(yaml_file, Loader=yaml.FullLoader)
    except FileNotFoundError:
        print('Well-position file was not found, check directory...')
        #error_log_file.write('\t'+'Well-position file was not found, check directory...')
        sys.exit(-1)
    
    # Now that all of the paths and configs are loaded, move to executing the schedule
    execution_log_file = open(execution_log_path,'a')
    execution_log_file.write('=========================================================\n')
    execution_log_file.write('Started schedule: ' + schedule_path + ' at ' + current_date.strftime("%Y%m%d-%H:%M:%S") + '\n')
    execution_log_file.flush()
    
    # Prepare the output path for the rest of the script with the pre-sets loaded
    # Right now it saves to def_paths_files['datapath']+now.strftime("%Y%m%d")+'\\'
    # This is outside of the loop for scans that span multiple days
    output_path = format_file_path(default_paths['datapath'])
    output_return = directory_check(output_path)
    if output_return == 1:
        print('New data folder made!')
    
    #sys.exit(-1)
    # Start executing the schedule that was loaded above... hopefully...
    data_files_generated = []
    for index in range(0, len(schedule)):
        task_to_run = schedule[0][index]
        start_time = datetime.datetime.now()
        start_string = 'Started at: ' + start_time.strftime("%H:%M:%S.%f") + '  ' + task_to_run
        print(start_string)
        execution_log_file.write('\t' + start_string + '\n')
        execution_log_file.flush()
        if task_to_run == 'test_dde':
            test_connection_results = test_opus_connection(task_to_run, startup_values['code_test_mode'])
            if type(test_connection_results) == int:
                print('    Could not connect to DDE Client, ending script: Return ' + str(test_connection_results))
                #error_log_file.write('    Could not connect to DDE Client, ending script: Return ' + str(test_connection_results))
                # sys.exit('Error')
            else:
                print('\t'+'Found connection to DDE Client, continuing script')
                continue
                
        elif task_to_run == 'eject_tray':
            eject_tray_results = eject_tray_function(hts_values['inside_home_location'], hts_values['outside_location'], startup_values['code_test_mode'])
            if eject_tray_results != 1:
                print('    Could not move the stage, ending script: Return ' + str(eject_tray_results))
                #error_log_file.write('    Could not move the stage, ending script: Return ' + str(eject_tray_results))
                # sys.exit('Error')
            input('Insert the new sample tray and press enter to continue schedule execution!')
        
        elif task_to_run == 'inject_tray':
            inject_tray_results = inject_tray_function(hts_values['inside_home_location'], hts_values['outside_location'], startup_values['code_test_mode'])
            if inject_tray_results != 1:
                print('    Could not move the stage, ending script: Return ' + str(inject_tray_results))
                #error_log_file.write('    Could not move the stage, ending script: Return ' + str(inject_tray_results))
                # sys.exit('Error')
        
        elif task_to_run == 'move_tray':
            try:
                move_tray_results = move_tray_function(well_positions[schedule[1][index]], startup_values['code_test_mode'])
            except KeyError:
                print('Check the well-positions that are loaded in the dictionary...')
                #error_log_file.write('Check the well-positions that are loaded in the dictionary correctly...')
                # sys.exit('Error')
        
        elif task_to_run == 'back_measure':
            back_xpm = schedule[3][index]
            back_xpm_path = default_paths['xpmpath']
            back_macro = schedule[2][index]
            back_macro_path = default_paths['macrofilespath']
            backmeasure_results = run_opus_background_measure(back_xpm, back_xpm_path, back_macro, back_macro_path)
            if backmeasure_results == -1:
                print('OPUS not reached, test mode: ' + str(startup_values['code_test_mode']))
            else:
                print('\t'+'Background measurement executed: Cell ' + schedule[1][index])
            
        elif task_to_run == 'sample_measure':
            output_file_path = output_path
            output_file_name = format_file_name(schedule[3][index], output_file_path, schedule[4][index]+"_"+schedule[1][index])
            sample_xpm = schedule[3][index]
            sample_xpm_path = default_paths['xpmpath']
            sample_macro = schedule[2][index]
            sample_macro_path = default_paths['macrofilespath']
            samplemeasure_results = run_opus_sample_measurement(sample_xpm, sample_xpm_path, output_file_name, output_file_path, sample_macro, sample_macro_path)
            output_xy_name = output_file_name.replace('.','_') + '.txt'
            data_files_generated.append(output_xy_name)
            convert_to_xy = convert_to_xy_data(output_file_name, output_file_path, output_xy_name)
            execution_log_file.write('\t\t\t'+'Datafile Output: ' + output_file_path + output_xy_name + '\n')
            if samplemeasure_results == -1:
                print('OPUS not reached, test mode: ' + str(startup_values['code_test_mode']))
            else:
                print('\t'+'Sample measurement executed and converted to X-Y data: ' + schedule[1][index] + ' ' + schedule[4][index])
    
        else:    
            print('\tNo methods setup for : ' + schedule[0][index])
            execution_log_file.write('\t' + 'Nothing for : ' + schedule[0][index] + '\n')
            execution_log_file.flush()
            # error_log_file.write('\t' + 'Nothing for : ' + schedule[0][index] + '\n')
            continue
        
        end_time=datetime.datetime.now()
        end_string = 'Ended at: ' + end_time.strftime("%Y/%d/%m %H:%M:%S.%f") + '  ' + task_to_run
        print('\t'+end_string)
        execution_log_file.write('\t\t'+end_string+'\n')
        execution_log_file.flush()
    print('Completed the schedule: ' + schedule_path)
    execution_log_file.write('Completed the schedule: ' + schedule_path + '\n')
    execution_log_file.write('=========================================================\n')
    execution_log_file.flush()
    execution_log_file.close()
    # error_log_file.write('Completed the schedule: ' + schedule_path + '\n')
    # error_log_file.write('=========================================================\n')
    # error_log_file.close()
    return -1

"""
Todo: Check the home position and the outside location positions
      Get the positions of all wells on the silicon plates and make dictionary
      How to handle OPUS when it crashes?
      Add useful file format saving to the OPUS Macro
      Write a sample measurement OPUS Macro
      Write a reference measurement OPUS Macro


Build a loop/structure to handle taking task list from file and executing the movements
1. Test Connection
---Check if reference measurement is needed---
2. Open Tray command - move to eject/inject location
    First shift tray to inject/eject location inside
    Send outward to outside location
3. Wait for signal (temp stub for eventlistener)
4. Close Tray upon signal
---Reference Structure---
5. Move to reference location
6. Measure reference spectrum
---End Reference Structure---
---Measurement loop---
7. Move to measurement location
    link.request("COMMAND_LINE GetStagePosition (0, {NDV=''TangoPCIe_HTSXT''});")
    link.request("COMMAND_LINE Stage Control(0,{NOX=2.75E5,NOY=6.588E4});")
    link.request("COMMAND_LINE Stage Control(0,{NOX=0,NOY=6.588E4});")
8. Take measurement
9. Save the file in a useful format
---Finish Mesure Loop---
10. Open Tray command
11. Wait for signal (temp stub for eventlistener)
12. Close Tray command
13. Event finished signal


    elif task_to_run == 'connect_dde':
        if startup_values['code_test_mode'] == "0":
            link = None
            link = DDEClient("OPUS","System")
    
    elif task_to_run == 'disconnect_dde':
        if startup_values['code_test_mode'] == "0":
            time.sleep(5)
            DDEClient.__del__(link)
"""