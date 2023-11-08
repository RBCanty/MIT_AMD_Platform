"""
Extra functions to interface between Python and OPUS (Bruker)
Version 0.0.1
Updated on 12/7/2019
@author: Brent
Based off of functions from BrukerPy Library and Bruker Reference Files
"""
import datetime
import os
import time
import re
from dde_client import *

def directory_check(path_to_check):
    if not os.path.exists(path_to_check):
        os.makedirs(path_to_check)
        print('Made a new directory: ' + path_to_check)
        return 1
    else:
        return 0
        
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

def test_opus_connection(taskname, test_mode):
    # Start by checking the connection to OPUS
    # It is quick to do and should return the OPUS Version Date and Instrument
    if test_mode == "0":
        link = None
        link = DDEClient("OPUS","System")
        output1 = link.request("GET_VERSION")
        output2 = link.request("GET_BENCH")
        output = "OPUS Ver. " + str(output1.decode().strip("\n")) + \
                    " Inst. " + str(output2.decode().strip("\n"))
        DDEClient.__del__(link)
        link = None
        return output
    else:
        # print('OPUS not reached, test mode: ' + str(test_mode))
        return -1

def eject_tray_function(inside_home, outside_home, test_mode):
    if test_mode == "0":
        link = None
        link = DDEClient("OPUS","System")
        # First move to inside home position to prepare to eject the tray
        inside = link.request("COMMAND_LINE Stage Control(0,{NOX=" + str(inside_home[0]) + ",NOY=" + str(inside_home[1]) + "});")
        time.sleep(1)
        # Now move to outside position to wait for a new sample plate
        outside = link.request("COMMAND_LINE Stage Control(0,{NOX=" + str(outside_home[0]) + ",NOY=" + str(outside_home[1]) + "});")
        time.sleep(2)
        DDEClient.__del__(link)
        return 1
    else:
        return -1

def inject_tray_function(inside_home, outside_home, test_mode):
    if test_mode == "0":
        link = None
        link = DDEClient("OPUS","System")
        # Return the tray to the inside home position
        inside = link.request("COMMAND_LINE Stage Control(0,{NOX=" + str(inside_home[0]) + ",NOY=" + str(inside_home[1]) + "});")
        time.sleep(2)
        DDEClient.__del__(link)
        return 1
    else:
        return -1

def move_tray_function(position, test_mode):
    if test_mode == "0":
        link = None
        link = DDEClient("OPUS","System")
        movement = link.request("COMMAND_LINE Stage Control(0,{NOX=" + str(position[0]) + ",NOY=" + str(position[1]) + "});")
        time.sleep(2)
        DDEClient.__del__(link)
        return 1
    else:
        return -1

def center_tray_function(inside_home, outside_home, test_mode):
    if test_mode == "0":
        link = None
        link = DDEClient("OPUS","System")
        # First move to inside home position to prepare to eject the tray
        inside = link.request("COMMAND_LINE Stage Control(0,{NOX=" + str(inside_home[0]) + ",NOY=" + str(inside_home[1]) + "});")
        time.sleep(1)
        DDEClient.__del__(link)
        return 1
    else:
        return -1


























