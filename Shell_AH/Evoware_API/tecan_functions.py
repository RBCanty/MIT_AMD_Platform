# -*- coding: utf-8 -*-
"""
Compilation of Tecan functions to be used in methods, this is based on
the useful functions that the liquid handler can execute that seemed to
be the most useful, this is not a full listing of all functions

@author: Brent
"""

from pprint import pprint
import traceback
import os

try:
    from useful_other_functions import check_for_problems
except ModuleNotFoundError:
    from Evoware_API.useful_other_functions import check_for_problems


def well_name_to_string(well_list, loc_type, carr_library):
    def find_in_list_of_list(mylist, char):
        for sub_list in mylist:
            for items in sub_list:
                if char in items:
                    return [mylist.index(sub_list), sub_list.index(char)]
        return ['Error', "'{char}' is not in list test".format(char=char)]
    if 'Trough' in loc_type:
        nrows = 8
        ncols = 1
        nwells = nrows * ncols
    else:
        try:
            nrows = carr_library['labware'][loc_type]['nrows']
            ncols = carr_library['labware'][loc_type]['ncols']
            nwells = nrows * ncols
        except Exception:
            return ['Error',  traceback.format_exc()]
    
    wells = []
    targets = []
    col, row = 1, 1
    for i in range(0, nwells):
        wells.append(chr(ord('A') + row - 1) + str(col))
        targets.append(0)
        if row == nrows:
            row = 1
            col += 1
        else:
            row += 1
    first_piece = '{:02x}'.format(ncols).upper() + '{:02x}'.format(nrows).upper()
    n = 7
    well_groups = [wells[i:i + n] for i in range(0, nwells, n)]
    targets_groups = [targets[i:i + n] for i in range(0, nwells, n)]
    string_groups = ['0'] * len(well_groups)
    for single_well in well_list:
        well_index = find_in_list_of_list(well_groups, single_well)
        if well_index[0] == 'Error':
            return ['Error', well_index[1]]
        targets_groups[well_index[0]][well_index[1]] = 1
    for a in range(0, len(targets_groups)):
        group = targets_groups[a]
        dec_value = sum([2 ** i for i in range(0, len(group)) if group[i] > 0])
        ascii_value = chr(dec_value + 48)
        string_groups[a] = ascii_value
    final_string = first_piece + ''.join(string_groups)
    return ['Success', final_string]


def transfer_labware(current_bed, wellplate_library, transfer_command):
    # Values that get sent into this expression
    # 0. Carrier Grid Position Source (1-67) [String]
    # 1. Destination Grid Position (1-67) [String]
    # 2. Move back home? (1 = true, 0 = false) [Integer]
    # 3. Lid Handling? (0 = No, 1 = Yes) [Integer]
    # 4. Speed? (0 = max speed, 1 = speed taught in training) [Integer]
    # 5. Number of the RoMa (0 = RoMa 1, 1 = RoMa 2...) [Integer]
    # 6. Cover? (0 = cover at the source, 1 = uncover at destination) [Integer]
    # 7. Lid grid location (0, or 1-67) [String]
    # 8. Labware type (from the Evo CFG file) [String]
    # 9. Vector Name ("Narrow", "Wide", ....) [String]
    # 10. "" (unused) [String]
    # 11. "" (unused) [String]
    # 12. Carrier Name Source (name from carrier names) [String]
    # 13. Carrier Name Lid (name from carrier names or "") [String]
    # 14 .Carrier Name Destination (name from carrier names) [String]
    # 15. Source Location Carrier Site (#) [String]
    # 16. Lid Location Carrier Site (# or '(Not defined)') [String]
    # 17. Destination Location Carrier Site (#) [String]
    # Start by defining some default values

    source_location = transfer_command[2]
    source_carrier_name = current_bed[source_location[0]][source_location[2]]['carrier_name']
    destination_location = transfer_command[3]
    destination_carrier_name = current_bed[destination_location[0]][destination_location[2]]['carrier_name']
    labware_type = wellplate_library[transfer_command[4]]['labware_type']

    # Check on the initial site to make sure that it is accessible
    problem_return = check_for_problems(current_bed, source_location, labware_type)
    if problem_return != 'Acceptable':
        return ['Error', 'Unable to access the initial labware at site %s' % (source_location)]
    # Then check on the destination location as well to make sure it is accessible
    problem_return = check_for_problems(current_bed, destination_location, labware_type)
    if problem_return != 'Acceptable':
        return ['Error', 'Unable to access the destination location at site %s' % destination_location]
    
    try:
        transfer_return = r'Transfer_Rack("%s","%s",%s,%s,%s,%s,%s,"%s","%s","%s","","","%s","","%s","%s","(Not defined)","%s");' % (
                            source_location[0], destination_location[0], 0, 0, 1, 0, 0, 0, labware_type, 
                            'Narrow', source_carrier_name, destination_carrier_name, source_location[1], destination_location[1])
        #print(transfer_return)
    except:
        return ['Error', 'Something was not defined correctly for transfer labware!']
    return ['Success', transfer_return]


def special_transfer_labware(current_bed, labware_type, transfer_command):
    # Values that get sent into this expression
    # 0. Carrier Grid Position Source (1-67) [String]
    # 1. Destination Grid Position (1-67) [String]
    # 2. Move back home? (1 = true, 0 = false) [Integer]
    # 3. Lid Handling? (0 = No, 1 = Yes) [Integer]
    # 4. Speed? (0 = max speed, 1 = speed taught in training) [Integer]
    # 5. Number of the RoMa (0 = RoMa 1, 1 = RoMa 2...) [Integer]
    # 6. Cover? (0 = cover at the source, 1 = uncover at destination) [Integer]
    # 7. Lid grid location (0, or 1-67) [String]
    # 8. Labware type (from the Evo CFG file) [String]
    # 9. Vector Name ("Narrow", "Wide", ....) [String]
    # 10. "" (unused) [String]
    # 11. "" (unused) [String]
    # 12. Carrier Name Source (name from carrier names) [String]
    # 13. Carrier Name Lid (name from carrier names or "") [String]
    # 14 .Carrier Name Destination (name from carrier names) [String]
    # 15. Source Location Carrier Site (#) [String]
    # 16. Lid Location Carrier Site (# or '(Not defined)') [String]
    # 17. Destination Location Carrier Site (#) [String]
    # Start by defining some default values
    transfer_return = [None, None, 0, 0, 1, 0, 0, 0, None, "Narrow", "", "", None, "", None, None, "(Not defined)",
                       None]
    transfer_types = [1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    source_location = transfer_command[2]
    transfer_return[0], transfer_return[15], transfer_return[12] = source_location
    transfer_return[12] = current_bed[source_location[0]][source_location[2]]['carrier_name']
    destination_location = transfer_command[3]
    transfer_return[1], transfer_return[17], transfer_return[14] = destination_location
    transfer_return[14] = current_bed[destination_location[0]][destination_location[2]]['carrier_name']
    transfer_return[8] = labware_type
    transfer_string = 'Transfer_Rack('
    for i in range(0, len(transfer_return)):
        if transfer_types[i] == 1:
            string_add = '"' + str(transfer_return[i]) + '"'
        else:
            string_add = str(transfer_return[i])
        if i == len(transfer_return) - 1:
            transfer_string += string_add + ');'
        else:
            transfer_string += string_add + ','
    if None in transfer_return:
        return ['Error', 'Something was not defined correctly!']
    return ['Success', transfer_string]


def special_lid_transfer_labware(current_bed, wellplate_library, transfer_command):
    # Values that get sent into this expression
    # 0. Carrier Grid Position Source (1-67) [String]
    # 1. Destination Grid Position (1-67) [String]
    # 2. Move back home? (1 = true, 0 = false) [Integer]
    # 3. Lid Handling? (0 = No, 1 = Yes) [Integer]
    # 4. Speed? (0 = max speed, 1 = speed taught in training) [Integer]
    # 5. Number of the RoMa (0 = RoMa 1, 1 = RoMa 2...) [Integer]
    # 6. Cover? (0 = cover at the source, 1 = uncover at destination) [Integer]
    # 7. Lid grid location (0, or 1-67) [String]
    # 8. Labware type (from the Evo CFG file) [String]
    # 9. Vector Name ("Narrow", "Wide", ....) [String]
    # 10. "" (unused) [String]
    # 11. "" (unused) [String]
    # 12. Carrier Name Source (name from carrier names) [String]
    # 13. Carrier Name Lid (name from carrier names or "") [String]
    # 14 .Carrier Name Destination (name from carrier names) [String]
    # 15. Source Location Carrier Site (#) [String]
    # 16. Lid Location Carrier Site (# or '(Not defined)') [String]
    # 17. Destination Location Carrier Site (#) [String]
    # Start by defining some default values
    source_location = transfer_command[2]
    destination_location = transfer_command[3]
    lid_location = transfer_command[4]
    if transfer_command[5] == 'source':
        lid_method = 0
        carrier_name_lid = current_bed[lid_location[0]][lid_location[2]]['carrier_name']
    elif transfer_command[5] == 'destination':
        lid_method = 1
        carrier_name_lid = current_bed[lid_location[0]][lid_location[2]]['carrier_name']
    else:
        return ['Error', 'The lid location %s is not a valid location' % transfer_command[5]]
    move_back_home = 0
    lid_handling = 1
    speed = 1
    roma_number = 0
    labware_type = wellplate_library[str(transfer_command[6])]['labware_type']
    carrier_name_source = current_bed[source_location[0]][source_location[2]]['carrier_name']
    carrier_name_destination = current_bed[destination_location[0]][destination_location[2]]['carrier_name']

    transfer_return = 'Transfer_Rack("%s","%s",%s,%s,%s,%s,%s,"%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s")' \
                      % (source_location[0], destination_location[0], move_back_home, lid_handling, speed,
                         roma_number, lid_method, lid_location[0], labware_type, "Narrow", '', '', carrier_name_source,
                         carrier_name_lid, carrier_name_destination, source_location[1], lid_location[1],
                         destination_location[1])
    return ['Success', transfer_return]


def wash_tips(tipmask=255, waste_vol=str(3), cleaner_vol=str(3)):
    # Values that get added into this command
    # 1 - tipMask [Integer] --- Bit encoded tips
    # 2 - Waste grid position [Integer]
    # 3 - Waste site x position [Integer]
    # 4 - Cleaner grid position [Integer]
    # 5 - Cleaner site position [Integer]
    # 6 - Waste Volume [String]
    # 7 - Waste delay [Integer]
    # 8 - Cleaner Volume [String]
    # 9 - Cleaner Delay [Integer]
    # 10 - Air-gap [Float]
    # 11 - Air-gap Speed [Integer]
    # 12 - Retract Speed [Integer]
    # 13 - FastWash [Integer] --- 1 or 0
    # 14 - Low Volume [Integer] --- 1 or 0
    # 15 - At Frequency [Integer] --- 50-1500 Hz for Active Tips
    # 16 - Arm [Integer] --- 0 or 1
    wash_position = [1, 1]
    cleaner_position = [1, 0]
    waste_delay = 500
    cleaner_delay = 500
    airgap = [50, 70]
    retract_speed = 30
    fast_wash = 1
    low_volume = 0
    frequency = 1000
    arm = 0
    wash_return = 'Wash(%s,%s,%s,%s,%s,"%s",%s,"%s",%s,%s,%s,%s,%s,%s,%s,%s);' % (
        tipmask, *wash_position, *cleaner_position, waste_vol, waste_delay, cleaner_vol,
        cleaner_delay, *airgap, retract_speed, fast_wash, low_volume, frequency, arm)
    return ['Success', wash_return]


def mix_liha(current_bed, cl_dict, location, tip_list, well_list, volume_list,
             liquid_class='Adjusted water free dispense', cycles=5):
    # Values that get added into this command
    # 1. tipMask - Bit Encoded Integer
    # 2. liquid class - String
    # 3 - 14. Volume - expression
    # 15. Grid Site - Integer
    # 16. Carrier Site - Integer [Site - 1]
    # 17. Tip spacing
    # 18. Well Selection - Bit encoded string
    # 19. Number of mixing cycles - Integer
    # 20. Number of loop options - Integer (If non-zero 21-23 are needed)
    #       21. loopname
    #       22. action
    #       23. difference
    # 24. Arm - Number of the LiHa Arm
    
    
    vol_list = [0] * 8
    for index, tip in enumerate(tip_list):
        vol_list[tip-1] = volume_list[index]
    
    tipmask = sum([2 ** (int(tip) - 1) for tip in tip_list])
    labware_type = current_bed[location[0]][location[2]]['labware_types'][location[1] - 1]
    tipspacing = 1
    well_selection = well_name_to_string(well_list, labware_type, cl_dict)
    if well_selection[0] == 'Error':
        return ['Error', well_selection[1]]
    mix_return = 'Mix(%s,"%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s",%s,%s,%s,"%s",%s,%s,%s);' % (
        tipmask, liquid_class, *vol_list, str(0), str(0), str(0), str(0), location[0], location[1]-1, tipspacing,
        well_selection[1], cycles, str(0), str(0))
    return ['Success', mix_return]


def move_liha(location=None, travel_type=None):
    # Values that get added into this command
    # 1. Tip mask - Bit Encoded Integer
    # 2. Grid Location - Integer
    # 3. Carrier Site - Integer [Site - 1]
    # 4. Spacing - Tip Spacing
    # 5. Well Selection - Bit-Encoded Well-Selection
    # 6. zMove - Integer
    #   0 = positioning with global z-travel
    #   1 = positioning with local z-travel
    #   2 = x-move
    #   3 = y-move
    #   4 = z-move
    # 7. zTarget - Integer
    #   0 = z-travel
    #   1 = z-dispense
    #   2 = z-start
    #   3 = z-max
    #   4 = global z-travel
    # 8. offset - floating point number -1000 - 1000
    # 9. speed - floating point number 0.1 - 400
    # 10. no of loop options - integer requires 11-13 to be repeated n times
    #       11. loopname - string
    #       12. action - integer
    #       13. difference - integer
    # 14. Arm - integer of LiHa Arm
    # The default will just send the LiHa to the waste station to get it out of the way
    if location is None:
        location = [1, 2, 30]
    if travel_type is None:
        travel_type = [0, 4]
    tipmask = str(255)
    spacing = str(1)
    well_encoding = '0108¯1'
    zmove = str(travel_type[0])
    ztarget = str(travel_type[1])
    offset = str(0)
    speed = str(10)
    move_return = 'MoveLiha(%s,%s,%s,%s,"%s",%s,%s,%s,%s,%s,%s);' % (
        tipmask, str(location[0]), str(location[1] - 1), spacing,
        well_encoding, zmove, ztarget, offset, speed, str(0), str(0))
    return ['Success', move_return]


def liha_detect_liquid(current_bed, cl_dict, command, liquidclass='Adjusted Water Liquid Detection'):
    # Values that get sent into this expression
    # 1. tipMask [Integer] - Bit encoded LiHa tip expression
    # 2. liquidClass [String] - Name from the LC files from EVO
    # 3. Carrier Grid Position [Integer] - Value is between 1-67
    # 4. Carrier Site Position [Integer] - Value is Site-1
    # 5. Tip Spacing [Integer]
    # 6. WellSelection [String] - Bit coded well-selection string
    # 7. Loop Options [Integer] - 0-7 (If non-zero additional options 8-10 are needed)
    #   8. LoopName [String]
    #   9. Action [Integer] - 0=Column, 1=Row, 2=Well, 3=Rack
    #   10. Difference [Integer] - Value to vary by in action 9
    # 8/11. Arm Number [Integer]
    labware, tip_list, well_list = command
    tipmask = sum([2 ** (int(tip) - 1) for tip in tip_list])
    labware_grid = labware[0]
    labware_site = labware[1] - 1
    tip_spacing = str(1)
    labware_type = current_bed[labware_grid][labware[2]]['labware_types'][labware_site]
    well_selection = well_name_to_string(well_list, labware_type, cl_dict)
    if well_selection[0] == 'Error':
        return ['Error', well_selection[1]]
    loop_options = str(0)
    liha_arm = str(0)
    detect_return = 'Detect_Liquid(%s,"%s",%s,%s,%s,"%s",%s,%s);' % (tipmask, liquidclass,
                                                                     labware_grid, labware_site, tip_spacing,
                                                                     well_selection[1], loop_options, liha_arm)
    return ['Success', detect_return]


def liha_aspirate_dispense(current_bed, cl_dict, command, liquidclass='Adjusted water free dispense'):
    # Values that get sent into this expression
    # 1 - Tip mask [Integer] - Bit encoded tips
    # 2 - Liquid Class [String] 'DMSO free dispense'
    # 3-14 - Volume 12 values, 0 for tips not installed, strings for the rest
    # 15 - Carrier Grid Position [Integer]
    # 16 - Carrier Site Positon [Integer]
    # 17 - Tip Spacing [Integer] 1
    # 18 - [String] bit encoded well selection
    # 19 - number of loop options (0-7 required integer)
    #       20 - if 19 != 0 then this is the loopName [String]
    #       21 - if 19 != 0 then this is the action [integer]
    #       22 - if 19 != 0 then this is the difference
    # 20/23 - Arm number (0 or 1) [Integer]

    origin_labware, dest_labware, tip_list, vol_list, well_list = command
    tipmask = sum([2 ** (int(tip) - 1) for tip in tip_list])
    volume = [0] * 12
    for i in range(0, len(tip_list)):
        volume[tip_list[i] - 1] = vol_list[i]
    volume_string = ','.join(['"%s"' % '{:.1f}'.format(el) for el in volume[0:8]]) + ',' + ','.join(
        '%s' % str(el) for el in volume[8:])
    asp_grid = origin_labware[0]
    asp_site = origin_labware[1] - 1
    tip_spacing = str(1)
    labware_type = current_bed[asp_grid][origin_labware[2]]['labware_types'][asp_site]
    if origin_labware[3] == 'XX':
        asp_tips = [chr(ord('A') + el - 1) + '1' for el in tip_list]
    else:
        asp_tips = origin_labware[3]
    well_selection = well_name_to_string(asp_tips, labware_type, cl_dict)
    if well_selection[0] == 'Error':
        return ['Error', well_selection[1], '']
    loops = str(0)
    liha_arm = str(0)
    aspirate_return = 'Aspirate(%s,"%s",%s,%s,%s,%s,"%s",%s,%s);' % (
        str(tipmask), liquidclass, volume_string, str(asp_grid), str(asp_site),
        tip_spacing, well_selection[1], loops, liha_arm)
    dest_grid = dest_labware[0]
    dest_site = dest_labware[1] - 1
    dest_labware_type = current_bed[dest_grid][dest_labware[2]]['labware_types'][dest_site]
    well_selection = well_name_to_string(well_list, dest_labware_type, cl_dict)
    if well_selection[0] == 'Error':
        return ['Error', well_selection[1], '']

    dispense_return = 'Dispense(%s,"%s",%s,%s,%s,%s,"%s",%s,%s);' % (
        str(tipmask), liquidclass, volume_string, str(dest_grid), str(dest_site),
        tip_spacing, well_selection[1], loops, liha_arm)
    return ['Success', aspirate_return, dispense_return]


def liha_aspirate(current_bed, cl_dict, command, liquidclass='Adjusted water free dispense'):
    # Values that get sent into this expression
    # 1 - Tip mask [Integer] - Bit encoded tips
    # 2 - Liquid Class [String] 'DMSO free dispense'
    # 3-14 - Volume 12 values, 0 for tips not installed, strings for the rest
    # 15 - Carrier Grid Position [Integer]
    # 16 - Carrier Site Positon [Integer]
    # 17 - Tip Spacing [Integer] 1
    # 18 - [String] bit encoded well selection
    # 19 - number of loop options (0-7 required integer)
    #       20 - if 19 != 0 then this is the loopName [String]
    #       21 - if 19 != 0 then this is the action [integer]
    #       22 - if 19 != 0 then this is the difference
    # 20/23 - Arm number (0 or 1) [Integer]

    origin_labware, tip_list, vol_list = command
    tipmask = sum([2 ** (int(tip) - 1) for tip in tip_list])
    volume = [0] * 12
    for i in range(0, len(tip_list)):
        volume[tip_list[i] - 1] = vol_list[i]
    volume_string = ','.join(['"%s"' % '{:.1f}'.format(el) for el in volume[0:8]]) + ',' + ','.join(
        '%s' % str(el) for el in volume[8:])
    asp_grid = origin_labware[0]
    asp_site = origin_labware[1] - 1
    tip_spacing = str(1)
    labware_type = current_bed[asp_grid][origin_labware[2]]['labware_types'][asp_site]
    if origin_labware[3] == 'XX':
        asp_tips = [chr(ord('A') + el - 1) + '1' for el in tip_list]
    else:
        asp_tips = origin_labware[3]
    well_selection = well_name_to_string(asp_tips, labware_type, cl_dict)
    if well_selection[0] == 'Error':
        return ['Error', well_selection[1], '']
    loops = str(0)
    liha_arm = str(0)
    aspirate_return = 'Aspirate(%s,"%s",%s,%s,%s,%s,"%s",%s,%s);' % (
        str(tipmask), liquidclass, volume_string, str(asp_grid), str(asp_site),
        tip_spacing, well_selection[1], loops, liha_arm)
    return ['Success', aspirate_return]


def liha_dispense(current_bed, cl_dict, command, liquidclass='Adjusted water free dispense'):
    # Values that get sent into this expression
    # 1 - Tip mask [Integer] - Bit encoded tips
    # 2 - Liquid Class [String] 'DMSO free dispense'
    # 3-14 - Volume 12 values, 0 for tips not installed, strings for the rest
    # 15 - Carrier Grid Position [Integer]
    # 16 - Carrier Site Positon [Integer]
    # 17 - Tip Spacing [Integer] 1
    # 18 - [String] bit encoded well selection
    # 19 - number of loop options (0-7 required integer)
    #       20 - if 19 != 0 then this is the loopName [String]
    #       21 - if 19 != 0 then this is the action [integer]
    #       22 - if 19 != 0 then this is the difference
    # 20/23 - Arm number (0 or 1) [Integer]

    dest_labware, tip_list, vol_list, well_list = command
    tipmask = sum([2 ** (int(tip) - 1) for tip in tip_list])
    volume = [0] * 12
    for i in range(0, len(tip_list)):
        volume[tip_list[i] - 1] = vol_list[i]
    volume_string = ','.join(['"%s"' % '{:.1f}'.format(el) for el in volume[0:8]]) + ',' + ','.join(
        '%s' % str(el) for el in volume[8:])
    tip_spacing = str(1)
    loops = str(0)
    liha_arm = str(0)

    dest_grid = dest_labware[0]
    dest_site = dest_labware[1] - 1
    dest_labware_type = current_bed[dest_grid][dest_labware[2]]['labware_types'][dest_site]
    well_selection = well_name_to_string(well_list, dest_labware_type, cl_dict)
    if well_selection[0] == 'Error':
        return ['Error', well_selection[1], '']
    dispense_return = 'Dispense(%s,"%s",%s,%s,%s,%s,"%s",%s,%s);' % (
        str(tipmask), liquidclass, volume_string, str(dest_grid), str(dest_site),
        tip_spacing, well_selection[1], loops, liha_arm)
    return ['Success', dispense_return]


def start_heater_shaker(location, rpms, shaking_type='0 Circle anticlockwise'):
    # Values that are needed in this string
    # 1 - Name of the thing [String] --- inhecoMTC
    # 2 - Command Name [String] --- inhecoMTC_StartShake
    # 3 - Module command string
    #   3a - Location of the heater shaker [Integer]
    #   3b - RPMs [Integer]
    #   3c - Type of Rotation [String]
    # 4 - Not sure.... [String]
    # 5 - Unused

    # Start by doing some quick sanity checks!
    # RPMs are bounded by [100, 2000 RPM]
    lower_limit = 100
    upper_limit = 2000
    try:
        if float(rpms) < lower_limit or float(rpms) > upper_limit:
            return ['Error', 'The requested RPMs (%s) is outside of the valid limits %s -> %s' % (rpms, lower_limit, upper_limit)]
    except TypeError:
        return ['Error', 'The requested RPM input (%s) cannot be converted to a numerical value' % rpms]

    location_id = {'323': {1: 5, 2: 6},
                   '322': {1: 3, 2: 4},
                   '324': {1: 1, 2: 2}}
    inheco_id = location_id[str(location[2])][location[1]]
    start_return = 'FACTS("inhecoMTC","inhecoMTC_StartShake","%s","0","");' % (
        ','.join([str(inheco_id), str(rpms), shaking_type]))
    return ['Success', start_return]


def set_heater_shaker_temperature(location, operation='off', temperature=20):
    # Values that are needed in this string
    # 1 - Name of the thing [String] --- inhecoMTC
    # 2 - Command Name [String] --- inhecoMTC_SetTemperature
    # 3 - Module command string
    #   3a - Location of the heater shaker [Integer]
    #   3b - Turn on or off --- true or flase [String]
    #   3c - Set temperature * 10 [String]
    #   3d - Time to wait for [String]
    # 4 - Not sure... [String]
    # 5 - Unused

    # In the event that the heater-shaker was not found by means of the find
    # open heater shaker command we will check the temperature bounds
    # Low T Heater Shaker bounds: +4C to +80C
    # High T Heater Shaker bounds: 5C above ambient (20C?) to +125C
    lt_limits = (4, 80)
    lt_heatshake_ids = ['322', '324']
    ht_limits = (20, 125)
    ht_heatshake_ids = ['323']
    if str(location[2]) in lt_heatshake_ids:
        if temperature not in range(*lt_limits):
            return ['Error', 'Issue setting the temperature (%s) for location %s, limits are %s -> %s' %
                    (temperature, location, lt_limits[0], lt_limits[1])]
    elif str(location[2]) in ht_heatshake_ids:
        if temperature not in range(*ht_limits):
            return ['Error', 'Issue setting the temperature (%s) for location %s, limits are %s -> %s' %
                    (temperature, location, ht_limits[0], ht_limits[1])]
    else:
        return ['Error', 'Really bad issue that terminated in unknown error, check inputs for location %s' % location]

    location_id = {'323': {1: 5, 2: 6},
                   '322': {1: 3, 2: 4},
                   '324': {1: 1, 2: 2}}
    try:
        inheco_id = location_id[str(location[2])][location[1]]
    except KeyError:
        return ['Error', 'Issue locating the heater-shaker location at %s' % location]

    if operation == 'on':
        operation_string = 'True'
    else:
        operation_string = 'False'
    temperature_set = '{:.0f}'.format(temperature * 10)
    temperature_return = 'FACTS("inhecoMTC","inhecoMTC_SetTemperature","%s","0","");' % ','.join(
        [str(inheco_id), operation_string, temperature_set, '0'])
    return ['Success', temperature_return]


def stop_heater_shaker(location):
    # Values that are needed in this string
    # 1 - Name of the thing [String] --- inhecoMTC
    # 2 - Command Name [String] --- inhecoMTC_StopShake
    # 3 - Module command string
    #   3a - Location of the heater shaker [Integer]
    #   3b - Turn the thing off
    # 4 - Not sure probably not used
    location_id = {'323': {1: 5, 2: 6},
                   '322': {1: 3, 2: 4},
                   '324': {1: 1, 2: 2}}
    inheco_id = location_id[str(location[2])][location[1]]
    stop_return = 'FACTS("inhecoMTC","inhecoMTC_StopShake","%s","0","");' % str(inheco_id)
    return ['Success', stop_return]


def start_timer(timer_number):
    # Values that are included in this string
    # Internal Timer value [String between 1 and 100]
    if 1 > int(timer_number) or int(timer_number) > 100:
        return ['Error' , '%s is not a valid timer id, needs to be between 1 and 100' % str(timer_number)]
    timer_return = 'StartTimer("%s");' % str(timer_number)
    return ['Success', timer_return]


def wait_for_timer(timer_number, time_to_wait):
    # Values that are included in this string
    # Internal Timer Value [String between 1 and 100]
    # Time Span [String between 0.02 and 86400]
    if 1 > int(timer_number) or int(timer_number) > 100:
        return ['Error', '%s is not a valid timer id, needs to be between 1 and 100' % str(timer_number)]
    if float(time_to_wait) < 0.02 or float(time_to_wait) > 86400:
        return ['Error', '%s is not a valid time to wait, needs to be between 0.02 and 86400' % str(time_to_wait)]
    wait_timer_return = 'WaitTimer("%s","%s");' % (str(timer_number), str(time_to_wait))
    return ['Success', wait_timer_return]


def comment(comment_string):
    # Values that are included in this string
    # String to become comment string...
    return ['Success', 'Comment("%s");' % str(comment_string)]


def set_variable(variable_name, variable_value, variable_type, variable_scope):
    # Values that are included in this string
    # 1. Name of the variable (string2)
    # 2. Value assigned to variable (expression)
    # 3. queryFlag (integer; 1 for user, 0 for no user query)
    # 4. queryString (string; message for user query)
    # 5. checkLimits (integer; 1 to validate user entry, 0 to not)
    # 6. lowerLimit (float; limit for user validate)
    # 7. upperLimit (float; limit for user validate)
    # 8. Type (integer; 0 = Numeric, 1 = String)
    # 9. Scope (integer; 0 = run, 1 = instance, 2 = script)
    # 10. initMode (integer; 0 = fixed value, 1 = user query, 2 = file import)
    # 11. QueryAtStart (integer; 1 = prompt at start, 0 to not)
    queryflag = 0
    querystring = ''
    checklimits = 0
    lowerlimit = '1.000000'
    upperlimit = '10.000000'
    initmode = 0
    queryatstart = 0
    set_variable_return = 'Variable(%s,"%s",%s,"%s",%s,%s,%s,%s,%s,%s,%s);' % (variable_name, variable_value, queryflag,
                                                                               querystring, checklimits, lowerlimit,
                                                                               upperlimit, variable_type,
                                                                               variable_scope, initmode, queryatstart)
    return ['Success', set_variable_return]


def export_variable(variable_names, output_filepath, write_header='1', overwrite='1'):
    # Values that are included in this string
    # 1. Variable names to export (list -> %s#%s)
    # 2. Output filepath
    # 3. Write Header? (1 = yes, 0 = no)
    # 4. Overwrite file? (1 = yes, 0 = no)
    export_variable_return = 'ExportVariable(%s,"%s",%s,%s);' % ('#'.join(variable_names),
                                                                 os.path.abspath(output_filepath),
                                                                 write_header, overwrite)
    return ['Success', export_variable_return]


def notification(sshot_bool, recieve_group, subject, message, action):
    # Values that are included in this string
    # 1. Attach screenshot to email? [integer 0 or 1]
    # 2. Group to Receive the email message [String, groups made in configuaration tool]
    # 3. Email Subject [String]
    # 4. Email Message [String]
    # 5. Action [Integer, 0=Now, 1=On Error, 2=Stop Sending Email on Error]
    notification_return = 'Notification(%s,"%s","%s","%s",%s);' % (str(sshot_bool), str(recieve_group),
                                                                   str(subject), str(message), str(action))
    return ['Success', notification_return]


def execute_application(application, responsevariable, scope, options=2):
    # Values that are included in this string
    # 1. Application [String] --- Path and name of the executable
    # 2. Options [Integer] --- 1) Release COM Port, 2) Wait for termination of app, 4) Store and return in variable
    # 3. ResponseVariable [String] --- Name of the return variable
    # 4. Scope [Integer] --- 0) Run, 1) Instance, 2) Script
    exec_app_return = 'Execute("%s",%s,"%s",%s);' % (str(application), str(options), str(responsevariable), str(scope))
    return ['Success', exec_app_return]


def execute_vbs_script(vbs_script_name, options=0):
    # Values that are included in this string
    # 1. Filename [String] --- Path to the VBS script to execute
    # 2. Action [Integer] --- 0) wait for end of vbs, 1) script continues, 2) script waits for previous vbs script
    # EXECUTE_VBSCRIPT("C:\Scripts\datetime_output.vbs",0);
    execute_vbs_return = 'EXECUTE_VBSCRIPT("%s",%s);' % (vbs_script_name, options)
    return ['Success', execute_vbs_return]


def user_prompt(message, sound=0, closetime=-1):
    # Values that are included in this string
    # 1. Text [String]
    # 2. Sound - 0=No Sound, 1=Beep Once, 2=Beep Three Times, 3=Beep Every 3 Seconds
    # 3. Close Time [integer, time when prompt closes, -1=Never]
    user_prompt_return = 'UserPrompt("%s",%s,%s);' % (str(message), sound, closetime)
    return ['Success', user_prompt_return]


def on_error_go_to(error_handle_bool, line):
    # Values that are included in this string
    # 1. Activate Error Handling? (0 = activate, 1 = deactivate)
    # 2. Scripting line to go to
    on_error_return = 'OnError(%s,%s);' % (error_handle_bool, line)
    return ['Success', on_error_return]


def find_mca_tips(current_bed, library, tip_type='DiTi_200uL_Nested'):
    if tip_type == 'DiTi_200uL_Nested':
        carrier_to_find = 328
        pass
    else:
        return ['Error', '%s tip type is not yet implemented correctly','']

    for grid in current_bed.keys():
        for carrier in current_bed[grid].keys():
            if carrier == carrier_to_find:
                return_location = [grid, '', carrier]
                break
        else:
            continue
        break
    else:
        return ['Error', '%s the specified carrier was not found','']

    for consumable in library['consumables'].keys():
        #pprint(library[consumable])
        
        if library['consumables'][consumable]['labware_type'] == tip_type:
            if library['consumables'][consumable]['location'][0] != 'liquid_handler':
                continue
            if library['consumables'][consumable]['location'][1] != [carrier_to_find]:
                continue
            if not library['consumables'][consumable]['tip_locations']:
                return ['Error', 'No tips available','']
            else:
                return_location[1] = int(library['consumables'][consumable]['tip_locations'][0])
                break
    else:
        return ['Error', 'Something is incorrectly defined in the database','']
    return ['Success', return_location, consumable]


def mca_get_tips(tips_location):
    # Values that are included in this string
    # 1. Carrier grid position of tips [Integer]
    # 2. Carrier site position of tips [Integer - 1]
    # 3. Blowout Airgap [String formatted float number Ex. 1.000000 uL]
    # 4. Number of rows [Integer]
    # 5. Number of columns [Integer]
    # 6. There is something else that goes here... [0 for now integer]
    # 7. There is something else that goes here... [0 for now integer]
    airgap_vol = '{:.6f}'.format(float(0))
    nrows = 8
    ncols = 12
    get_ditis_return = 'MCAGetDitis(%s,%s,%s,%s,%s,%s,%s);' % (str(tips_location[0]), str(tips_location[1] - 1),
                                                               airgap_vol, nrows, ncols, str(0), str(0))
    return ['Success', get_ditis_return]


def mca_aspirate(location, volume, liquid_class='Adjusted Water free dispense DiTi 200 Clear', compartment=1):
    # Values that are included in this string
    # 1. Liquid Class [String]
    # 2. Volume [String]
    # 3. Carrier grid position [Integer]
    # 4. Carrier site position [Integer - 1]
    # 5. Bit encoded well-selection (Also need to know the compartment that this happens in)
    # 6. There is something else that goes here... [0 for now integer]
    # 7. There is something else that goes here... [0 for now integer]
    if int(compartment) == 1:
        bit_encoded_tips = '0C08¯¯¯¯¯¯¯¯¯¯¯¯¯O'
    elif int(compartment) == 2:
        bit_encoded_tips = '0C080000000000000¯¯¯¯¯¯¯¯¯¯¯¯¯7'
    else:
        return ['Error', 'Compartment %s has not been included yet (only 1 and 2), need to know bit encoded tips' % str(compartment)]

    mca_aspirate_return = 'MCAAspirate("%s","%s",%s,%s,"%s",%s,%s);' % (liquid_class, str(volume),
                                                                        str(location[0]), str(location[1] - 1),
                                                                        bit_encoded_tips, str(0), str(0))
    return ['Success', mca_aspirate_return]


def mca_dispense(location, volume, liquid_class='Adjusted Water free dispense DiTi 200 Clear', compartment=1):
    # Values that are included in this string
    # 1. Liquid Class [String]
    # 2. Volume [String]
    # 3. Carrier grid position [Integer]
    # 4. Carrier site position [Integer - 1]
    # 5. Bit encoded well-selection (Also need to know the compartment that this happens in)
    # 6. There is something else that goes here... [0 for now integer]
    # 7. There is something else that goes here... [0 for now integer]
    if int(compartment) == 1:
        bit_encoded_tips = '0C08¯¯¯¯¯¯¯¯¯¯¯¯¯O'
    elif int(compartment) == 2:
        bit_encoded_tips = '0C080000000000000¯¯¯¯¯¯¯¯¯¯¯¯¯7'
    else:
        return ['Error', 'Compartment %s has not been included yet (only 1 and 2), need to know bit encoded tips' % str(compartment)]

    mca_dispense_return = 'MCADispense("%s","%s",%s,%s,"%s",%s,%s);' % (liquid_class, str(volume),
                                                                        str(location[0]), str(location[1] - 1),
                                                                        bit_encoded_tips, str(0), str(0))
    return ['Success', mca_dispense_return]


def mca_mix(location, volume, cycles=20, liquid_class='Adjusted Water free dispense DiTi 200 Clear Fast', compartment=1):
    # Values that are included in this string
    # 1. Liquid Class [String]
    # 2. Volume [String]
    # 3. Carrier Grid Position [Integer]
    # 4. Carrier Site Position [Integer - 1]
    # 5. Bit encoded well-selection 
    # 6. Number of cycles
    # 7. There is something else that goes here... [0 for now integer]
    # 8. There is something else that goes here... [0 for now integer]
    if int(compartment) == 1:
        bit_encoded_tips = '0C08¯¯¯¯¯¯¯¯¯¯¯¯¯O'
    elif int(compartment) == 2:
        bit_encoded_tips = '0C080000000000000¯¯¯¯¯¯¯¯¯¯¯¯¯7'
    else:
        return ['Error', 'Compartment %s has not been included yet (only 1 and 2), need to know bit encoded tips' % str(compartment)]
    mca_mix_return = 'MCAMix("%s","%s",%s,%s,"%s",%s,%s,%s);' % (liquid_class, str(volume),
                                                                 str(location[0]), str(location[1] - 1),
                                                                 bit_encoded_tips, str(cycles), str(0), str(0))
    return ['Success', mca_mix_return]


def mca_droptips(current_bed, position_drop=1, waste_site=17):
    # Values that are included in this string
    # 1. Carrier grid position [Integer]
    # 2. Carrier site position [Integer]
    # 3. Return tips to source or to position [0 is source, 1 is position]
    # 4. There is something else that goes here... [0 for now integer]
    # 5. There is something else that goes here... [0 for now integer]
    # 6. There is something else that goes here... [0 for now integer]
    carrier_to_find = 328
    for grid in current_bed:
        for carrier in current_bed[grid]:
            if carrier == carrier_to_find:
                carrier_grid = grid
                break
        else:
            continue
        break
    else:
        return ['Error']
    tip_waste_location = [carrier_grid, waste_site]
    mca_droptips_return = 'MCADropDitis(%s,%s,%s,%s,%s,%s);' % (
        str(tip_waste_location[0]), str(tip_waste_location[1] - 1),
        str(position_drop), str(0), str(0), str(0))
    return ['Success', mca_droptips_return]


def mca_move(location=None, travel_type=None):
    # Values that get added into this command
    # 1. Tip mask - Bit Encoded Integer [weirdly defined...]
    # 2. Grid Location - Integer
    # 3. Carrier Site - Integer [Site - 1]
    # 4. Spacing - Tip Spacing
    # 5. Well Selection - Bit-Encoded Well-Selection [basically always one value so...]
    # 6. zMove - Integer
    #   0 = positioning with global z-travel
    #   1 = positioning with local z-travel
    #   2 = x-move
    #   3 = y-move
    #   4 = z-move
    # 7. zTarget - Integer
    #   0 = z-travel
    #   1 = z-dispense
    #   2 = z-start
    #   3 = z-max
    #   4 = global z-travel
    # 8. offset - floating point number -1000 - 1000
    # 9. speed - floating point number 0.1 - 400
    # 10. no of loop options - integer requires 11-13 to be repeated n times
    #       11. loopname - string
    #       12. action - integer
    #       13. difference - integer
    # 14. Arm - integer of LiHa Arm
    # The default will just send the LiHa to the waste station to get it out of the way
    if location is None:
        location = [48, 17, 328]
    if travel_type is None:
        travel_type = [0, 4]
    tipmask = str(255)
    spacing = str(1)
    well_encoding = r'0C08¯¯¯¯¯¯¯¯¯¯¯¯¯O'
    zmove = str(travel_type[0])
    ztarget = str(travel_type[1])
    offset = str(0)
    speed = str(10)
    move_return = 'MCAMove(%s,%s,%s,%s,"%s",%s,%s,%s,%s,%s,%s);' % (
        tipmask, str(location[0]), str(location[1] - 1), spacing,
        well_encoding, zmove, ztarget, offset, speed, str(0), str(0))
    return ['Success', move_return]


def roma_vector(vector_name, location, direction, command_dict):
    # Values that get added into this command
    # 1. Vector Name - Defined in EvoWare, annoying to lookup....
    # 2. Grid Location - Integer
    # 3. Carrier Site - Integer [Site - 1]
    # 4. Direction - 0 safe to end, 1 end to safe
    # 5. Back? - 1 does round trip, 0 is oneway
    # 6. begin action
    # 7. end action
    # 8. speed
    # 9. romaNo

    if 'back' in command_dict.keys():
        back = command_dict['back']
    else:
        back = 0 
    if 'endAction' in command_dict.keys():
        endAction = command_dict['endAction']
    else:
        endAction = 2
    if 'beginAction' in command_dict.keys():
        beginAction = command_dict['beginAction']
    else:
        beginAction = 2
    if 'speed' in command_dict.keys():
        speed = command_dict['speed']
    else:
        speed = 0
    romaNo = 0
    
    move_return = 'Vector("%s","%s","%s",%s,%s,%s,%s,%s,%s)' % (vector_name, location[0],
                                                                location[1], direction,
                                                                back, beginAction,
                                                                endAction, speed,  romaNo)
    return ['Success', move_return]


def restock_mca_tips(tip_location, current_bed):
    # These are the two important SBS base sites on the nested tip carrier
    # base_site = [25, 26]
    # Find the waste chute location on the platform bed
    # carrier_to_find = 331
    # waste_site = 1
    # for grid in current_bed:
    #    for carrier in current_bed[grid]:
    #        if carrier == carrier_to_find:
    #            carrier_grid = grid
    #            break
    #    else:
    #        continue
    #    break
    # waste_location = [carrier_grid, waste_site, carrier_to_find]
    # restock_return = []
    # for site in base_site:
    #    transfer_command = ['TRANSFER_LABWARE', duplicate_details['Well_Plate_ID'], duplicate_plate_initial_location,
    #                        child_well_plate_home]
    #    method_strings.append(transfer_labware(updated_carriers, full_library['well_plates'], transfer_command))

    return ['Error', 'Not Implemented']


def large_waste_dump(starting_location, current_bed, labware_type='DiTi 200ul Nested MCA96'):
    # This is basically a transfer_rack command except there is only one destination
    # for now, this can be updated if needed later
    # Values that get sent into this expression
    # 0. Carrier Grid Position Source (1-67) [String]
    # 1. Destination Grid Position (1-67) [String]
    # 2. Move back home? (1 = true, 0 = false) [Integer]
    # 3. Lid Handling? (0 = No, 1 = Yes) [Integer]
    # 4. Speed? (0 = max speed, 1 = speed taught in training) [Integer]
    # 5. Number of the RoMa (0 = RoMa 1, 1 = RoMa 2...) [Integer]
    # 6. Cover? (0 = cover at the source, 1 = uncover at destination) [Integer]
    # 7. Lid grid location (0, or 1-67) [String]
    # 8. Labware type (from the Evo CFG file) [String]
    # 9. Vector Name ("Narrow", "Wide", ....) [String]
    # 10. "" (unused) [String]
    # 11. "" (unused) [String]
    # 12. Carrier Name Source (name from carrier names) [String]
    # 13. Carrier Name Lib (name from carrier names or "") [String]
    # 14 .Carrier Name Destination (name from carrier names) [String]
    # 15. Source Location Carrier Site (#) [String]
    # 16. Lid Location Carrier Site (# or '(Not defined)') [String]
    # 17. Destination Location Carrier Site (#) [String]
    # Start by defining some default values   
    carrier_to_find = 331
    waste_site = 1
    for grid in current_bed:
        for carrier in current_bed[grid]:
            if carrier == carrier_to_find:
                carrier_grid = grid
                break
        else:
            continue
        break
    else:
        return 'Error'
    source = starting_location
    source_carrier_name = current_bed[source[0]][source[2]]['carrier_name']
    destination = [carrier_grid, waste_site, carrier_to_find]
    destination_carrier_name = current_bed[destination[0]][destination[2]]['carrier_name']
    vector_name = "Narrow"
    waste_string = 'Transfer_Rack("%s","%s",%s,%s,%s,%s,%s,"%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s");' % (
        str(source[0]), str(destination[0]), str(0), str(0), str(1), str(0), str(0), '', str(labware_type),
        str(vector_name), '', '', str(source_carrier_name), '', str(destination_carrier_name), str(source[1]),
        str('(Not defined)'), str(destination[1]))
    return ['Success', waste_string]


def tevac_string_return(input_command, extra_input):
    if input_command[0] == 'vacuum_on':
        if input_command[1] == 'front':
            tevac_command = 'Te-VacS_ApplyVacuumFront'
        elif input_command[1] == 'back':
            tevac_command = 'Te-VacS_ApplyVacuumRear'
        else:
            return ['Error', 'Location %s is not recognized as a TeVac location' % input_command[1]]
    elif input_command[0] == 'vent':
        if input_command[1] == 'front':
            tevac_command = 'Te-VacS_VentFront'
        elif input_command[1] == 'back':
            tevac_command = 'Te-VacS_VentRear'
        else:
            return ['Error', 'Location %s is not recognized as a TeVac location' % input_command[1]]
    elif input_command[0] == 'deactivate':
        tevac_command = 'Te-VacS_DeactivateSystem'
    else:
        return ['Error', 'Input command not recognized %s' % input_command[0]]
    return_command = 'FACTS("Te-VacS","%s",%s,"0","");' % (tevac_command, extra_input)
    return ['Success', return_command]

# import json, sys, copy
# from worktable_cleaner import initial_worktable_prep, worktable_export
# carriers_labware_filepath = r'.\lh_dictionaries\carriers_and_labware.json'
# empty_worktable_filepath = r'.\current_platform\empty_layout_20200821_final.ewt'
# updated_library_path = r'.\current_platform\updated_library_test.json'
# with open(updated_library_path, 'r') as jsonfile:
#    full_library = json.load(jsonfile)
# with open(carriers_labware_filepath, 'r') as jsonfile:
#    carriers_labware = json.load(jsonfile)
# start_carriers = initial_worktable_prep(carriers_labware_filepath, empty_worktable_filepath, updated_library_path)
# updated_carriers = copy.deepcopy(start_carriers)
#
# returned_val = large_waste_dump([48, 1, 328], updated_carriers)
#
# print(returned_val)
