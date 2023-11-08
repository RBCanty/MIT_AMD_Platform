# -*- coding: utf-8 -*-
"""
Updated on 11/02/2021 --- MongoDB containing version of prepare internal standard

@author: Brent
"""
import copy
import datetime
import json
import math
import traceback
import sys
import pymongo
import yaml
from pprint import pprint

if __name__ == "__main__":
    from platform_library_functions import query_collection, update_database, query_document, \
        update_database_document, move_database_document
    from tecan_functions import mca_move, find_mca_tips, tevac_string_return, start_timer, wait_for_timer, \
        mca_droptips, mca_get_tips, mca_mix, mca_aspirate, mca_dispense
    from tecan_functions import transfer_labware, wash_tips, liha_aspirate_dispense, \
        move_liha, special_transfer_labware, special_lid_transfer_labware, large_waste_dump
    from useful_other_functions import solvent_rinse_location, solvent_finder, \
        liha_grouping, find_accessible_open_location, directory_check, chemical_well_origin, priority_sort
    from worktable_cleaner import initial_worktable_prep, worktable_export
else:
    from Evoware_API.platform_library_functions import query_collection, update_database, query_document, \
        update_database_document, move_database_document
    from Evoware_API.tecan_functions import mca_move, find_mca_tips, tevac_string_return, start_timer, wait_for_timer, \
        mca_droptips, mca_get_tips, mca_mix, mca_aspirate, mca_dispense
    from Evoware_API.tecan_functions import transfer_labware, wash_tips, liha_aspirate_dispense, \
        move_liha, special_transfer_labware, special_lid_transfer_labware, large_waste_dump
    from Evoware_API.useful_other_functions import solvent_rinse_location, solvent_finder, \
        liha_grouping, find_accessible_open_location, directory_check, chemical_well_origin, priority_sort
    from Evoware_API.worktable_cleaner import initial_worktable_prep, worktable_export


def prepare_standard_wells(queue_name, operation_number):
    # ----------------------------------------------------------------------------
    # ----------------------------------------------------------------------------
    # Open the configuration file before starting anything
    success_statements = []
    config_filepath = r'C:\Users\admin\PycharmProjects\AMD_Control_Platform\venv\Shell_AH\Evoware_API\config.yaml'
    try:
        with open(config_filepath, 'r') as yamlfile:
            config = yaml.load(yamlfile, Loader=yaml.Loader)
    except:
        return ['Error', 'Problem loading the configuration file %s' % config_filepath, success_statements]
    
    # Grab the important paths from the configuration file
    try:
        carriers_labware_filepath = config['paths/files']['carriers_labware_path']
        empty_worktable_filepath = config['paths/files']['empty_worktable_path']
        username = config['mongo_details']['mongo_user']
        password = config['mongo_details']['mongo_pass']
        if config['startup']['mongo_testing']:
            platform_mongodb_address = config['mongo_details']['mongo_address_testing']
            platform_mongodb_port = config['mongo_details']['mongo_port_testing']
        else:
            platform_mongodb_address = config['mongo_details']['mongo_address']
            platform_mongodb_port = config['mongo_details']['mongo_port']
    except Exception as e:
        return ['Error', str(e), success_statements]

    # Start by loading the specific queue element name and check whether it is real
    try:
        platform_mongodb = pymongo.MongoClient(platform_mongodb_address, platform_mongodb_port)
    except Exception as e:
        return ['Error', str(e), success_statements]
    queue_document, queue_return_statements = query_document('queue', queue_name, 'queue_name', platform_mongodb)
    if type(queue_document) == str:
        platform_mongodb.close()
        return ['Error', queue_return_statements, success_statements]
    else:
        success_statements.append('Found element %s in the queue collection' % queue_name)

    # Now we check to make sure that the plate is ready to be prepared
    try:
        plate_operations = queue_document['operations']
        queue_containers = queue_document['containers']
        target_operation = 'prepare_standard_wells'
        if operation_number == '1':
            if plate_operations[operation_number]['completed'] != 'no':
                success_statements.append('Operation %s is not ready to be executed, status %s' % (
                                           operation_number, plate_operations[operation_number]['completed']))
                return ['Error', 'Queue completion error', success_statements]
        else:
            try:
                previous_operation = str(int(operation_number) - 1)
            except Exception as e:
                success_statements.append(['Operation number key is invalid', e])
                return ['Error', success_statements]
            if plate_operations[previous_operation]['completed'] == 'error':
                return ['Error', 'Operation %s has an error status, not ready for operation %s' % (previous_operation, operation_number), success_statements]
            elif plate_operations[previous_operation]['completed'] == 'no':
                return ['Error','Previous operation %s not marked as complete in database, not ready for operation %s' % (
                                                        previous_operation, operation_number), success_statements]
            elif plate_operations[previous_operation]['completed'] == 'yes' and \
                    plate_operations[operation_number]['completed'] == 'no':
                success_statements.append('Operation %s is in valid position in queue document' % operation_number)
            else:
                return ['Error', 'The operation %s is already marked as complete' % operation_number, success_statements]
        if plate_operations[operation_number]['operation'] != target_operation:
            return ['Error','Expected operation %s does not match queue document operation %s' % (
                                        plate_operations[operation_number]['operation'], target_operation), success_statements]
        else:
            success_statements.append('Expected operation matches queue document operation %s' % target_operation)
        additional_prep_details = plate_operations[operation_number]['details']
        target_container = plate_operations[operation_number]['container']
        wellplate_name = queue_containers[target_container]['container_name']
        target_well_plate_type = queue_containers[target_container]['plate_type']
        success_statements.append('The %s operation is in a valid position and ready to be executed' % target_operation)
    except Exception as e:
        return ['Error', str(e), success_statements]
    # TODO check to make sure the required details provided are there or provide default values (??)

    # Here there is an optional statement from config to send to the testing folder
    try:
        current_time = datetime.datetime.now()
        if config['startup']['testing']:
            schedule_path_output = config['paths/files']['schedule_testing_output_path'] % current_time.strftime('%Y%m%d')
        else:
            folder_outpath = config['paths/files']['schedule_output_path'] % current_time.strftime('%Y%m%d')
            directory_check(folder_outpath)
            schedule_path_output = folder_outpath + r'\%s_%s_%s_%s.esc' % (
                queue_name, operation_number, 'plate_prep', current_time.strftime('%Y%m%d%H%M%S'))
    except Exception as e:
        return ['Error', str(e), success_statements]
    
    # Then load all of the documents from the platform library which includes empty
    # wellplates, reagents, solvents, and consumables. This could probably be made
    # more efficient later when it is needed but for now it will be fast enough
    full_library_mongo = {}
    important_collections = ['wellplates', 'solvents', 'reagents', 'consumables']
    for important_collection in important_collections:
        message = {'auth': {'port': '27017', 'user': username, 'password': password},
                   'request': {'request_type': 'query_collection', 'collection': important_collection}}
        document_return, extra_statements = query_collection(message['request']['collection'], platform_mongodb)
        if document_return == 'Error':
            return ['Error', extra_statements, success_statements]
        full_library_mongo[important_collection] = document_return
    initial_full_library = copy.deepcopy(full_library_mongo)
    platform_mongodb.close()
    
    # Then build the current worktable dictionary for building output files, we are
    # also making a duplicate of the worktable so it can be updated in this code
    try:
        with open(carriers_labware_filepath, 'r') as jsonfile:
            carriers_labware = json.load(jsonfile)
    except Exception:
        return ['Error', 'Carriers and labware file not accessible', success_statements]
    return_statement, start_carriers = initial_worktable_prep(carriers_labware_filepath, empty_worktable_filepath, full_library_mongo)
    if return_statement != 'Success':
        return ['Error', start_carriers, success_statements]
    updated_carriers = copy.deepcopy(start_carriers)
    success_statements.append('Loaded the carrier and labware library file!')

    # We start by finding the initial wellplate that we want to mirror for the standards wells
    collections_to_check = ['wellplates', 'reagents']
    for collection in collections_to_check:
        for wellplate in full_library_mongo[collection].keys():
            if full_library_mongo[collection][wellplate]['container_name'] == wellplate_name:
                wellplate_key = [collection, wellplate]
                success_statements.append('Found the target wellplate named %s in the database' % wellplate_name)
                break
        else:
            continue
        break
    else:
        return ['Error', 'No wellplate named %s found in the database' % wellplate_name, success_statements]
    wellplate_definition = full_library_mongo[wellplate_key[0]][wellplate_key[1]]
    initial_wellplate_definition = copy.deepcopy(wellplate_definition)     
    initial_location = initial_wellplate_definition['location'][1]
       
    internal_standard_solution = additional_prep_details['internal_standard_solution']
    standard_well_volume = additional_prep_details['total_standard_well_volume']
    contents_to_prepare = queue_document['containers'][plate_operations[operation_number]['container']]['contents']
    plate_to_prepare = {'Incoming': {}, 'Plate_ID': queue_name}
    for plate_well_key in contents_to_prepare.keys():
        if 'A' in plate_well_key:
            standard_plate_well_key = plate_well_key.replace('A', 'H')
            plate_to_prepare['Incoming'][standard_plate_well_key] = {'plate_well': standard_plate_well_key,
                                                 'reagents': contents_to_prepare[plate_well_key]['reagents'],
                                                 'solvent': [internal_standard_solution],
                                                 'total_volume': standard_well_volume}
            for detail_key in contents_to_prepare[plate_well_key].keys():
                if detail_key not in plate_to_prepare['Incoming'][standard_plate_well_key].keys():
                    plate_to_prepare['Incoming'][standard_plate_well_key][detail_key] = contents_to_prepare[plate_well_key][detail_key]
    
    # TODO: Which reagents do we not bother adding into the well-plate?
    # Track down all of the reagents that we need for the reaction wellplate, these
    # are found first because the volume of the solvent added depends on the reagents
    reagents_needed = {}
    for well_location in plate_to_prepare['Incoming']:
        reagents_to_find = plate_to_prepare['Incoming'][well_location]['reagents']
        if reagents_to_find[0] == '':
            continue
        for reagent in reagents_to_find:
            if reagent in ['', ':']:
                continue
            if float(reagent[1]) == 0:
                continue
            if reagent[0] not in reagents_needed.keys():
                reagents_needed[reagent[0]] = {'destination': [[well_location, reagent[1]]]}
            else:
                reagents_needed[reagent[0]]['destination'].append([well_location, reagent[1]])
    reagent_locations = {}
    problematic_reagents = []
    for reagent in reagents_needed.keys():
        reagent_information_return = chemical_well_origin(full_library_mongo, reagent, reagents_needed[reagent])
        if reagent_information_return[0] == 'Error':
            problematic_reagents.append(reagent_information_return[1])
            continue
        elif reagent_information_return[0] == 'Success':
            reagent_information = reagent_information_return[1]
        else:
            problematic_reagents.append(reagent_information_return)
        reagent_locations[reagent] = reagent_information
        full_library_mongo['reagents'][reagent_information['plate_id']]['contents'][
            reagent_information['well_location']]['volume_ul'] -= \
            reagent_information['final_well_details']['volume_ul']
    # Now that we have the volume of each of the reagents needed we can get the solvents needed
    # and then update the library after finding the needed solvent
    solvents_needed = {}
    for well_location in plate_to_prepare['Incoming']:
        reaction_total_volume = plate_to_prepare['Incoming'][well_location]['total_volume']
        reagents_to_find = plate_to_prepare['Incoming'][well_location]['reagents']
        reagent_volume_needed = 0
        for reagent in reagents_to_find:
            try:
                destinations = reagent_locations[reagent[0]]['destination']
            except KeyError:
                continue
            for item in destinations:
                if item[0] == well_location:
                    volume_needed = item[2]
                    break
            else:
                return ['Error', 'A coding error prevented the solvent well_location and the reagent_well location connection',
                        success_statements]
            reagent_volume_needed += volume_needed
        solvent_volume_needed = reaction_total_volume - reagent_volume_needed
        solvents_to_find = plate_to_prepare['Incoming'][well_location]['solvent'][0].split()
        for solvent in solvents_to_find:
            if solvent not in solvents_needed.keys():
                solvents_needed[solvent] = {'destination': [[well_location, solvent_volume_needed]]}
            else:
                solvents_needed[solvent]['destination'].append([well_location, solvent_volume_needed])
    solvent_locations = {}
    problematic_solvents = []
    for solvent in solvents_needed.keys():
        solvent_information_return = solvent_finder(full_library_mongo, solvent, solvents_needed[solvent])
        if solvent_information_return[0] == 'Error':
            problematic_solvents.append(solvent_information_return[1])
            continue
        elif solvent_information_return[0] == 'Success':
            solvent_information = solvent_information_return[1]
        else:
            problematic_solvents.append(solvent_information_return)
        solvent_locations[solvent] = solvent_information
        full_library_mongo['solvents'][solvent_information['container_id']]['contents'][
            solvent_information['well_location']]['volume_ul'] -= \
            solvent_information['final_solvent_details']['volume_ul']

    if len(problematic_reagents) > 0 or len(problematic_solvents) > 0:
        success_statements.append(['Reagents', problematic_reagents])
        success_statements.append(['Solvents', problematic_solvents])
        return ['Error', 'Not enough material to prepare the wellplate', success_statements]
    # With the solvents and reagents we can group them by originating location
    return_statement, reagent_order = priority_sort(
        [[reagent, reagent_locations[reagent]['plate_id'], len(reagent_locations[reagent]['destination']),
          reagent_locations[reagent]['container_name']] for reagent in reagent_locations.keys()])
    if return_statement == 'Error':
        return ['Error', reagent_order, success_statements]
    return_statement, solvent_order = priority_sort(
        [[solvent, solvent_locations[solvent]['container_id'], len(solvent_locations[solvent]['destination']),
          solvent_locations[solvent]['container_name']] for solvent in solvent_locations.keys()])
    if return_statement == 'Error':
        return ['Error', solvent_order, success_statements]
    success_statements.append('Found all of the reagents and solvents on the platform')

    # This gives us all of the chemicals that we need to distribute into wells so now
    # we can build up the schedule of tasks needed to distribute the chemicals, we
    # need to get the locations that we are placing the plates at and also find an
    # empty reaction plate that we can use for the method
    simplified_schedule, method_strings = [], []
    collections_to_check = ['wellplates', 'reagents']
    for collection in collections_to_check:
        for wellplate in full_library_mongo[collection].keys():
            if full_library_mongo[collection][wellplate]['container_name'] == wellplate_name:
                wellplate_key = [collection, wellplate]
                success_statements.append('Found the target wellplate named %s in the database' % wellplate_name)
                break
        else:
            continue
        break
    else:
        return ['Error', 'No wellplate named %s found in the database' % wellplate_name, success_statements]
    wellplate_definition = full_library_mongo[wellplate_key[0]][wellplate_key[1]]
    initial_wellplate_definition = copy.deepcopy(wellplate_definition)     
    initial_location = initial_wellplate_definition['location'][1]
    # Now we need to find a suitable location for the preparation to take place at
    # TODO: Determine if this plate type is even relevant still??
    if target_well_plate_type == 'Special High T Plate':
        preparation_location = [20, 2, 323]
    else:
        allowed_prep_carriers = [32, 38, 44] #TODO Add other locations such as the heater-shakers!!
        if initial_location[0] in allowed_prep_carriers:
            preparation_location = initial_location
        else:
            return_statement, preparation_location = find_accessible_open_location(updated_carriers)
            if return_statement == 'Error':
                return ['Error', preparation_location, success_statements]
    success_statements.append('Found a suitable prep location at location: %s' % preparation_location)
    if initial_wellplate_definition['location'][1] != preparation_location:
        simplified_schedule.append(['TRANSFER_LABWARE', initial_wellplate_definition['container_name'],
                                    initial_wellplate_definition['location'][1], preparation_location,
                                    str(initial_wellplate_definition['_id']), 'wellplates'])
    else:
        success_statements.append('Target well-plate %s is already in a suitable preparation location %s' % (wellplate_name, preparation_location))
    
    # Here is an optional statement to prepare the backing lines of the LiHa for dispensing
    if additional_prep_details['tip_prep_rinse'] == True:
        simplified_schedule.append(['TIP_RINSE', 10, 10])
    else:
        simplified_schedule.append(['TIP_RINSE', 3, 3])
    simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 250])
    simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 250])
    # Now the aspirate and dispense functions are added to the simplified schedule
    solvent_index = 0
    for platform_location in solvent_order:
        solvent_definition = full_library_mongo['solvents'][platform_location[0]]
        for chemical in platform_location[1]:
            solvent_origin = solvent_definition['location']
            origin_labware_type = carriers_labware['labware'][solvent_definition['labware_type']]
            destination_labware_type = carriers_labware['labware'][initial_wellplate_definition['labware_type']]
            return_statement, liha_groups = liha_grouping(solvents_needed[chemical[0]]['destination'], origin_labware_type,
                                        destination_labware_type)
            if return_statement == 'Error':
                return ['Error', liha_groups, success_statements]
            for group_no, liha_group in enumerate(liha_groups.keys()):
                group = liha_groups[liha_group]
                wells, volumes, tips = map(list, zip(*group))
                simplified_schedule.append(['LiHa_ASPIRATE_DISPENSE', solvent_origin[1], preparation_location,
                                            tips, volumes, wells, 'Adjusted water free dispense'])
                if additional_prep_details['empty_wells'] is True and solvent_index == 0:
                    continue
                else:
                    continue
                    simplified_schedule.append(['MOVE_LiHa'])
                    simplified_schedule.append(['TIP_RINSE', 3, 3])
                    simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
                    simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
                if additional_prep_details['washing_frequency'] != 'None':
                    simplified_schedule.append(['TIP_RINSE', 3, 3])
                    simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
                    simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
            solvent_index += 1
    simplified_schedule.append(['MOVE_LiHa'])
    simplified_schedule.append(['TIP_RINSE', 5, 5])
    simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
    simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
    
    acceptable_bed_locations = [32, 38, 44]
    location_occupancies = []
    liha_tip_number = 1
    for platform_location in reagent_order:
        reagent_tray_definition = full_library_mongo['reagents'][platform_location[0]]
        current_plate_location = reagent_tray_definition['location'][1]
        success_statements.append('Currently the reagent tray "%s" is located at location: %s' % (
            reagent_tray_definition['container_name'], current_plate_location))
        if reagent_tray_definition['location'][1][0] not in acceptable_bed_locations:
            bad_locations = []
            for i in range(0, 5):
                return_statement, reagent_well_plate_home = find_accessible_open_location(updated_carriers, exceptions=[
                                                                        preparation_location] + bad_locations)
                if return_statement == 'Error':
                    return ['Error', reagent_well_plate_home, success_statements]
                bad_location = []
                for location in location_occupancies:
                    if reagent_well_plate_home == location[0]:
                        if reagent_tray_definition['labware_type'] != location[1]:
                            bad_location = [location[0]]
                        else:
                            break
                if len(bad_location) > 0:
                    bad_locations += bad_location
                else:
                    break
            if reagent_tray_definition['location'] != reagent_well_plate_home:
                simplified_schedule.append(['TRANSFER_LABWARE', reagent_tray_definition['container_name'],
                                            reagent_tray_definition['location'][1], reagent_well_plate_home,
                                            str(reagent_tray_definition['_id']), 'reagents'])
                current_plate_location = reagent_well_plate_home
        location_occupancies.append([current_plate_location, reagent_tray_definition['labware_type']])
        for chemical in platform_location[1]:
            # reagent_origin = reagent_tray_definition['location']
            reagent_origin_well = reagent_locations[chemical[0]]['initial_well_details']['plate_well']
            for destination in reagent_locations[chemical[0]]['destination']:
                tips = [liha_tip_number]
                volumes = [destination[2]]
                wells = [destination[0]]
                simplified_schedule.append(['LiHa_ASPIRATE_DISPENSE', current_plate_location + [[reagent_origin_well]],
                                            preparation_location, tips, volumes, wells, 'Adjusted water free dispense'])
            if liha_tip_number == 8:
                simplified_schedule.append(['TIP_RINSE', 5, 5])
                simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
                simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
                liha_tip_number = 1
            else:
                liha_tip_number += 1
        simplified_schedule.append(['MOVE_LiHa'])
        if reagent_tray_definition['location'][1] != current_plate_location:
            simplified_schedule.append(['TRANSFER_LABWARE', reagent_tray_definition['container_name'],
                                        current_plate_location, reagent_tray_definition['location'][1],
                                        str(reagent_tray_definition['_id']), 'reagents'])
    simplified_schedule.append(['TIP_RINSE', 5, 5])
    simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
    simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
    simplified_schedule.append(['TIP_RINSE', 3, 4])
    
    if initial_wellplate_definition['location'][1] != preparation_location:
        simplified_schedule.append(['TRANSFER_LABWARE', initial_wellplate_definition['container_name'],
                                    preparation_location, initial_wellplate_definition['location'][1],
                                    str(initial_wellplate_definition['_id']), 'wellplates'])
    
    # Now all of the used locations on the bed needs to be populated to avoid errors
    # when EVOWARE attempts to execute the method
    for task in simplified_schedule:
        if 'TRANSFER_LABWARE' in task[0]:
            updated_carriers[task[2][0]][task[2][2]]['labware_labels'][task[2][1] - 1] = 'location_%s_%s_%s' % (
                str(task[2][0]), str(task[2][1]), str(task[2][2]))
            updated_carriers[task[2][0]][task[2][2]]['labware_types'][task[2][1] - 1] = \
                full_library_mongo[task[5]][task[4]]['labware_type']
            updated_carriers[task[3][0]][task[3][2]]['labware_labels'][task[3][1] - 1] = 'location_%s_%s_%s' % (
                str(task[3][0]), str(task[3][1]), str(task[3][2]))
            updated_carriers[task[3][0]][task[3][2]]['labware_types'][task[3][1] - 1] = \
                full_library_mongo[task[5]][task[4]]['labware_type'] 
    
    # Build the schedule from the method strings
    # Build the schedule from the method strings
    return_statement, filled_worktable_export = worktable_export(carriers_labware_filepath, empty_worktable_filepath, updated_carriers)
    if return_statement == 'Error':
        return ['Error', filled_worktable_export, success_statements]
    for item in simplified_schedule:
        if item[0] == 'TRANSFER_LABWARE':
            return_statement, method_string = transfer_labware(updated_carriers, full_library_mongo[item[5]], item)
            method_strings.append(method_string)
            updated_carriers[item[3][0]][item[3][2]]['labware_labels'][item[3][1] - 1] = copy.deepcopy(
                updated_carriers[item[2][0]][item[2][2]]['labware_labels'][item[2][1] - 1])
            updated_carriers[item[2][0]][item[2][2]]['labware_labels'][item[2][1] - 1] = ''
            updated_carriers[item[3][0]][item[3][2]]['labware_types'][item[3][1] - 1] = copy.deepcopy(
                updated_carriers[item[2][0]][item[2][2]]['labware_types'][item[2][1] - 1])
            updated_carriers[item[2][0]][item[2][2]]['labware_types'][item[2][1] - 1] = ''
        elif item[0] == 'SPECIAL_TRANSFER_LABWARE':
            return_statement, method_string = special_transfer_labware(updated_carriers, item[1], item)
            method_strings.append(method_string)
        elif item[0] == 'SPECIAL_LID_TRANSFER':
            return_statement, method_string = special_lid_transfer_labware(updated_carriers, full_library_mongo[item[7]], item)
            method_strings.append(method_string)
        elif item[0] == 'MCA_GET_TIPS':
            return_statement, method_string = mca_get_tips(item[1])
            method_strings.append(method_string)
        elif item[0] == 'MCA_MIX':
            return_statement, method_string = mca_mix(*item[1::])
            method_strings.append(method_string)
        elif item[0] == 'TEVAC_FUNCTION':
            return_statement, method_string = tevac_string_return(*item[1::])
            method_strings.append(method_string)
            #if item[1][0] == 'vacuum_on':
            #    method_strings.append('ROMA(0,55,75,0,0,0,150,1,0);')
            #    method_strings.append('Vector("TeVacS Custom_Narrow_1","12","10",0,1,2,2,0,0);')
        elif item[0] == 'MCA_ASPIRATE':
            return_statement, method_string = mca_aspirate(*item[1::])
            method_strings.append(method_string)
        elif item[0] == 'MCA_DISPENSE':
            return_statement, method_string = mca_dispense(*item[1::])
            method_strings.append(method_string)
        elif item[0] == 'MCA_DROP_TIPS':
            return_statement, method_string = mca_droptips(updated_carriers)
            method_strings.append(method_string)
        elif item[0] == 'MCA_MOVE':
            return_statement, method_string = mca_move()
            method_strings.append(method_string)
        elif item[0] == 'MCA_SOFT_MIX':
            method_strings.append('BeginLoop("%s","soft_mixing");' % item[5])
            if item[4] == 'top':
                return_statement, method_string = mca_aspirate(item[1], item[3], 'Adjusted Water free dispense DiTi 200 Clear Fast Top')
                method_strings.append(method_string)
                return_statement, method_string = mca_dispense(item[2], item[3], 'Adjusted Water free dispense DiTi 200 Clear Fast Bottom')
                method_strings.append(method_string) 
            elif item[4] == 'bottom':
                return_statement, method_string = mca_aspirate(item[1], item[3], 'Adjusted Water free dispense DiTi 200 Clear Fast Bottom')
                method_strings.append(method_string)
                return_statement, method_string = mca_dispense(item[2], item[3], 'Adjusted Water free dispense DiTi 200 Clear Fast Top')
                method_strings.append(method_string)
            method_strings.append('EndLoop();')
        elif item[0] == 'MCA_LAYER':
            if item[4] == 'top':
                return_statement, method_string = mca_aspirate(item[1], item[3], 'Adjusted Water free dispense DiTi 200 Clear Fast Top')
                method_strings.append(method_string)
                return_statement, method_string = mca_dispense(item[2], item[3], 'Adjusted Water free dispense DiTi 200 Clear Fast Bottom')
                method_strings.append(method_string) 
            elif item[4] == 'bottom':
                return_statement, method_string = mca_aspirate(item[1], item[3], 'Adjusted Water free dispense DiTi 200 Clear Fast Bottom')
                method_strings.append(method_string)
                return_statement, method_string = mca_dispense(item[2], item[3], 'Adjusted Water free dispense DiTi 200 Clear Fast Top')
                method_strings.append(method_string)
        elif item[0] == 'START_TIMER':
            return_statement, method_string = start_timer(item[1])
            method_strings.append(method_string)
        elif item[0] == 'WAIT_FOR_TIMER':
            return_statement, method_string = wait_for_timer(*item[1::])
            method_strings.append(method_string)
        elif item[0] == 'LARGE_WASTE_DUMP':
            return_statement, method_string = large_waste_dump(*item[1::], updated_carriers)
            method_strings.append(method_string)
        elif 'LiHa_ASPIRATE_DISPENSE' in item[0]:
            return_statement, aspirate_val, dispense_val = liha_aspirate_dispense(updated_carriers, carriers_labware, item[1:-1], item[-1])
            if return_statement == 'Error':
                return ['Error', aspirate_val, success_statements]
            method_strings.append(aspirate_val)
            method_strings.append(dispense_val)
            continue
        elif 'TIP_RINSE' in item[0]:
            return_statement, method_string = wash_tips()
            method_strings.append(method_string)
        elif 'MOVE_LiHa' in item[0]:
            return_statement, method_string = move_liha()
            method_strings.append(method_string)
        elif 'SOLVENT_RINSE' in item[0]:
            return_statement, solvent_rinse_origin = solvent_rinse_location(item[1], full_library_mongo)
            if return_statement == 'Error':
                return ['Error', solvent_rinse_origin, success_statements]
            return_statement, aspirate_val, dispense_val = liha_aspirate_dispense(updated_carriers, carriers_labware,
                                                                [solvent_rinse_origin, [1, 2, 30], [1, 2, 3, 4, 5, 6, 7, 8],
                                                                [int(item[2])] * 8, ['A1', 'B1', 'C1', 'D1', 'E1', 'F1', 'G1', 'H1']])
            if return_statement == 'Error':
                return ['Error', aspirate_val, success_statements]
            method_strings.append(aspirate_val)
            method_strings.append(dispense_val)
            continue
        else:
            return ['Error', 'LiHa Operation %s not defined in the method file' % item, method_strings]
        if return_statement == 'Error':
            return ['Error', method_string, success_statements]
    with open(schedule_path_output, 'w', encoding='latin1') as outfile:
        for line in filled_worktable_export:
            outfile.write(line + '\n')
        for line in method_strings:
            outfile.write(line + '\n')
    
    # Now that the schedule was generated, the database needs to be updated
    # TODO: Include an output to send the current library for roll-backs to use
    # when the liquid handler encounters an error
    
    # TODO: How to break-down the internal standard into components to add into the well-plate
    # internal_standard_solution
    
    for internal_standard in solvent_locations:
        internal_standard_definition = solvent_locations[solvent]['final_solvent_details']
        internal_standard_destinations = solvent_locations[solvent]['destination']
    internal_standard_name = internal_standard_definition['internal_standard']
    internal_standard_concentration = internal_standard_definition['internal_standard_conc']
    final_wellplate_definition = copy.deepcopy(wellplate_definition)
    for destination in internal_standard_destinations:
        plate_well_key = destination[0]
        internal_standard_volume = destination[1]
        prepared_well = plate_to_prepare['Incoming'][plate_well_key]
        
        reagent_list = []
        for reagent in prepared_well['reagents']:
            reagent_information = reagent_locations[reagent[0]]
            reagent_list.append([reagent_information['final_well_details']['chemical_name'], 
                                 reagent[1], reagent_information['final_well_details']['chemical_smiles']])
        prepared_well['reagents'] = reagent_list
        solvent_list = []
        for solvent in prepared_well['solvent']:
            solvent_information = solvent_locations[solvent]
            solvent_list.append([solvent_information['final_solvent_details']['chemical_name'],
                                solvent_information['final_solvent_details']['chemical_smiles']])
        prepared_well['solvent'] = solvent_list

        prepared_well['internal_standard'] = internal_standard_name
        prepared_well['internal_standard_concentration'] = internal_standard_concentration
        prepared_well['internal_standard_volume'] = internal_standard_volume
        final_wellplate_definition['contents'][plate_well_key] = prepared_well
    final_wellplate_definition['location'] = initial_wellplate_definition['location']
    full_library_mongo['wellplates'][str(wellplate_definition['_id'])] = final_wellplate_definition

    platform_mongodb = pymongo.MongoClient(platform_mongodb_address, platform_mongodb_port)
    for library_key in full_library_mongo.keys():
        previous = initial_full_library[library_key]
        updated = full_library_mongo[library_key]
        return_statements = update_database(previous, updated, library_key, platform_mongodb)
    updated_queue_document = copy.deepcopy(queue_document)
    if __name__ == "__main__":
        updated_queue_document['operations'][operation_number]['completed'] = 'yes'
    flag_return, return_statements = update_database_document(queue_document, updated_queue_document, 'queue',
                                                             platform_mongodb)
    platform_mongodb.close()
    success_statements.append('Updated the platform library to be up-to-date following this run')
    return ['Success', schedule_path_output, success_statements]
 
if __name__ == "__main__":
    queuename = 'tztz_test_plate'
    prep_return = prepare_standard_wells(queuename, '5')
    pprint(prep_return)