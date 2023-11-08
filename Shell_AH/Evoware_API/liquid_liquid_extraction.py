# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 15:25:47 2021

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
from mergedeep import merge

if __name__ == "__main__":
    from platform_library_functions import query_collection, update_database, query_document, \
        update_database_document, move_database_document
    from tecan_functions import mca_move, find_mca_tips, tevac_string_return, start_timer, wait_for_timer, \
        mca_droptips, mca_get_tips, mca_mix, mca_aspirate, mca_dispense
    from tecan_functions import transfer_labware, wash_tips, liha_aspirate_dispense, \
        move_liha, special_transfer_labware, special_lid_transfer_labware, large_waste_dump
    from useful_other_functions import solvent_rinse_location, solvent_finder, \
        liha_grouping, find_accessible_open_location, directory_check, find_labware
    from worktable_cleaner import initial_worktable_prep, worktable_export
else:
    from Evoware_API.platform_library_functions import query_collection, update_database, query_document, \
        update_database_document, move_database_document
    from Evoware_API.tecan_functions import mca_move, find_mca_tips, tevac_string_return, start_timer, wait_for_timer, \
        mca_droptips, mca_get_tips, mca_mix, mca_aspirate, mca_dispense
    from Evoware_API.tecan_functions import transfer_labware, wash_tips, liha_aspirate_dispense, \
        move_liha, special_transfer_labware, special_lid_transfer_labware, large_waste_dump
    from Evoware_API.useful_other_functions import solvent_rinse_location, solvent_finder, \
        liha_grouping, find_accessible_open_location, directory_check, find_labware
    from Evoware_API.worktable_cleaner import initial_worktable_prep, worktable_export


def liquid_liquid_extraction(queue_name, operation_number):
    try:
        # ----------------------------------------------------------------------------
        # ----------------------------------------------------------------------------
        # Open the configuration file before starting anything
        success_statements = []
        config_filepath = r'C:\Users\admin\PycharmProjects\AMD_Control_Platform\venv\Shell_AH\Evoware_API\config.yaml'
        with open(config_filepath, 'r') as yamlfile:
            config = yaml.load(yamlfile, Loader=yaml.Loader)
        
        # Grab the important paths from the configuration file
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
    
        # Start by loading the specific queue element name and check whether it is real
        platform_mongodb = pymongo.MongoClient(platform_mongodb_address, platform_mongodb_port)
        queue_document, queue_return_statements = query_document('queue', queue_name, 'queue_name', platform_mongodb)
        if type(queue_document) == str:
            platform_mongodb.close()
            return ['Error', queue_return_statements, success_statements]
        else:
            success_statements.append('Found element %s in the queue collection' % queue_name)
    
        # Now we check to make sure that the plate is ready to be prepared
        plate_operations = queue_document['operations']
        queue_containers = queue_document['containers']
        target_operation = 'liquid_liquid_extraction'
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
        target_well_plate_type = queue_containers[target_container]['plate_type']
        success_statements.append('The %s operation is in a valid position and ready to be executed' % target_operation)
        # TODO check to make sure the required details provided are there or provide default values (??)
    
        # Here there is an optional statement from config to send to the testing folder
        current_time = datetime.datetime.now()
        if config['startup']['testing']:
            schedule_path_output = config['paths/files']['schedule_testing_output_path'] % current_time.strftime('%Y%m%d')
        else:
            folder_outpath = config['paths/files']['schedule_output_path'] % current_time.strftime('%Y%m%d')
            directory_check(folder_outpath)
            schedule_path_output = folder_outpath + r'\%s_%s_%s_%s.esc' % (
                queue_name, operation_number, 'liquid_extraction', current_time.strftime('%Y%m%d%H%M%S'))
        
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
        with open(carriers_labware_filepath, 'r') as jsonfile:
            carriers_labware = json.load(jsonfile)
        return_statement, start_carriers = initial_worktable_prep(carriers_labware_filepath, empty_worktable_filepath, full_library_mongo)
        if return_statement != 'Success':
            return ['Error', start_carriers, success_statements]
        updated_carriers = copy.deepcopy(start_carriers)
        success_statements.append('Loaded the carrier and labware library file!')
        
        # Since a wellplate will be made during the method, get both plate types
        initial_plate_details = {'plate_type': target_well_plate_type,
                                 'container_name': queue_containers[plate_operations[operation_number]['container']][
                                     'container_name'],
                                 'queue_container_name': plate_operations[operation_number]['container']}
        final_plate_details = {'plate_type': queue_containers[additional_prep_details['target_container']]['plate_type'],
                               'container_name': queue_containers[additional_prep_details['target_container']][
                                   'container_name'],
                               'queue_container_name': additional_prep_details['target_container']}
        if initial_plate_details['container_name'] is None:
            return ['Error', 'No wellplate name specified for "%s"' % plate_operations[operation_number]['container'],
                    success_statements]
        
        # Now find the initial well-plate on the platform
        return_statement, extra_initial_plate_details = find_labware(updated_carriers, full_library_mongo, 
                                                                     initial_plate_details['container_name'],
                                                                     initial_plate_details['plate_type'],
                                                                     initial_plate_details['queue_container_name'])
        if return_statement != 'Success':
            return ['Error', extra_initial_plate_details, success_statements]
        merge(initial_plate_details, extra_initial_plate_details)
        initial_plate_details['wellplate_key'] = ['wellplates', str(initial_plate_details['_id'])]
        success_statements.append('Found the initial wellplate %s at %s' % (initial_plate_details['container_name'], 
                                                                            initial_plate_details['location'][1]))
        
        # Check to make sure the the initial plate is in an accessible location
        simplified_schedule, method_strings = [], []
        allowed_prep_carriers = [32, 38, 44]
        if initial_plate_details['location'][1][0] in allowed_prep_carriers:
            initial_plate_details['prep_location'] = initial_plate_details['location'][1]
        else:
            return_statement, open_location = find_accessible_open_location(updated_carriers)
            if return_statement == 'Error':
                return ['Error', open_location, success_statements]
            initial_plate_details['prep_location'] = open_location
        success_statements.append('Found a suitable prep location at location: %s' % initial_plate_details['prep_location'])
        if initial_plate_details['location'][1] != initial_plate_details['prep_location']:
            simplified_schedule.append(['TRANSFER_LABWARE', initial_plate_details['container_name'],
                                        initial_plate_details['location'][1], initial_plate_details['prep_location'],
                                        str(initial_plate_details['_id']), 'wellplates'])
        success_statements.append('Found a suitable prep location plate %s at %s' % (
            initial_plate_details['container_name'], initial_plate_details['prep_location']))
    
        # Now find the final well-plate to use on the platform      
        return_statement, extra_final_plate_details = find_labware(updated_carriers, full_library_mongo, 
                                                                       final_plate_details['container_name'], 
                                                                       final_plate_details['plate_type'],
                                                                   final_plate_details['queue_container_name'])
        if return_statement != 'Success':
            return ['Error', extra_final_plate_details, success_statements]
        merge(final_plate_details, extra_final_plate_details)
        final_plate_details['wellplate_key'] = ['wellplates', str(final_plate_details['_id'])]
        success_statements.append('Found the filtrate plate %s at %s' % (final_plate_details['container_name'], 
                                                                         final_plate_details['location'][1]))
        
        # Check to make sure the the final plate is in an accessible location
        allowed_prep_carriers = [32, 38, 44]
        if final_plate_details['location'][1][0] in allowed_prep_carriers:
            final_plate_details['prep_location'] = final_plate_details['location'][1]
        else:
            return_statement, open_location = find_accessible_open_location(updated_carriers,
                                                                                 exceptions = [initial_plate_details['prep_location']])
            if return_statement == 'Error':
                return ['Error', open_location, success_statements]
            final_plate_details['prep_location'] = open_location
        success_statements.append('Found a suitable prep location at location: %s' % final_plate_details['prep_location'])
        if final_plate_details['location'][1] != final_plate_details['prep_location']:
            simplified_schedule.append(['TRANSFER_LABWARE', final_plate_details['container_name'],
                                        final_plate_details['location'][1], final_plate_details['prep_location'],
                                        str(final_plate_details['_id']), 'wellplates'])
        success_statements.append('Found a suitable prep location plate %s at %s' % (
            final_plate_details['container_name'], final_plate_details['prep_location']))
        
        if int(additional_prep_details['extraction_iterations']) == 0:
            return ['Error', 'No extraction iterations were specified, nothing to do?', success_statements]
        
        if float(additional_prep_details['other_phase_volume_to_add']) > 0:
            solvent_to_find = additional_prep_details['other_phase_solvent']
            if 'internal_standards' not in additional_prep_details.keys():
                additional_prep_details['internal_standards'] = False
            if additional_prep_details['internal_standards'] == True:
                relevant_wells = [plate_well_key for plate_well_key in initial_plate_details['contents'].keys() 
                                  if 'internal_standard' not in initial_plate_details['contents'][plate_well_key]]
            else:
                relevant_wells = list(initial_plate_details['contents'].keys())
            
            volume_to_add = float(additional_prep_details['other_phase_volume_to_add'])
            
            solvents = {solvent_to_find: {'destination': []}}
            for relevant_well in relevant_wells:
                well_details = initial_plate_details['contents'][relevant_well]
                existing_water_volume = 0
                dmf_dmso_solvent_flag = 0
                for solvent in well_details['solvents']:
                    if solvent[0] in ['Water', 'water']:
                        existing_water_volume += solvent[1]
                    if solvent[0] in ['dmf', 'dmso', 'DMF', 'DMSO']:
                        dmf_dmso_solvent_flag = 1
                if dmf_dmso_solvent_flag == 1:
                    continue
                remaining_volume_to_add = volume_to_add - existing_water_volume
                if remaining_volume_to_add <= 0:
                    continue
                solvents[solvent_to_find]['destination'].append([relevant_well, volume_to_add, '', 1])
            
            #pprint(solvents)
            #solvents = {solvent_to_find: {'destination': [[relevant_well, volume_to_add]
            #                                              for relevant_well in relevant_wells]}}
            solvent_information_return = solvent_finder(full_library_mongo, solvent_to_find, solvents[solvent_to_find])
            if return_statement == 'Error':
                return ['Error', solvent_information_return, success_statements]
            if solvent_information_return[0] == 'Error':
                success_statements.append(['Solvents', solvent_information_return[1]])
                return ['Error', 'Issue finding solvents on the liquid handler', success_statements]
            elif solvent_information_return[0] == 'Success':
                solvent_information = solvent_information_return[1]
            else:
                success_statements.append(['Solvents', solvent_information_return])
                return ['Error', 'Unknown issue finding solvents on the liquid handler', success_statements]
            full_library_mongo['solvents'][solvent_information['container_id']]['contents'][
                solvent_information['well_location']]['volume_ul'] -= \
                solvent_information['final_solvent_details']['volume_ul']
            solvent_information['library_details'] = full_library_mongo['solvents'][solvent_information['container_id']]
            solvent_origin = solvent_information['library_details']['location']
            origin_labware_type = carriers_labware['labware'][solvent_information['library_details']['labware_type']]
            destination_labware_type = carriers_labware['labware'][initial_plate_details['plate_type']]
            simplified_schedule.append(['TIP_RINSE', 10, 10])
            simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
            simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
            return_statement, liha_groups = liha_grouping(solvents[solvent_to_find]['destination'], origin_labware_type,
                                        destination_labware_type)
            if return_statement == 'Error':
                return ['Error', liha_groups, success_statements]
            # This is only true because the TeVaC dispense liquid class dispenses
            # does a free dispense at the moment, this can be changed as needed
            
            # Now we need to check to see if we need to handle an acid-base quenching
            acids_list = ['sulfuric_acid', 'nitric_acid']
            bases_list = ['water_naoh']
            for group_no, liha_group in enumerate(liha_groups.keys()):
                group = liha_groups[liha_group]
                wells, volumes, tips = map(list, zip(*group))
                acids_check = 0
                bases_check = 0
                if solvent_to_find in acids_list:
                    acids_check = 1
                if solvent_to_find in bases_list:
                    bases_check = 1
                for well in wells:
                    well_details = initial_plate_details['contents'][well]
                    for solvent in well_details['solvents']:
                        if solvent[0] in acids_list:
                            acids_check = 1
                        if solvent[0] in bases_list:
                            bases_check = 1
                if acids_check == 1 and bases_check == 1:
                    liquid_class = 'Slow Quench Dispense'
                else:
                    liquid_class = 'Adjusted water free dispense Filter Prep'
                
                simplified_schedule.append( ['LiHa_ASPIRATE_DISPENSE', solvent_origin[1], initial_plate_details['prep_location'],
                                             tips, volumes, wells, liquid_class])
            simplified_schedule.append(['MOVE_LiHa'])
            simplified_schedule.append(['TIP_RINSE', 3, 3])
            simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
            simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
        
        # Now that the other phase solvent has been added, we are ready for the extractions
        for extraction_iteration in range(0, additional_prep_details['extraction_iterations']):
            if extraction_iteration == 0:
                volume_to_add = additional_prep_details['iteration1_volume']
            else:
                volume_to_add = additional_prep_details['filter_volume_per_iteration']
            solvent_to_find = additional_prep_details['liha_solvent_to_dispense']
            if additional_prep_details['internal_standards'] == True:
                relevant_wells = [plate_well_key for plate_well_key in initial_plate_details['contents'].keys() 
                                  if 'internal_standard' not in initial_plate_details['contents'][plate_well_key]]
            else:
                relevant_wells = list(initial_plate_details['contents'].keys())
            solvents = {solvent_to_find: {'destination': []}}
            #print(relevant_wells)
            for relevant_well in relevant_wells:
                well_details = initial_plate_details['contents'][relevant_well]
                evaporated_volume = 0
                dmf_dmso_volume = 0
                for solvent in well_details['solvents']:
                    if solvent[0].lower() not in ['dmf', 'dmso', 'water']:
                        evaporated_volume += solvent[1]
                    if solvent[0] in ['dmf', 'dmso']:
                        dmf_dmso_volume += solvent[1]
                if dmf_dmso_volume > 0:
                    needed_volume = 100 + additional_prep_details['filter_volume_per_iteration'] - dmf_dmso_volume
                    if needed_volume <= 0:
                        continue
                    solvents[solvent_to_find]['destination'].append([relevant_well, needed_volume, '', 1])
                else:
                    needed_volume = volume_to_add #- (well_details['total_volume'] - evaporated_volume)
                    #print(well_details['total_volume'], evaporated_volume, dmf_dmso_volume, needed_volume)
                    if needed_volume <= 0:
                        continue
                    solvents[solvent_to_find]['destination'].append([relevant_well, needed_volume, '', 1])
            #pprint(solvents)
            #solvents = {solvent_to_find: {'destination': [[relevant_well, volume_to_add]
            #                                              for relevant_well in relevant_wells]}}
            solvent_information_return = solvent_finder(full_library_mongo, solvent_to_find, solvents[solvent_to_find])
            if solvent_information_return[0] == 'Error':
                success_statements.append(['Solvents', solvent_information_return[1]])
                return ['Error', 'Issue finding solvents on the liquid handler', success_statements]
            elif solvent_information_return[0] == 'Success':
                solvent_information = solvent_information_return[1]
            else:
                success_statements.append(['Solvents', solvent_information_return])
                return ['Error', 'Unknown issue finding solvents on the liquid handler', success_statements]
            full_library_mongo['solvents'][solvent_information['container_id']]['contents'][
                solvent_information['well_location']]['volume_ul'] -= \
                solvent_information['final_solvent_details']['volume_ul']
            solvent_information['library_details'] = full_library_mongo['solvents'][solvent_information['container_id']]
            solvent_origin = solvent_information['library_details']['location']
            origin_labware_type = carriers_labware['labware'][solvent_information['library_details']['labware_type']]
            destination_labware_type = carriers_labware['labware'][initial_plate_details['plate_type']]    
            if extraction_iteration == 0 and additional_prep_details['other_phase_volume_to_add'] <= 0:
                simplified_schedule.append(['TIP_RINSE', 10, 10])
                simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
                simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
            return_statement, liha_groups = liha_grouping(solvents[solvent_to_find]['destination'], origin_labware_type,
                                        destination_labware_type)
            if return_statement == 'Error':
                return ['Error', liha_groups, success_statements]
            
            # This is only true because the TeVaC dispense liquid class dispenses
            # does a free dispense at the moment, this can be changed as needed
            for group_no, liha_group in enumerate(liha_groups.keys()):
                group = liha_groups[liha_group]
                wells, volumes, tips = map(list, zip(*group))
                simplified_schedule.append(
                    ['LiHa_ASPIRATE_DISPENSE', solvent_origin[1], initial_plate_details['prep_location'], tips, volumes, wells,
                     'Adjusted water free dispense Filter Prep'])
            simplified_schedule.append(['MOVE_LiHa'])
            
            if extraction_iteration == 0:
                # Find suitable MCA tips on the platform
                return_statement, suitable_mca_tips, consumable_id = find_mca_tips(updated_carriers, full_library_mongo['consumables'])
                if return_statement == 'Error':
                    return ['Error', suitable_mca_tips, success_statements]
                simplified_schedule.append(['MCA_GET_TIPS', suitable_mca_tips])
                simplified_schedule.append(['MCA_MOVE'])
                simplified_schedule.append(['LARGE_WASTE_DUMP', suitable_mca_tips])
                full_library_mongo['consumables'][consumable_id]['tip_locations'].remove(str(suitable_mca_tips[1]))
            
            if additional_prep_details['number_of_soft_mixes'] != 0:
                # TODO figure out how to do the soft mixing with the liquid handler
                if additional_prep_details['important_layer'] == 'bottom':
                    simplified_schedule.append(['MCA_SOFT_MIX', initial_plate_details['prep_location'], 
                                                initial_plate_details['prep_location'], 150, 'bottom', 
                                                additional_prep_details['number_of_soft_mixes']])
                elif additional_prep_details['important_layer'] == 'top':
                    simplified_schedule.append(['MCA_SOFT_MIX', initial_plate_details['prep_location'], 
                                                initial_plate_details['prep_location'], 150, 'top',
                                                additional_prep_details['number_of_soft_mixes']])
                else:
                    return ['Error', 'Random important layer...', success_statements]
                pass
            simplified_schedule.append(['MCA_MOVE'])
            simplified_schedule.append(['START_TIMER', 1])
            simplified_schedule.append(['WAIT_FOR_TIMER', 1, additional_prep_details['settling_time']])
            if additional_prep_details['important_layer'] == 'bottom':
                simplified_schedule.append(['MCA_LAYER', initial_plate_details['prep_location'], 
                                final_plate_details['prep_location'], additional_prep_details['filter_volume_per_iteration'],
                                'bottom'])
            elif additional_prep_details['important_layer'] == 'top':
                simplified_schedule.append(['MCA_LAYER', initial_plate_details['prep_location'], 
                                final_plate_details['prep_location'], additional_prep_details['filter_volume_per_iteration'],
                                'top'])
            simplified_schedule.append(['MCA_MOVE'])
            
        
        simplified_schedule.append(['MCA_DROP_TIPS'])
        success_statements.append('Completed operations that use the MCA')
        simplified_schedule.append(['MOVE_LiHa'])
        simplified_schedule.append(['TIP_RINSE', 5, 5])
        simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
        simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
        simplified_schedule.append(['TIP_RINSE', 3, 4])
        success_statements.append('Completed operations that use the LiHa')
        
        if initial_plate_details['location'][1] != initial_plate_details['prep_location']:
            simplified_schedule.append(['TRANSFER_LABWARE', initial_plate_details['container_name'],
                                        initial_plate_details['prep_location'], initial_plate_details['location'][1],
                                        str(initial_plate_details['_id']), 'wellplates'])
        if final_plate_details['location'][1] != final_plate_details['prep_location']:
            simplified_schedule.append(['TRANSFER_LABWARE', final_plate_details['container_name'],
                                        final_plate_details['prep_location'], final_plate_details['location'][1],
                                        str(final_plate_details['_id']), 'wellplates'])
        
        # TODO Update all of the carriers with the sites that were used and the wellplate type that was present
        for task in simplified_schedule:
            if task[0] in ['TRANSFER_LABWARE', 'SPECIAL_LID_TRANSFER']:
                updated_carriers[task[2][0]][task[2][2]]['labware_labels'][task[2][1] - 1] = 'location_%s_%s_%s' % (
                    str(task[2][0]), str(task[2][1]), str(task[2][2]))
                updated_carriers[task[2][0]][task[2][2]]['labware_types'][task[2][1] - 1] = \
                    full_library_mongo[task[-1]][str(task[-2])]['labware_type']
                updated_carriers[task[3][0]][task[3][2]]['labware_labels'][task[3][1] - 1] = 'location_%s_%s_%s' % (
                    str(task[3][0]), str(task[3][1]), str(task[3][2]))
                updated_carriers[task[3][0]][task[3][2]]['labware_types'][task[3][1] - 1] = \
                    full_library_mongo[task[-1]][str(task[-2])]['labware_type']
        # Build the schedule from the method strings
        return_statement, filled_worktable_export = worktable_export(carriers_labware_filepath, empty_worktable_filepath, updated_carriers)
        if return_statement == 'Error':
            return ['Error', filled_worktable_export, success_statements]
        for item in simplified_schedule:
            if item[0] == 'TRANSFER_LABWARE':
                return_statement, method_string = transfer_labware(updated_carriers, full_library_mongo[item[5]], item)
                if return_statement == 'Error':
                    return ['Error', method_string, success_statements]
                method_strings.append(method_string)
                updated_carriers[item[3][0]][item[3][2]]['labware_labels'][item[3][1] - 1] = copy.deepcopy(
                    updated_carriers[item[2][0]][item[2][2]]['labware_labels'][item[2][1] - 1])
                updated_carriers[item[2][0]][item[2][2]]['labware_labels'][item[2][1] - 1] = ''
                updated_carriers[item[3][0]][item[3][2]]['labware_types'][item[3][1] - 1] = copy.deepcopy(
                    updated_carriers[item[2][0]][item[2][2]]['labware_types'][item[2][1] - 1])
                updated_carriers[item[2][0]][item[2][2]]['labware_types'][item[2][1] - 1] = ''
            elif item[0] == 'SPECIAL_TRANSFER_LABWARE':
                return_statement, method_string = special_transfer_labware(updated_carriers, item[1], item)
                if return_statement == 'Error':
                    return ['Error', method_string, success_statements]
                method_strings.append(method_string)
            elif item[0] == 'SPECIAL_LID_TRANSFER':
                return_statement, method_string = special_lid_transfer_labware(updated_carriers, full_library_mongo[item[7]], item)
                if return_statement == 'Error':
                    return ['Error', method_string, success_statements]
                method_strings.append(method_string)
            elif item[0] == 'MCA_GET_TIPS':
                return_statement, method_string = mca_get_tips(item[1])
                if return_statement == 'Error':
                    return ['Error', method_string, success_statements]
                method_strings.append(method_string)
            elif item[0] == 'MCA_MIX':
                return_statement, method_string = mca_mix(*item[1::])
                if return_statement == 'Error':
                    return ['Error', method_string, success_statements]
                method_strings.append(method_string)
            elif item[0] == 'TEVAC_FUNCTION':
                return_statement, method_string = tevac_string_return(*item[1::])
                if return_statement == 'Error':
                    return ['Error', method_string, success_statements]
                method_strings.append(method_string)
                #if item[1][0] == 'vacuum_on':
                #    method_strings.append('ROMA(0,55,75,0,0,0,150,1,0);')
                #    method_strings.append('Vector("TeVacS Custom_Narrow_1","12","10",0,1,2,2,0,0);')
            elif item[0] == 'MCA_ASPIRATE':
                return_statement, method_string = mca_aspirate(*item[1::])
                if return_statement == 'Error':
                    return ['Error', method_string, success_statements]
                method_strings.append(method_string)
            elif item[0] == 'MCA_DISPENSE':
                return_statement, method_string = mca_dispense(*item[1::])
                if return_statement == 'Error':
                    return ['Error', method_string, success_statements]
                method_strings.append(method_string)
            elif item[0] == 'MCA_DROP_TIPS':
                return_statement, method_string = mca_droptips(updated_carriers)
                if return_statement == 'Error':
                    return ['Error', method_string, success_statements]
                method_strings.append(method_string)
            elif item[0] == 'MCA_MOVE':
                return_statement, method_string = mca_move()
                if return_statement == 'Error':
                    return ['Error', method_string, success_statements]
                method_strings.append(method_string)
            elif item[0] == 'MCA_SOFT_MIX':
                method_strings.append('BeginLoop("%s","soft_mixing");' % item[5])
                if item[4] == 'top':
                    if initial_plate_details['plate_type'] == 'Paradox Thermal Plate':
                        liquid_class = 'Adjusted Water free dispense DiTi Paradox %s'
                    elif initial_plate_details['plate_type'] == '96 Well DeepWell':
                        liquid_class = 'Adjusted Water free dispense Diti DeepWell %s'
                    else:
                        return ['Error', 'Input wellplate not acceptable for LLE' % initial_plate_details['plate_type'],
                                success_statements]
                    return_statement, method_string = mca_aspirate(item[1], item[3], liquid_class % 'Top')
                    if return_statement == 'Error':
                        return ['Error', method_string, success_statements]
                    method_strings.append(method_string)
                    return_statement, method_string = mca_dispense(item[2], item[3], liquid_class % 'Bottom')
                    if return_statement == 'Error':
                        return ['Error', method_string, success_statements]
                    method_strings.append(method_string) 
                elif item[4] == 'bottom':
                    if initial_plate_details['plate_type'] == 'Paradox Thermal Plate':
                        liquid_class = 'Adjusted Water free dispense DiTi Paradox %s'
                    elif initial_plate_details['plate_type'] == '96 Well DeepWell':
                        liquid_class = 'Adjusted Water free dispense Diti DeepWell %s'
                    else:
                        return ['Error', 'Input wellplate not acceptable for LLE' % initial_plate_details['plate_type'],
                                success_statements]
                    return_statement, method_string = mca_aspirate(item[1], item[3], liquid_class % 'Bottom')
                    if return_statement == 'Error':
                        return ['Error', method_string, success_statements]
                    method_strings.append(method_string)
                    return_statement, method_string = mca_dispense(item[2], item[3], liquid_class % 'Top')
                    if return_statement == 'Error':
                        return ['Error', method_string, success_statements]
                    method_strings.append(method_string)
                method_strings.append('EndLoop();')
            elif item[0] == 'MCA_LAYER':
                if item[4] == 'top':
                    if initial_plate_details['plate_type'] == 'Paradox Thermal Plate':
                        liquid_class = 'Adjusted Water free dispense DiTi Paradox %s'
                    elif initial_plate_details['plate_type'] == '96 Well DeepWell':
                        liquid_class = 'Adjusted Water free dispense Diti DeepWell %s'
                    else:
                        return ['Error', 'Input wellplate not acceptable for LLE' % initial_plate_details['plate_type'],
                                success_statements]
                    return_statement, method_string = mca_aspirate(item[1], item[3], liquid_class % 'Top')
                    if return_statement == 'Error':
                        return ['Error', method_string, success_statements]
                    method_strings.append(method_string)
                    return_statement, method_string = mca_dispense(item[2], item[3], liquid_class % 'Bottom')
                    if return_statement == 'Error':
                        return ['Error', method_string, success_statements]
                    method_strings.append(method_string) 
                elif item[4] == 'bottom':
                    if initial_plate_details['plate_type'] == 'Paradox Thermal Plate':
                        liquid_class = 'Adjusted Water free dispense DiTi Paradox %s'
                    elif initial_plate_details['plate_type'] == '96 Well DeepWell':
                        liquid_class = 'Adjusted Water free dispense Diti DeepWell %s'
                    else:
                        return ['Error', 'Input wellplate not acceptable for LLE' % initial_plate_details['plate_type'],
                                success_statements]
                    return_statement, method_string = mca_aspirate(item[1], item[3], liquid_class % 'Bottom')
                    if return_statement == 'Error':
                        return ['Error', method_string, success_statements]
                    method_strings.append(method_string)
                    return_statement, method_string = mca_dispense(item[2], item[3], liquid_class % 'Top')
                    if return_statement == 'Error':
                        return ['Error', method_string, success_statements]
                    method_strings.append(method_string)
            elif item[0] == 'START_TIMER':
                return_statement, method_string = start_timer(item[1])
                if return_statement == 'Error':
                    return ['Error', method_string, success_statements]
                method_strings.append(method_string)
            elif item[0] == 'WAIT_FOR_TIMER':
                return_statement, method_string = wait_for_timer(*item[1::])
                if return_statement == 'Error':
                    return ['Error', method_string, success_statements]
                method_strings.append(method_string)
            elif item[0] == 'LARGE_WASTE_DUMP':
                return_statement, method_string = large_waste_dump(*item[1::], updated_carriers)
                if return_statement == 'Error':
                    return ['Error', method_string, success_statements]
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
                if return_statement == 'Error':
                    return ['Error', method_string, success_statements]
                method_strings.append(method_string)
            elif 'MOVE_LiHa' in item[0]:
                return_statement, method_string = move_liha()
                if return_statement == 'Error':
                    return ['Error', method_string, success_statements]
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
        contents_to_prepare = queue_document['containers'][plate_operations[operation_number]['container']]['contents']
        contents_to_prepare = initial_plate_details['contents']
        contents_formatted = {}
        for well_key in initial_plate_details['contents'].keys():
            if well_key in contents_to_prepare.keys():
                if 'internal_standards' in additional_prep_details.keys() and additional_prep_details['internal_standards'] == True:
                    potential_volume = (additional_prep_details['extraction_iterations']-1) * \
                                       additional_prep_details['filter_volume_per_iteration'] + \
                                       additional_prep_details['iteration1_volume'] + \
                                       contents_to_prepare[well_key]['total_volume']
                    volume_dispensed = additional_prep_details['extraction_iterations'] * \
                                       additional_prep_details['filter_volume_per_iteration']
                    internal_standard_volume = volume_dispensed
                    internal_standard_identity = solvent_information['final_solvent_details']['internal_standard']
                    internal_standard_concentration =solvent_information['final_solvent_details']['internal_standard_conc']
                    contents_formatted[well_key] = {'plate_well': well_key,
                                                 'reagents': contents_to_prepare[well_key]['reagents'],
                                                 'solvents': [[additional_prep_details['liha_solvent_to_dispense'], volume_dispensed]],
                                                 'total_volume': internal_standard_volume,
                                                 'potential_volume': potential_volume,
                                                 'internal_standard': internal_standard_identity,
                                                 'internal_standard_volume': internal_standard_volume,
                                                 'internal_standard_concentration': internal_standard_concentration,
                                                 'target_product': contents_to_prepare[well_key]['target_product'],
                                                 'confirmation': 'none'}
                else:
                    accumulated_solvents = []
                    previous_well_solvents = contents_to_prepare[well_key]['solvents']
                    for previous_solvent in previous_well_solvents:
                        if previous_solvent[0] in ['dmf', 'dmso']:
                            accumulated_solvents.append(previous_solvent)
                    for destination in solvent_information['destination']:
                        if destination[0] == well_key:
                            volume_dispensed = additional_prep_details['iteration1_volume'] + \
                                                (additional_prep_details['extraction_iterations']-1) * \
                                                additional_prep_details['filter_volume_per_iteration']
                            accumulated_solvents.append([solvent_information['initial_solvent_details']['chemical_name'],
                                                         volume_dispensed])
                    volume_extracted = additional_prep_details['extraction_iterations'] * \
                                       additional_prep_details['filter_volume_per_iteration']
                    total_accumulated_volume = sum([item[1] for item in accumulated_solvents])
                    new_solvent_ratio = volume_extracted / total_accumulated_volume
                    final_well_solvents = [[item[0], item[1] * new_solvent_ratio] for item in accumulated_solvents]
                    contents_formatted[well_key] = {'plate_well': well_key,
                                                     'reagents': contents_to_prepare[well_key]['reagents'],
                                                     'solvents': final_well_solvents,
                                                     'total_volume': volume_extracted,
                                                     'target_product': contents_to_prepare[well_key]['target_product'],
                                                     'confirmation': 'none'}
            else:
                contents_formatted[well_key] = initial_plate_details['contents'][well_key]
                contents_formatted[well_key]['potential_volume'] = initial_plate_details['contents'][well_key]['total_volume']
            for key in contents_to_prepare[well_key].keys():
                if key not in contents_formatted[well_key].keys():
                    contents_formatted[well_key][key] = contents_to_prepare[well_key][key]
        pprint(contents_formatted)

        full_library_mongo['wellplates'][str(final_plate_details['_id'])]['location'][1] = final_plate_details['location'][1]
        full_library_mongo['wellplates'][str(final_plate_details['_id'])]['contents'] = contents_formatted
        full_library_mongo['wellplates'][str(final_plate_details['_id'])]['container_category'] = additional_prep_details[
            'target_container']
        success_statements.append('The completed extraction well-plate is located at: %s' % final_plate_details['location'][1])
        
        if __name__ == "__main__":
            pass
            return ['testing!', success_statements]


        platform_mongodb = pymongo.MongoClient(platform_mongodb_address, platform_mongodb_port)
        for library_key in full_library_mongo.keys():
            previous = initial_full_library[library_key]
            updated = full_library_mongo[library_key]
            return_statements = update_database(previous, updated, library_key, platform_mongodb)
    
        updated_queue_document = copy.deepcopy(queue_document)
        if __name__ == "__main__":
            updated_queue_document['operations'][operation_number]['completed'] = 'yes'
        reaction_plate_name = final_plate_details['container_name']
        target_container = additional_prep_details['target_container']
        updated_queue_document['containers'][target_container]['container_name'] = reaction_plate_name
        updated_queue_document['containers'][target_container]['contents'] = contents_formatted
        flag_return, return_statements = update_database_document(queue_document, updated_queue_document, 'queue',
                                                                 platform_mongodb)
        platform_mongodb.close()
        success_statements.append('Updated the platform library to be up-to-date following this run')
        return ['Success', schedule_path_output, success_statements]
    except Exception:
        return ['Error', traceback.format_exc(), ['']]


if __name__ == "__main__":
    queuename = 'scaffold_set_run_20220312_1_20220423'
    prep_return = liquid_liquid_extraction(queuename, '9')
    pprint(prep_return)
