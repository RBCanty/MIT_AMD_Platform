# -*- coding: utf-8 -*-
"""
Updated on 7/28/2021 --- MongoDB containing version of copy wellplate

@author: Brent
"""

import copy
import datetime
import json
import math
import pymongo
import yaml
from pprint import pprint
import sys
from mergedeep import merge
import traceback

if __name__ == "__main__":
    from supporting_functions import get_queue_information, convert_to_method_strings
    from platform_library_functions import update_database, update_database_document, query_document
    from useful_other_functions import solvent_finder, find_accessible_open_location, liha_grouping, find_labware
    from worktable_cleaner import initial_worktable_prep, worktable_export
    from tecan_functions import find_mca_tips
else:
    from Evoware_API.supporting_functions import get_queue_information, convert_to_method_strings
    from Evoware_API.platform_library_functions import update_database, update_database_document, query_document
    from Evoware_API.useful_other_functions import solvent_finder, find_accessible_open_location, liha_grouping, find_labware
    from Evoware_API.worktable_cleaner import initial_worktable_prep, worktable_export
    from Evoware_API.tecan_functions import find_mca_tips


def copy_wellplate(queue_name, operation_number):
    success_statements = []
    try:
        # ----------------------------------------------------------------------------
        # First we start with opening the queue document and getting operation info.
        queue_information = get_queue_information(queue_name, operation_number, 'copy_wellplate')
        if queue_information[0] != 'Success':
            return ['Error', queue_information[1], success_statements]
        else:
            return_statement, full_library_mongo, queue_document, additional_details, statements = queue_information
            
        success_statements.extend(statements)
        initial_full_library = copy.deepcopy(full_library_mongo)
        additional_prep_details = additional_details['additional_prep_details']
        
        # Then build the current worktable dictionary for building output files, we are
        # also making a duplicate of the worktable so it can be updated in this code
        return_statement, start_carriers, carriers_labware = initial_worktable_prep(additional_details['carriers_labware_filepath'], 
                                                                                    additional_details['empty_worktable_filepath'], 
                                                                                    full_library_mongo)
        if return_statement != 'Success':
            return ['Error', start_carriers, success_statements]
        updated_carriers = copy.deepcopy(start_carriers)
        success_statements.append('Loaded the carrier and labware library file!')
        # ----------------------------------------------------------------------------
        pprint(updated_carriers)
        # We need to find the initial wellplate for copying
        simplified_schedule = []
        occupied_locations = []
        allowed_prep_carriers = [32, 38, 44]
        inert_heater_shaker_carrier = [25]
        initial_plate_details = {'plate_type': additional_details['target_well_plate_type'],
                                 'container_name': additional_details['wellplate_name'],
                                 'queue_container_name': additional_details['target_container']}
        if initial_plate_details['container_name'] is None:
            return ['Error', 'No initial wellplate specified for "%s"' % additional_details['target_container'], success_statements]
        return_statement, extra_initial_plate_details = find_labware(updated_carriers, full_library_mongo, 
                                                                       initial_plate_details['container_name'], 
                                                                       initial_plate_details['plate_type'])
        if return_statement != 'Success':
            return ['Error', extra_initial_plate_details, success_statements]
        merge(initial_plate_details, extra_initial_plate_details)
        success_statements.append('Found the initial wellplate %s at %s' % (initial_plate_details['container_name'], 
                                                                            initial_plate_details['location'][1]))
        # Check to make sure the the initial plate is in an accessible location
        # TODO: Check to see if the wellplate is on a heater-shaker and if it is turn off shaking if the queue says so
        if initial_plate_details['location'][1][0] in allowed_prep_carriers:
            initial_plate_details['prep_location'] = initial_plate_details['location'][1]
        elif initial_plate_details['location'][1][0] in inert_heater_shaker_carrier:
            # TODO Add heater shaker details here
            initial_plate_details['prep_location'] = initial_plate_details['location'][1]
        else:
            return_statement, open_location = find_accessible_open_location(updated_carriers)
            if return_statement == 'Error':
                return ['Error', open_location, success_statements]
            initial_plate_details['prep_location'] = open_location
        occupied_locations.append(initial_plate_details['prep_location'])
        success_statements.append('Found a suitable prep location at location: %s' % initial_plate_details['prep_location'])
        if initial_plate_details['location'][1] != initial_plate_details['prep_location']:
            simplified_schedule.append(['TRANSFER_LABWARE', initial_plate_details['container_name'],
                                        initial_plate_details['location'][1], initial_plate_details['prep_location'],
                                        str(initial_plate_details['_id']), initial_plate_details['collection']])
        success_statements.append('Found a suitable prep location plate %s at %s' % (
            initial_plate_details['container_name'], initial_plate_details['prep_location']))
        
        # Now we work on the final well-plate to use on the platform (filtrate plate)
        final_plate_details = {'plate_type': additional_details['queue_containers'][additional_prep_details['target_container']]['plate_type'],
                               'container_name': additional_details['queue_containers'][additional_prep_details['target_container']]['container_name'],
                               'queue_container_name': additional_prep_details['target_container']}
        return_statement, extra_final_plate_details = find_labware(updated_carriers, full_library_mongo, 
                                                                   final_plate_details['container_name'], 
                                                                   final_plate_details['plate_type'])
        if return_statement != 'Success':
            return ['Error', extra_final_plate_details, success_statements]
        merge(final_plate_details, extra_final_plate_details)
        final_plate_details['wellplate_key'] = ['wellplates', str(final_plate_details['_id'])]
        success_statements.append('Found the filtrate plate %s at %s' % (final_plate_details['container_name'], 
                                                                         final_plate_details['location'][1]))
        # Check to make sure the the final plate is in an accessible location
        if final_plate_details['location'][1][0] in allowed_prep_carriers:
            final_plate_details['prep_location'] = final_plate_details['location'][1]
        else:
            return_statement, accessible_location = find_accessible_open_location(updated_carriers, exceptions = occupied_locations)
            if return_statement == 'Error':
                return ['Error', accessible_location, success_statements]
            final_plate_details['prep_location'] = accessible_location
        success_statements.append('Found a suitable prep location at location: %s' % final_plate_details['prep_location'])
        if final_plate_details['location'][1] != final_plate_details['prep_location']:
            simplified_schedule.append(['TRANSFER_LABWARE', final_plate_details['container_name'],
                                        final_plate_details['location'][1], final_plate_details['prep_location'],
                                        str(final_plate_details['_id']), 'wellplates'])
        success_statements.append('Found a suitable prep location plate %s at %s' % (
            final_plate_details['container_name'], final_plate_details['prep_location']))
        
        # Change the below to be consistent with taking an aliquot (volumes instead of dilutions)
        # With the plate tracking out of the way, we need to build the rest of the schedule
        if 'vol_to_add_initial' in additional_details.keys() and additional_prep_details['vol_to_add_initial'] != 0:
            solvent_to_find = additional_prep_details['makeup_solvent']
            relevant_wells = list(initial_plate_details['contents'].keys())
            solvents = {solvent_to_find: {'destination': [[relevant_well, additional_prep_details['vol_to_add_initial'], '', 1]
                                                          for relevant_well in relevant_wells]}}
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
                solvent_information['well_location']]['volume_ul'] = \
                solvent_information['final_solvent_details']['volume_ul']
            solvent_information['library_details'] = full_library_mongo['solvents'][solvent_information['container_id']]
    
            solvent_origin = solvent_information['library_details']['location']
            origin_labware_type = carriers_labware['labware'][solvent_information['library_details']['labware_type']]
            destination_labware_type = carriers_labware['labware'][initial_plate_details['plate_type']]
    
            simplified_schedule.append(['TIP_RINSE', 10, 10])
            simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 250])
            simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 250])
            simplified_schedule.append(['TIP_RINSE', 3, 3])
            return_statement, liha_groups = liha_grouping(solvents[solvent_to_find]['destination'], origin_labware_type,
                                        destination_labware_type)
            if return_statement == 'Error':
                return ['Error', liha_groups, success_statements]
            for group_no, liha_group in enumerate(liha_groups.keys()):
                group = liha_groups[liha_group]
                wells, volumes, tips = map(list, zip(*group))
                simplified_schedule.append(
                    ['LiHa_ASPIRATE_DISPENSE', solvent_origin[1], initial_plate_details['prep_location'], tips, volumes, wells])
            simplified_schedule.append(['MOVE_LiHa'])
            simplified_schedule.append(['TIP_RINSE', 5, 5])
            simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 100])
            simplified_schedule.append(['TIP_RINSE', 3, 4])
            success_statements.append('Added the LiHa solvent dispensing functions for initial plate prep')
        
        mca_volume = None
        
        if 'dilution_factor' in additional_prep_details.keys() and additional_prep_details['dilution_factor'] > 1:
            target_volume = additional_prep_details['target_volume']
            dilution_volume_needed = target_volume * (1 - (1 / additional_prep_details['dilution_factor']))
            mca_volume = dilution_volume_needed
            solvent_to_find = additional_prep_details['makeup_solvent']
            relevant_wells = list(initial_plate_details['contents'].keys())
            solvents = {solvent_to_find: {'destination': [[relevant_well, dilution_volume_needed, '', 1]
                                                          for relevant_well in relevant_wells]}}
            return_statement, solvent_information = solvent_finder(full_library_mongo, solvent_to_find, solvents[solvent_to_find])
            if return_statement == 'Error':
                return ['Error', solvent_information, success_statements]
            full_library_mongo['solvents'][solvent_information['container_id']]['contents'][
                solvent_information['well_location']]['volume_ul'] -= \
                solvent_information['final_solvent_details']['volume_ul']
            solvent_information['library_details'] = full_library_mongo['solvents'][solvent_information['container_id']]
    
            solvent_origin = solvent_information['library_details']['location']
            origin_labware_type = carriers_labware['labware'][solvent_information['library_details']['labware_type']]
            destination_labware_type = carriers_labware['labware'][final_plate_details['plate_type']]
    
            simplified_schedule.append(['TIP_RINSE', 10, 10])
            simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 250])
            simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 250])
            return_statement, liha_groups = liha_grouping(solvents[solvent_to_find]['destination'], origin_labware_type,
                                        destination_labware_type)
            if return_statement == 'Error':
                return ['Error', liha_groups, success_statements]
            for group_no, liha_group in enumerate(liha_groups.keys()):
                group = liha_groups[liha_group]
                wells, volumes, tips = map(list, zip(*group))
                simplified_schedule.append(
                    ['LiHa_ASPIRATE_DISPENSE', solvent_origin[1], final_plate_details['prep_location'], tips, volumes, wells])
            simplified_schedule.append(['MOVE_LiHa'])
            simplified_schedule.append(['TIP_RINSE', 5, 5])
            simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 100])
            success_statements.append('Added the LiHa solvent dispensing functions for initial plate prep')
        
        if 'aliquot_volume' in  additional_prep_details.keys() and additional_prep_details['aliquot_volume'] != 0:
            target_volume = additional_prep_details['target_volume']
            aliquot_volume = additional_prep_details['aliquot_volume']
            mca_volume = aliquot_volume
            additional_volume = target_volume - aliquot_volume
            if additional_volume < 0:
                return ['Error', 'Problem with volume specified: %s' % additional_volume, success_statements]
            solvent_to_find = additional_prep_details['solvent_to_use']
            relevant_wells = list(initial_plate_details['contents'].keys())
            solvents = {solvent_to_find: {'destination': [[relevant_well, additional_volume, '', 1]
                                                          for relevant_well in relevant_wells]}}
            return_statement, solvent_information = solvent_finder(full_library_mongo, solvent_to_find, solvents[solvent_to_find])
            if return_statement == 'Error':
                return ['Error', solvent_information, success_statements]
            full_library_mongo['solvents'][solvent_information['container_id']]['contents'][
                solvent_information['well_location']]['volume_ul'] = \
                solvent_information['final_solvent_details']['volume_ul']
            solvent_information['library_details'] = full_library_mongo['solvents'][solvent_information['container_id']]
            solvent_origin = solvent_information['library_details']['location']
            origin_labware_type = solvent_information['library_details']['labware_type']
            destination_labware_type = final_plate_details['plate_type']
    
            simplified_schedule.append(['TIP_RINSE', 10, 10])
            simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 250])
            simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 250])
            return_statement, liha_groups = liha_grouping(solvents[solvent_to_find]['destination'], 
                                                          carriers_labware,
                                                          origin_labware_type,
                                                          destination_labware_type)
            if return_statement == 'Error':
                return ['Error', liha_groups, success_statements]
            for group_no, liha_group in enumerate(liha_groups.keys()):
                group = liha_groups[liha_group]
                wells, volumes, tips = map(list, zip(*group))
                simplified_schedule.append(
                    ['LiHa_ASPIRATE_DISPENSE', solvent_origin[1], final_plate_details['prep_location'], tips, volumes, wells])
            simplified_schedule.append(['MOVE_LiHa'])
            simplified_schedule.append(['TIP_RINSE', 5, 5])
            simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 100])
            success_statements.append('Added the LiHa solvent dispensing functions for initial plate prep')
        
        if not mca_volume:
            return ['Error', 'Invalid MCA volume, workflow not specified', success_statements]
        # Now we make the transfer to the other well-plate with the MCA
        # Find suitable MCA tips on the platform
        return_statement, suitable_mca_tips, consumable_id = find_mca_tips(updated_carriers, full_library_mongo)
        if return_statement == 'Error':
            return ['Error', suitable_mca_tips, success_statements]
        simplified_schedule.append(['MCA_GET_TIPS', suitable_mca_tips])
        simplified_schedule.append(['MCA_MOVE'])
        simplified_schedule.append(['LARGE_WASTE_DUMP', suitable_mca_tips])
        full_library_mongo['consumables'][consumable_id]['tip_locations'].remove(str(suitable_mca_tips[1]))
        
        mca_asp_disp_divs = math.ceil(mca_volume / 200)
        
        if 'aliquot_volume' in additional_prep_details.keys():
            simplified_schedule.append(['EXECUTE_VBS_SCRIPT', 'C:\Scripts\datetime_output.vbs'])
            success_statements.append('Checking time for aliquot>%s>%s' % (queue_name, final_plate_details['container_name']))
        
        for mca_div in range(0, mca_asp_disp_divs):
            simplified_schedule.append(
                ['MCA_ASPIRATE', initial_plate_details['prep_location'], mca_volume / mca_asp_disp_divs])
            simplified_schedule.append(
                ['MCA_DISPENSE', final_plate_details['prep_location'], mca_volume / mca_asp_disp_divs])
        simplified_schedule.append(['MCA_MIX', final_plate_details['prep_location'], 200])
        simplified_schedule.append(['MCA_DROP_TIPS'])
        success_statements.append('Completed operations that use the MCA')
        
        # TODO: Add back the shaking here!
        if final_plate_details['location'][1] != final_plate_details['prep_location']:
            simplified_schedule.append(['TRANSFER_LABWARE', final_plate_details['container_name'],
                                        final_plate_details['prep_location'], final_plate_details['location'][1],
                                        str(final_plate_details['_id']), 'wellplates'])
        else:
            pass
            # full_library_mongo['wellplates'][final_plate_details['_id']]['location'][1] = final_plate_details['location'][1]
        # simplified_schedule.append(['TRANSFER_WELLPLATE'])    
        
        # Now we need to rebuild the worktable so we can write a file to run
        return_statement, filled_worktable_export = worktable_export(additional_details['carriers_labware_filepath'], 
                                                                     additional_details['empty_worktable_filepath'], 
                                                                     updated_carriers,
                                                                     full_library_mongo,
                                                                     simplified_schedule)
        if return_statement == 'Error':
            return ['Error', filled_worktable_export, success_statements]
        
        # We need to convert the simplified schedule into Tecan method strings
        return_statement, method_strings, updated_carriers = convert_to_method_strings(simplified_schedule, updated_carriers, full_library_mongo, carriers_labware)
        if return_statement == 'Error':
            return ['Error', method_strings, success_statements]
        
        
        with open(additional_details['schedule_path_output'], 'w', encoding='latin1') as outfile:
            for line in filled_worktable_export:
                outfile.write(line + '\n')
            for line in method_strings:
                outfile.write(line + '\n')
    
        contents_to_prepare = full_library_mongo['wellplates'][initial_plate_details['_id']]['contents']
        contents_formatted = {}
        for key in contents_to_prepare.keys():
            if 'dilution_factor' in additional_prep_details.keys():
                contents_formatted[key] = {'plate_well': key,
                                           'reagents': contents_to_prepare[key]['reagents'],
                                           'solvents': contents_to_prepare[key]['solvents'],
                                           'total_volume': contents_to_prepare[key]['total_volume'],
                                           'dilution_factor': additional_prep_details['dilution_factor']}
            elif 'aliquot_volume' in additional_prep_details.keys():
                contents_formatted[key] = {'plate_well': key,
                                           'reagents': contents_to_prepare[key]['reagents'],
                                           'solvents': contents_to_prepare[key]['solvents'],
                                           'total_volume': contents_to_prepare[key]['total_volume'],
                                           'dilution_factor': additional_prep_details['aliquot_volume']/contents_to_prepare[key]['total_volume']}
            for detail_key in contents_to_prepare[key].keys():
                if detail_key not in contents_formatted[key].keys():
                    contents_formatted[key][detail_key] = contents_to_prepare[key][detail_key]

        full_library_mongo['wellplates'][final_plate_details['_id']]['contents'] = contents_formatted
        full_library_mongo['wellplates'][final_plate_details['_id']]['container_category'] = additional_prep_details[
            'target_container']
        
        if __name__ == "__main__":
            pass
            return ['Testing', success_statements]
        
        platform_mongodb = pymongo.MongoClient(additional_details['platform_mongodb_address'], 
                                               additional_details['platform_mongodb_port'])
        for library_key in full_library_mongo.keys():
            previous = initial_full_library[library_key]
            updated = full_library_mongo[library_key]
            return_statements = update_database(previous, updated, library_key, platform_mongodb)
    
        updated_queue_document = copy.deepcopy(queue_document)
        
        if 'queues_to_update' in additional_prep_details.keys() and additional_prep_details:
            for queue_to_update in additional_prep_details['queues_to_update']:
                if queue_to_update == queue_name:
                    continue
                document_return, return_statements = query_document('queue', queue_to_update, 'queue_name', platform_mongodb)
                if not type(document_return) == dict:
                    return ['Error', document_return, success_statements]
                original_queue = copy.deepcopy(document_return)
                container_details = document_return['containers'][additional_prep_details['target_container']]
                container_details['container_name'] = final_plate_details['container_name']
                container_details['contents'] = contents_formatted
                if additional_prep_details['target_container'] == 'aliquot_plate_1':
                    document_return['status'] = 'idle'
                flag_return, return_statements = update_database_document(original_queue, document_return, 'queue',
                                                                  platform_mongodb)
                # TODO: If error...
                
        if __name__ == "__main__":
            updated_queue_document['operations'][operation_number]['completed'] = 'yes'
        updated_queue_document['containers'][additional_prep_details['target_container']]['container_name'] = \
            final_plate_details['container_name']
        updated_queue_document['containers'][additional_prep_details['target_container']]['contents'] = contents_formatted
        flag_return, return_statements = update_database_document(queue_document, updated_queue_document, 'queue',
                                                                  platform_mongodb)
        platform_mongodb.close()
    
        return ['Success', additional_details['schedule_path_output'], success_statements]
    
    except Exception:
        return ['Error', traceback.format_exc(), success_statements]

if __name__ == "__main__":
    queuename = 'real_dual_cat_test_20221118_1'
    prep_return = copy_wellplate(queuename, '19')
    pprint(prep_return)
