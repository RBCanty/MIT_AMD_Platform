# -*- coding: utf-8 -*-
"""
Updated on 4/11/2022 --- MongoDB containing version of filter plate sequence

@author: Brent
"""

import copy
import math
import traceback
import pymongo
from mergedeep import merge
from pprint import pprint
from datetime import datetime
import sys
from copy import deepcopy

if __name__ == "__main__":
    from supporting_functions import get_queue_information, convert_to_method_strings
    from platform_library_functions import update_database, update_database_document, move_database_document
    from tecan_functions import find_mca_tips
    from useful_other_functions import solvent_finder_origin, liha_grouping, find_accessible_open_location, \
        check_for_problems, find_labware, evaporation_estimation
    from worktable_cleaner import initial_worktable_prep, worktable_export
else:
    from Evoware_API.supporting_functions import get_queue_information, convert_to_method_strings
    from Evoware_API.platform_library_functions import update_database, update_database_document, move_database_document
    from Evoware_API.tecan_functions import find_mca_tips
    from Evoware_API.useful_other_functions import solvent_finder_origin, liha_grouping, find_accessible_open_location,\
        check_for_problems, find_labware, evaporation_estimation
    from Evoware_API.worktable_cleaner import initial_worktable_prep, worktable_export


def filter_plate_sequence(queue_name, operation_number):
    success_statements = []
    try:
        # ----------------------------------------------------------------------------
        # First we start with opening the queue document and getting operation info.
        queue_information = get_queue_information(queue_name, operation_number, 'filter_wellplate')
        if queue_information[0] != 'Success':
            return ['Error', 'Problem: %s' % queue_information[1], success_statements]
        else:
            return_statement, full_library_mongo, queue_document, additional_details, statements = queue_information

        success_statements.extend(statements)
        initial_full_library = copy.deepcopy(full_library_mongo)
        additional_prep_details = additional_details['additional_prep_details']

        # Then build the current worktable dictionary for building output files, we are
        # also making a duplicate of the worktable so it can be updated in this code
        return_statement, start_carriers, carriers_labware = initial_worktable_prep(
            additional_details['carriers_labware_filepath'], additional_details['empty_worktable_filepath'],
            full_library_mongo)
        if return_statement != 'Success':
            return ['Error', 'Problem: %s' % start_carriers, success_statements]
        updated_carriers = copy.deepcopy(start_carriers)
        # pprint(updated_carriers, width=120)
        success_statements.append('Loaded the carrier and labware library file!')
        # ----------------------------------------------------------------------------

        # We need to start by looking for the wellplates that are needed for filtration,
        # starting with the initial wellplate
        simplified_schedule = []
        occupied_locations = []
        allowed_prep_carriers = [32, 38, 44]
        initial_plate_details = {'plate_type': additional_details['target_well_plate_type'],
                                 'container_name': additional_details['wellplate_name'],
                                 'queue_container_name': additional_details['target_container']}
        if initial_plate_details['container_name'] is None:
            return ['Error', 'Problem: no wellplate specified: %s' % additional_details['target_container'],
                    success_statements]
        return_statement, extra_initial_plate_details = find_labware(updated_carriers, full_library_mongo,
                                                                     initial_plate_details['container_name'],
                                                                     initial_plate_details['plate_type'],
                                                                     initial_plate_details['queue_container_name'])
        if return_statement != 'Success':
            return ['Error', 'Problem: %s' % extra_initial_plate_details, success_statements]
        merge(initial_plate_details, extra_initial_plate_details)
        success_statements.append('Found the initial wellplate %s at: %s' % (initial_plate_details['container_name'],
                                                                             initial_plate_details['location'][1]))
        # Check to make sure the initial plate is in an accessible location
        if initial_plate_details['location'][1][0] in allowed_prep_carriers:
            initial_plate_details['prep_location'] = initial_plate_details['location'][1]
        else:
            return_statement, open_location = find_accessible_open_location(updated_carriers)
            if return_statement == 'Error':
                return ['Error', 'Problem: %s' % open_location, success_statements]
            initial_plate_details['prep_location'] = open_location
        occupied_locations.append(initial_plate_details['prep_location'])
        success_statements.append('Found a prep location for %s at: %s' % (initial_plate_details['container_name'],
                                                                           initial_plate_details['prep_location']))
        if initial_plate_details['location'][1] != initial_plate_details['prep_location']:
            simplified_schedule.append(['TRANSFER_LABWARE', initial_plate_details['container_name'],
                                        initial_plate_details['location'][1], initial_plate_details['prep_location'],
                                        str(initial_plate_details['_id']), initial_plate_details['collection']])

        # Now we work on the final well-plate to use on the platform (filtrate plate)
        final_plate_queue_details = additional_details['queue_containers'][additional_prep_details['target_container']]
        final_plate_details = {'plate_type': final_plate_queue_details['plate_type'],
                               'container_name': final_plate_queue_details['container_name'],
                               'queue_container_name': additional_prep_details['target_container']}
        return_statement, extra_final_plate_details = find_labware(updated_carriers, full_library_mongo,
                                                                   final_plate_details['container_name'],
                                                                   final_plate_details['plate_type'],
                                                                   final_plate_details['queue_container_name'])
        if return_statement != 'Success':
            return ['Error', 'Problem: %s' % extra_final_plate_details, success_statements]
        merge(final_plate_details, extra_final_plate_details)
        final_plate_details['wellplate_key'] = ['wellplates', str(final_plate_details['_id'])]
        success_statements.append('Found the filtrate plate %s at %s' % (final_plate_details['container_name'],
                                                                         final_plate_details['location'][1]))
        # The final well-plate will be put into the TeVac so no need to figure
        # out if it is bed position accessible or not

        # Now we are going to look for a filter plate that is ready to be used
        return_statement, filter_plate_details = find_labware(updated_carriers, full_library_mongo,
                                                              None, '96 Well Filtration Plate',
                                                              'solid_filter_plate')
        if return_statement != 'Success':
            return ['Error', 'Problem: %s' % filter_plate_details, success_statements]
        full_library_mongo['consumables'][str(filter_plate_details['_id'])]['consumable_status'] = 'used'
        success_statements.append('Found a suitable filter plate %s at %s' % (
            filter_plate_details['container_name'], filter_plate_details['location'][1]))
        if filter_plate_details['location'][1][0] in allowed_prep_carriers:
            filter_plate_details['prep_location'] = filter_plate_details['location'][1]
        else:
            return_statement, open_location = find_accessible_open_location(updated_carriers,
                                                                            exceptions=occupied_locations)
            if return_statement == 'Error':
                return ['Error', 'Problem: %s' % open_location, success_statements]
            filter_plate_details['prep_location'] = open_location
        occupied_locations.append(filter_plate_details['prep_location'])
        success_statements.append('Found a suitable prep location for filter plate %s at %s' % (
            filter_plate_details['container_name'], filter_plate_details['prep_location']))
        # The filter-plate requires a fancy lid handling transfer to work, this
        # will happen during the TeVac assembly process

        # Now look at the additional prep details from the queue document to see if a HPLC wellpate needs to be found
        if 'hplc_container' in additional_prep_details.keys():
            hplc_wellplate_entry = additional_details['queue_containers'][additional_prep_details['hplc_container']]
            hplc_wellplate_details = {'plate_type': hplc_wellplate_entry['plate_type'],
                                      'container_name': hplc_wellplate_entry['container_name'],
                                      'queue_container_name': additional_prep_details['hplc_container']}
            return_statement, extra_details = find_labware(updated_carriers, full_library_mongo,
                                                           hplc_wellplate_details['container_name'],
                                                           hplc_wellplate_details['plate_type'],
                                                           hplc_wellplate_details['queue_container_name'],
                                                           exceptions=[final_plate_details['container_name']])
            if return_statement != 'Success':
                return ['Error', 'Problem: %s' % extra_details, success_statements]
            merge(hplc_wellplate_details, extra_details)
            hplc_wellplate_details['wellplate_key'] = ['wellplates', str(hplc_wellplate_details['_id'])]
            success_statements.append('Found the hplc plate %s at %s' % (hplc_wellplate_details['container_name'],
                                                                         hplc_wellplate_details['location'][1]))
            if hplc_wellplate_details['location'][1][0] in allowed_prep_carriers:
                hplc_wellplate_details['prep_location'] = hplc_wellplate_details['location'][1]
            else:
                return_statement, open_location = find_accessible_open_location(updated_carriers,
                                                                                exceptions=occupied_locations)
                if return_statement == 'Error':
                    return ['Error', 'Problem: %s' % open_location, success_statements]
                hplc_wellplate_details['prep_location'] = open_location
            occupied_locations.append(hplc_wellplate_details['prep_location'])
            success_statements.append(
                'Found a suitable prep location at location: %s' % hplc_wellplate_details['prep_location'])
            if hplc_wellplate_details['location'][1] != hplc_wellplate_details['prep_location']:
                simplified_schedule.append(['TRANSFER_LABWARE', hplc_wellplate_details['container_name'],
                                            hplc_wellplate_details['location'][1],
                                            hplc_wellplate_details['prep_location'],
                                            str(hplc_wellplate_details['_id']), 'wellplates'])
        else:
            hplc_wellplate_details = {}

        # With the wellplate bookkeeping out of the way, we can start working on the methods
        # --------------------------------------------------------------------------------------------------------------
        # Start by finding the TeVac on the liquid handler bed
        tevac_location = []
        for grid in updated_carriers.keys():
            for carrier_id in updated_carriers[grid].keys():
                if updated_carriers[grid][carrier_id]['carrier_name'] == 'TeVacS Custom':
                    tevac_location = [grid, carrier_id]
        if not tevac_location:
            return ['Error', 'Problem: TeVacS was not found on the current platform worktable', success_statements]

        # Now we need to figure out where the filtration is taking place
        if 'tevac_location' not in additional_prep_details.keys():
            return ['Error', 'TeVac filtration location not specified']
        tevac_site = additional_prep_details['tevac_location']
        small_wellplates = ['96 Well Microplate']
        large_wellplates = ['96 Well DeepWell']
        
        if tevac_site == 'front' and final_plate_details['plate_type'] in large_wellplates:
            problem_return = check_for_problems(updated_carriers, [tevac_location[0], 5, tevac_location[1]],
                                                final_plate_details['plate_type'])
            if problem_return != 'Acceptable':
                return ['Error', 'Problem: Tevac %s is inaccessible, collision hazard' % tevac_site, success_statements]
            tevac_filterplate_location = [tevac_location[0], 4, tevac_location[1]]
            tevac_filtrateplate_location = [tevac_location[0], 5, tevac_location[1]]
            tevac_block_location = [tevac_location[0], 6, tevac_location[1]]
            tevac_block_transfer_location = [tevac_location[0], 7, tevac_location[1]]
        elif tevac_site == 'back' and final_plate_details['plate_type'] in small_wellplates:
            problem_return = check_for_problems(updated_carriers, [tevac_location[0], 2, tevac_location[1]],
                                                final_plate_details['plate_type'])
            if problem_return != 'Acceptable':
                return ['Error', 'Problem: Tevac %s is inaccessible, collision hazard' % tevac_site, success_statements]
            tevac_filterplate_location = [tevac_location[0], 1, tevac_location[1]]
            tevac_filtrateplate_location = [tevac_location[0], 2, tevac_location[1]]
            tevac_block_location = [tevac_location[0], 3, tevac_location[1]]
            tevac_block_transfer_location = [tevac_location[0], 8, tevac_location[1]]
        else:
            return ['Error', 'Problem: Invalid Tevac position requested: %s' % tevac_site, success_statements]

        # Add the method_strings to assemble the tevac filter
        simplified_schedule.append(['SPECIAL_TRANSFER_LABWARE', 'TeVac Block', tevac_block_location,
                                    tevac_block_transfer_location, '', ''])
        simplified_schedule.append(['TRANSFER_LABWARE', final_plate_details['container_name'],
                                    final_plate_details['location'][1], tevac_filtrateplate_location,
                                    str(final_plate_details['_id']), 'wellplates'])
        simplified_schedule.append(['SPECIAL_TRANSFER_LABWARE', 'TeVac Block', tevac_block_transfer_location,
                                    tevac_block_location, '', ''])
        simplified_schedule.append(['SPECIAL_LID_TRANSFER', filter_plate_details['container_name'],
                                    filter_plate_details['location'][1], filter_plate_details['prep_location'],
                                    tevac_filterplate_location, 'destination', str(filter_plate_details['_id']),
                                    filter_plate_details['collection']])
        # --------------------------------------------------------------------------------------------------------------
        # Now we need to look at each well in the well-plate to figure out what
        # solvent needs to be added and where, starting with the HPLC plate
        solvent_details = {}
        if 'initial_solvent_volume' in additional_prep_details.keys():
            solvent_to_find = additional_prep_details['initial_filter_solvent']
        if 'hplc_container_solvent' in additional_prep_details.keys():
            hplc_well_solvent = additional_prep_details['hplc_container_solvent']
        else:
            hplc_well_solvent = 'dmso'
        if 'hplc_prep_volumes' in additional_prep_details.keys():
            hplc_prep_volumes = additional_prep_details['hplc_prep_volumes']
        else:
            hplc_prep_volumes = {'final_product_transfer': 100, 'extra_final_product_volume': 0,
                                 'intermediate_product_transfer': 50, 'extra_intermediate_product_volume': 75,
                                 'analytical_product_transfer': 20, 'extra_analytical_product_volume': 50}
        hplc_plate_solvents_needed = {}
        hplc_plate_tranfers_needed = {}
        problem_wells = []
        for well_key in initial_plate_details['contents'].keys():
            well_details = initial_plate_details['contents'][well_key]
            product_status = well_details['final_product']
            if 'hplc_container' in additional_prep_details.keys():
                if product_status[0] == 0:
                    extra_volume_needed = hplc_prep_volumes['extra_final_product_volume']
                    transfer_key = 'final_product_transfer'
                elif product_status[0] != 0 and product_status[1] == 'yes':
                    extra_volume_needed = hplc_prep_volumes['extra_intermediate_product_volume']
                    transfer_key = 'intermediate_product_transfer'
                elif product_status[0] != 0 and product_status[1] == 'no':
                    extra_volume_needed = hplc_prep_volumes['extra_analytical_product_volume']
                    transfer_key = 'analytical_product_transfer'
                else:
                    problem_wells.append('Unable to handle %s : %s' % (well_key, product_status))
                    continue
                if extra_volume_needed != 0:
                    if hplc_well_solvent not in hplc_plate_solvents_needed.keys():
                        hplc_plate_solvents_needed[hplc_well_solvent] = {'destinations': []}
                    hplc_plate_solvents_needed[hplc_well_solvent]['destinations'].append(
                        [well_key, extra_volume_needed, '', 1])
                hplc_plate_tranfers_needed[well_key] = {'total': hplc_prep_volumes[transfer_key], 'transferred': 0}
        if len(problem_wells) != 0:
            return ['Error', 'Problem: %s' % problem_wells, success_statements]

        simplified_schedule.append(['TIP_RINSE', 10, 10])
        simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 150])
        simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 150])

        problematic_solvents = []
        if hplc_plate_solvents_needed:
            for solvent_to_find in hplc_plate_solvents_needed:
                return_statement, solvent_info = solvent_finder_origin(full_library_mongo, [[solvent_to_find, 1]], {})
                if return_statement != 'Success':
                    problematic_solvents.append(solvent_info)
                solvent_information = full_library_mongo[solvent_info[solvent_to_find]['collection']][
                                                            solvent_info[solvent_to_find]['container_id']]
                if solvent_info[solvent_to_find]['collection'] == 'solvents':
                    solvent_origin = solvent_information['location'][1]
                elif solvent_info[solvent_to_find]['collection'] == 'reagents':
                    solvent_origin = [solvent_information['location'][1]+[solvent_information['origin_well_location']]]
                else:
                    return ['Error', 'Problem: Unsuitable collection %s for solvents' %
                            solvent_info[solvent_to_find]['collection'], success_statements]
                mongo_solvent_name = solvent_info[solvent_to_find]['chemical_name']
                solvent_destinations = hplc_plate_solvents_needed[solvent_to_find]['destinations']
                if mongo_solvent_name not in solvent_details.keys():
                    solvent_details[mongo_solvent_name] = solvent_info[solvent_to_find]
                for solvent_destination in solvent_destinations:
                    solvent_details[mongo_solvent_name]['volume_needed'] += solvent_destination[1]
                return_statement, liha_groups = liha_grouping(solvent_destinations, carriers_labware,
                                                              solvent_information['labware_type'],
                                                              hplc_wellplate_details['plate_type'])
                if return_statement == 'Error':
                    return ['Error', 'Problem: %s' % liha_groups, success_statements]
                for group_no, liha_group in enumerate(liha_groups.keys()):
                    group = liha_groups[liha_group]
                    wells, volumes, tips = map(list, zip(*group))
                    simplified_schedule.append(['LiHa_ASPIRATE_DISPENSE', solvent_origin,
                                                hplc_wellplate_details['prep_location'], tips, volumes, wells,
                                                'Adjusted water free dispense'])
                if len(list(liha_groups.keys())) != 0:
                    simplified_schedule.append(['MOVE_LiHa'])
                    simplified_schedule.append(['TIP_RINSE', 5, 5])
                    simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 150])
                    simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 150])
                    simplified_schedule.append(['TIP_RINSE', 3, 4])

        # Now we need to figure out the solvents that need to be added to the original wellplate
        original_plate_solvents_needed = {}
        solvent_problems = []
        current_time = datetime.now()
        if 'initial_solvent_volume' in additional_prep_details.keys():
            for well_key in list(initial_plate_details['contents'].keys()):
                well_details = initial_plate_details['contents'][well_key]
                previous_time = datetime.strptime(well_details['date_updated'], '%m/%d/%Y %H:%M:%S')
                elapsed_time = (current_time - previous_time).total_seconds()
                return_statement, solvent_return = evaporation_estimation(well_details['solvents'],
                                                                          initial_plate_details['labware_type'],
                                                                          elapsed_time, 80)
                product_status = well_details['final_product']
                if 'hplc_container' in additional_prep_details.keys() and product_status[0] == 0:
                    solvent_to_find = additional_prep_details['hplc_container_solvent']
                else:
                    solvent_to_find = additional_prep_details['initial_filter_solvent']
                if solvent_to_find not in original_plate_solvents_needed.keys():
                    original_plate_solvents_needed[solvent_to_find] = {'destinations': []}
                evaporated_volume = 0
                dmf_dmso_volume = 0
                for solvent in well_details['solvents']:
                    if solvent[0].lower() not in ['dmf', 'dmso', 'water']:
                        if 'CN(C)C=O' not in solvent[0] and 'CS(C)=O' not in solvent[0]:
                            evaporated_volume += solvent[1]
                    if solvent[0] in ['dmf', 'dmso']:
                        dmf_dmso_volume += solvent[1]
                    if 'CN(C)C=O' in solvent[0] or 'CS(C)=O' in solvent[0]:
                        dmf_dmso_volume += solvent[1]
                if dmf_dmso_volume > 0:
                    needed_volume = additional_prep_details['initial_solvent_volume'] - dmf_dmso_volume
                    if needed_volume <= 1:
                        continue
                    original_plate_solvents_needed[solvent_to_find]['destinations'].append(
                        [well_key, needed_volume, '', 1])
                else:
                    needed_volume = additional_prep_details['initial_solvent_volume'] - (
                                well_details['total_volume'] - evaporated_volume)
                    if needed_volume <= 1:
                        continue
                    original_plate_solvents_needed[solvent_to_find]['destinations'].append(
                        [well_key, needed_volume, '', 1])

            for solvent_to_find in original_plate_solvents_needed.keys():
                return_status, solvent_info = solvent_finder_origin(full_library_mongo, [[solvent_to_find, 1]], {})
                if return_status != 'Success':
                    solvent_problems.append(['Issue with %s, %s' % (solvent_to_find, solvent_info)])
                    continue
                mongo_solvent_name = solvent_info[solvent_to_find]['chemical_name']
                solvent_destinations = original_plate_solvents_needed[solvent_to_find]['destinations']
                if mongo_solvent_name not in solvent_details.keys():
                    solvent_details[mongo_solvent_name] = solvent_info[solvent_to_find]
                for solvent_destination in solvent_destinations:
                    solvent_details[mongo_solvent_name]['volume_needed'] += solvent_destination[1]
                if solvent_details[mongo_solvent_name]['volume_needed'] > solvent_details[mongo_solvent_name]['solvent_volume']:
                    solvent_problems.append('Problem: Not enough volume of %s' % mongo_solvent_name)
                solvent_information = full_library_mongo[solvent_info[solvent_to_find]['collection']][
                    solvent_info[solvent_to_find]['container_id']]
                return_statement, liha_groups = liha_grouping(solvent_destinations, carriers_labware,
                                                              solvent_information['labware_type'],
                                                              initial_plate_details['plate_type'])
                if return_statement == 'Error':
                    return ['Error', 'Problem: %s' % liha_groups, success_statements]
                if solvent_info[solvent_to_find]['collection'] == 'solvents':
                    solvent_origin = solvent_information['location'][1]
                else:
                    solvent_origin = solvent_information['location'][1] + \
                                     [solvent_info[solvent_to_find]['origin_well_location']]
                for group_no, liha_group in enumerate(liha_groups.keys()):
                    group = liha_groups[liha_group]
                    wells, volumes, tips = map(list, zip(*group))
                    simplified_schedule.append(['LiHa_ASPIRATE_DISPENSE', solvent_origin,
                                                initial_plate_details['prep_location'], tips, volumes, wells,
                                                'Adjusted water free dispense Filter Prep'])
                if len(list(liha_groups.keys())) != 0:
                    simplified_schedule.append(['MOVE_LiHa'])
                    simplified_schedule.append(['TIP_RINSE', 5, 5])
                    simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
                    simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
                    simplified_schedule.append(['TIP_RINSE', 3, 4])

            if solvent_problems:
                return ['Error', 'Problem: Issue with solvents: %s' % solvent_problems, success_statements]
            success_statements.append('Added the LiHa solvent dispensing functions for initial plate prep')
        else:
            simplified_schedule.append(['MOVE_LiHa'])
        
        # Find suitable MCA tips on the platform
        return_statement, suitable_mca_tips, consumable_id = find_mca_tips(updated_carriers, full_library_mongo)
        if return_statement == 'Error':
            return ['Error', suitable_mca_tips, success_statements]
        simplified_schedule.append(['MCA_GET_TIPS', suitable_mca_tips])
        simplified_schedule.append(['MCA_MOVE'])
        simplified_schedule.append(['LARGE_WASTE_DUMP', suitable_mca_tips])
        full_library_mongo['consumables'][consumable_id]['tip_locations'].remove(str(suitable_mca_tips[1]))
        if additional_prep_details['mca_mix']:
            mix_volume = 200
            simplified_schedule.append(['MCA_MIX', initial_plate_details['prep_location'], mix_volume, 50])
            simplified_schedule.append(['MCA_MOVE'])

        # Now it is time to turn on the Tevac
        simplified_schedule.append(['TEVAC_FUNCTION', ['vacuum_on', additional_prep_details['tevac_location']], 300])
        mca_volume = additional_prep_details['mca_filter_volume']
        mca_divs = math.ceil(mca_volume / 200)
        for mca_div in range(0, mca_divs):
            simplified_schedule.append(['MCA_ASPIRATE', initial_plate_details['prep_location'], mca_volume / mca_divs])
            simplified_schedule.append(['MCA_DISPENSE', tevac_filterplate_location, mca_volume / mca_divs,
                                        'Adjusted Water free dispense DiTi 200 Clear TeVaC'])
        simplified_schedule.append(['MCA_DROP_TIPS'])
        if additional_prep_details['tevac_filter_time'] > 0:
            simplified_schedule.append(['START_TIMER', 1])
            simplified_schedule.append(['WAIT_FOR_TIMER', 1, additional_prep_details['tevac_filter_time']])

        simplified_schedule.append(['TEVAC_FUNCTION', ['vent', additional_prep_details['tevac_location']], ''])
        simplified_schedule.append(['TEVAC_FUNCTION', ['deactivate'], ''])

        # Now it is time to disassemble the TeVac
        simplified_schedule.append(['SPECIAL_LID_TRANSFER', filter_plate_details['container_name'],
                                    filter_plate_details['prep_location'], filter_plate_details['location'][1],
                                    tevac_filterplate_location, 'source', str(filter_plate_details['_id']),
                                    filter_plate_details['collection']])
        simplified_schedule.append(['SPECIAL_TRANSFER_LABWARE', 'TeVac Block', tevac_block_location,
                                    tevac_block_transfer_location, '', ''])
        
        # If we need to make an HPLC plate, then we need to divert the final plate
        # to an accessible location on the liquid handler bed
        if 'hplc_container' in additional_prep_details.keys():
            pprint(final_plate_details)
            if final_plate_details['location'][1][0] not in allowed_prep_carriers:
                return_statement, open_location = find_accessible_open_location(updated_carriers,
                                                                                exceptions=occupied_locations)
                if return_statement == 'Error':
                    return ['Error', open_location, success_statements]
                final_plate_details['prep_location'] = open_location
            else:
                final_plate_details['prep_location'] = final_plate_details['location'][1]
            
            simplified_schedule.append(['TRANSFER_LABWARE', final_plate_details['container_name'],
                                        tevac_filtrateplate_location, final_plate_details['prep_location'],
                                        str(final_plate_details['_id']), 'wellplates'])
            simplified_schedule.append(['SPECIAL_TRANSFER_LABWARE', 'TeVac Block', tevac_block_transfer_location,
                                        tevac_block_location, '', ''])
            # Now we need to get the sample transferred to the HPLC wellplate
            # Find suitable MCA tips on the platform
            return_statement, suitable_mca_tips, consumable_id = find_mca_tips(updated_carriers, full_library_mongo)
            if return_statement == 'Error':
                return ['Error', suitable_mca_tips, success_statements]
            simplified_schedule.append(['MCA_GET_TIPS', suitable_mca_tips])
            simplified_schedule.append(['MCA_MOVE'])
            simplified_schedule.append(['LARGE_WASTE_DUMP', suitable_mca_tips])
            full_library_mongo['consumables'][consumable_id]['tip_locations'].remove(str(suitable_mca_tips[1]))

            mca_transfer_amount = min([hplc_plate_tranfers_needed[key]['total'] for key in
                                       hplc_plate_tranfers_needed.keys()])
            mca_divs = math.ceil(mca_transfer_amount / 200)
            for mca_div in range(0, mca_divs):
                simplified_schedule.append(['MCA_ASPIRATE', final_plate_details['prep_location'],
                                            mca_transfer_amount / mca_divs])
                simplified_schedule.append(['MCA_DISPENSE', hplc_wellplate_details['prep_location'],
                                            mca_transfer_amount / mca_divs])
            mix_volume = 200
            simplified_schedule.append(['MCA_MIX', hplc_wellplate_details['prep_location'], mix_volume, 50])
            simplified_schedule.append(['MCA_MOVE'])
            simplified_schedule.append(['MCA_DROP_TIPS'])

            for well_key in hplc_plate_tranfers_needed.keys():
                hplc_plate_tranfers_needed[well_key]['transferred'] += mca_transfer_amount

            # Now for the LiHa specific transfers that are needed
            liha_tip_number = 1
            simplified_schedule.append(['MOVE_LiHa'])
            simplified_schedule.append(['TIP_RINSE', 5, 5])
            simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
            simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
            for well_key in hplc_plate_tranfers_needed.keys():
                total_needed = hplc_plate_tranfers_needed[well_key]['total']
                amount_transferred = hplc_plate_tranfers_needed[well_key]['transferred']
                amount_left = total_needed - amount_transferred
                if amount_left <= 0:
                    continue
                tips = [liha_tip_number]
                volumes = [amount_left]
                wells = [well_key]
                simplified_schedule.append(['LiHa_ASPIRATE_DISPENSE', final_plate_details['prep_location'] + [[well_key]],
                                            hplc_wellplate_details['prep_location'] + [[well_key]], tips, volumes,
                                            wells])
                hplc_plate_tranfers_needed[well_key]['transferred'] += amount_left
                if liha_tip_number == 8:
                    simplified_schedule.append(['TIP_RINSE', 5, 5])
                    simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
                    simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
                    liha_tip_number = 1
                else:
                    liha_tip_number += 1
            if hplc_plate_tranfers_needed and liha_tip_number != 1:
                simplified_schedule.append(['MOVE_LiHa'])
                simplified_schedule.append(['TIP_RINSE', 5, 5])
                simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
                simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
                simplified_schedule.append(['TIP_RINSE', 4, 4])
                liha_tip_number = 1
            if final_plate_details['prep_location'] != final_plate_details['location'][1]:
                simplified_schedule.append(['TRANSFER_LABWARE', final_plate_details['container_name'],
                                            final_plate_details['prep_location'], final_plate_details['location'][1],
                                            str(final_plate_details['_id']), 'wellplates'])
            if hplc_wellplate_details['prep_location'] != hplc_wellplate_details['location'][1]:
                simplified_schedule.append(['TRANSFER_LABWARE', hplc_wellplate_details['container_name'],
                                            hplc_wellplate_details['prep_location'], hplc_wellplate_details['location'][1],
                                            str(hplc_wellplate_details['_id']), 'wellplates'])
        else:  # Otherwise, complete the transfer as normal
            simplified_schedule.append(['TRANSFER_LABWARE', final_plate_details['container_name'],
                                        tevac_filtrateplate_location, final_plate_details['location'][1],
                                        str(final_plate_details['_id']), 'wellplates'])
            simplified_schedule.append(['SPECIAL_TRANSFER_LABWARE', 'TeVac Block', tevac_block_transfer_location,
                                        tevac_block_location, '',  ''])

        if initial_plate_details['prep_location'] != initial_plate_details['location'][1]:
            simplified_schedule.append(['TRANSFER_LABWARE', initial_plate_details['container_name'],
                                        initial_plate_details['prep_location'], initial_plate_details['location'][1],
                                        str(initial_plate_details['_id']), 'wellplates'])
            
        # Now we need to rebuild the worktable, so we can write a file to run
        return_statement, filled_worktable_export = worktable_export(additional_details['carriers_labware_filepath'],
                                                                     additional_details['empty_worktable_filepath'],
                                                                     updated_carriers, full_library_mongo,
                                                                     simplified_schedule)
        if return_statement == 'Error':
            return ['Error', filled_worktable_export, success_statements]
        
        # We need to convert the simplified schedule into Tecan method strings
        return_statement, method_strings, updated_carriers = convert_to_method_strings(simplified_schedule,
                                                                                       updated_carriers,
                                                                                       full_library_mongo,
                                                                                       carriers_labware)
        if return_statement == 'Error':
            return ['Error', method_strings, success_statements]

        with open(additional_details['schedule_path_output'], 'w', encoding='latin1') as outfile:
            for line in filled_worktable_export:
                outfile.write(line + '\n')
            for line in method_strings:
                outfile.write(line + '\n')

        # TODO The below code needs to be cleaned and updated!

        contents_to_prepare = initial_plate_details['contents']
        final_contents_formatted = {}
        hplc_contents_formatted = {}
        current_time = datetime.now()
        for well_key in initial_plate_details['contents'].keys():
            well_details = initial_plate_details['contents'][well_key]
            product_status = well_details['final_product']
            if 'hplc_container' in additional_prep_details.keys():
                hplc_container_solvent = additional_prep_details['hplc_container_solvent']
                if product_status[0] == 0:
                    extra_volume_needed = hplc_prep_volumes['extra_final_product_volume']
                    volume_needed = hplc_prep_volumes['final_product_transfer']
                    initial_filter_solvent = additional_prep_details['hplc_container_solvent']
                elif product_status[0] != 0 and product_status[1] == 'yes':
                    extra_volume_needed = hplc_prep_volumes['extra_intermediate_product_volume']
                    volume_needed = hplc_prep_volumes['intermediate_product_transfer']
                    initial_filter_solvent = additional_prep_details['initial_filter_solvent']
                elif product_status[0] != 0 and product_status[1] == 'no':
                    extra_volume_needed = hplc_prep_volumes['extra_analytical_product_volume']
                    volume_needed = hplc_prep_volumes['analytical_product_transfer']
                    initial_filter_solvent = additional_prep_details['initial_filter_solvent']
                original_well_solvents = []
                for solvent in well_details['solvents']:
                    if solvent[0].lower() in ['dmf', 'dmso', 'water']:
                        original_well_solvents.append(solvent)
                needed_volume = additional_prep_details['initial_solvent_volume'] - sum([item[1] for item in original_well_solvents])
                if needed_volume > 0:
                    original_well_solvents.append([initial_filter_solvent, needed_volume])
                filter_fraction = additional_prep_details['mca_filter_volume'] / additional_prep_details['initial_solvent_volume']
                original_well_solvents = [[solvent[0], solvent[1] * filter_fraction] for solvent in original_well_solvents]
                transfer_fraction = volume_needed / additional_prep_details['mca_filter_volume']
                transferred_solvents = [[solvent[0], solvent[1] * transfer_fraction] for solvent in original_well_solvents]
                if extra_volume_needed > 0:
                    transferred_solvents.append([hplc_container_solvent, extra_volume_needed])

                remaining_hplc_volume = 0
                for solvent in transferred_solvents:
                    if solvent[0].lower() in ['dmf', 'dmso', 'water']:
                        remaining_hplc_volume += solvent[1]
                dilution_factor = (extra_volume_needed + remaining_hplc_volume) / remaining_hplc_volume

                hplc_contents_formatted[well_key] = {'plate_well': well_key,
                                                     'reagents': contents_to_prepare[well_key]['reagents'],
                                                     'solvents': transferred_solvents,
                                                     'total_volume': sum([item[1] for item in transferred_solvents]),
                                                     'target_product': contents_to_prepare[well_key]['target_product'],
                                                     'dilution_factor': dilution_factor,
                                                     'dilution_volume': extra_volume_needed,
                                                     'fraction_transferred': transfer_fraction,
                                                     'confirmation': 'none',
                                                     'date_updated': current_time.strftime('%m/%d/%Y %H:%M:%S')}
                for key in contents_to_prepare[well_key].keys():
                    if key not in hplc_contents_formatted[well_key].keys():
                        hplc_contents_formatted[well_key][key] = contents_to_prepare[well_key][key]

            # Also we need to update the filtrate plate as well
            initial_filter_solvent = additional_prep_details['initial_filter_solvent']
            original_well_solvents = []
            for solvent in well_details['solvents']:
                if solvent[0].lower() in ['dmf', 'dmso', 'water']:
                    original_well_solvents.append(solvent)
            needed_volume = additional_prep_details['initial_solvent_volume'] - sum([item[1] for item in original_well_solvents])
            if needed_volume > 0:
                original_well_solvents.append([initial_filter_solvent, needed_volume])
            filter_fraction = additional_prep_details['mca_filter_volume'] / additional_prep_details['initial_solvent_volume']
            original_well_solvents = [[solvent[0], solvent[1] * filter_fraction] for solvent in original_well_solvents]
            final_contents_formatted[well_key] = {'plate_well': well_key,
                                                  'reagents': contents_to_prepare[well_key]['reagents'],
                                                  'solvents': original_well_solvents,
                                                  'total_volume': sum([item[1] for item in original_well_solvents]),
                                                  'target_product': contents_to_prepare[well_key]['target_product'],
                                                  'confirmation': 'none',
                                                  'date_updated': current_time.strftime('%m/%d/%Y %H:%M:%S')}
            if hplc_contents_formatted and well_key in final_contents_formatted:
                updated_solvents = []
                for solvent in final_contents_formatted[well_key]['solvents']:
                    updated_solvents.append([solvent[0], solvent[1] * (1 - hplc_contents_formatted[well_key]['fraction_transferred'])])
                final_contents_formatted[well_key]['solvents'] = updated_solvents
                total_volume = sum([item[1] for item in updated_solvents])
                final_contents_formatted[well_key]['total_volume'] = total_volume

            for key in contents_to_prepare[well_key].keys():
                if key not in final_contents_formatted[well_key].keys():
                    final_contents_formatted[well_key][key] = contents_to_prepare[well_key][key]

        db_initial_plate = full_library_mongo['wellplates'][str(initial_plate_details['_id'])]
        # if initial_plate_details['prep_location'] != initial_plate_details['location'][1]:
        #    db_initial_plate['location'][1] = initial_plate_details['prep_location']
        db_initial_plate['date_updated'] = current_time.strftime('%m/%d/%Y')

        db_final_plate = full_library_mongo['wellplates'][str(final_plate_details['_id'])]
        db_final_plate['location'][1] = final_plate_details['location'][1]
        db_final_plate['contents'] = final_contents_formatted
        db_final_plate['container_category'] = additional_prep_details['target_container']
        db_final_plate['date_updated'] = current_time.strftime('%m/%d/%Y')

        db_filter_plate = full_library_mongo['consumables'][str(filter_plate_details['_id'])]
        db_filter_plate['contents'] = final_contents_formatted
        db_filter_plate['container_category'] = 'solid_filter_plate'
        db_filter_plate['date_updated'] = current_time.strftime('%m/%d/%Y')

        if hplc_contents_formatted:
            db_hplc_plate = full_library_mongo['wellplates'][str(hplc_wellplate_details['_id'])]
            db_hplc_plate['location'][1] = hplc_wellplate_details['location'][1]
            db_hplc_plate['contents'] = hplc_contents_formatted
            db_hplc_plate['container_category'] = additional_prep_details['hplc_container']
            db_hplc_plate['date_updated'] = current_time.strftime('%m/%d/%Y')

        if __name__ == "__main__":
            pass
            return ['Testing!!!!', success_statements]

        platform_mongodb = pymongo.MongoClient(additional_details['platform_mongodb_address'],
                                               additional_details['platform_mongodb_port'])
        for library_key in full_library_mongo.keys():
            previous = initial_full_library[library_key]
            updated = full_library_mongo[library_key]
            return_statements = update_database(previous, updated, library_key, platform_mongodb)

        updated_queue_document = copy.deepcopy(queue_document)
        if __name__ == "__main__":
            updated_queue_document['operations'][operation_number]['completed'] = 'yes'

        updated_queue_document['containers'][additional_prep_details['target_container']]['container_name'] = \
            final_plate_details['container_name']
        updated_queue_document['containers'][additional_prep_details['target_container']]['contents'] = final_contents_formatted
        updated_queue_document['containers']['solid_filter_plate']['container_name'] = filter_plate_details['container_name']
        updated_queue_document['containers']['solid_filter_plate']['contents'] = final_contents_formatted
        if hplc_contents_formatted:
            hplc_container_name = additional_prep_details['hplc_container']
            updated_queue_document['containers'][hplc_container_name]['container_name'] = hplc_wellplate_details['container_name']
            updated_queue_document['containers'][hplc_container_name]['contents'] = hplc_contents_formatted
        flag_return, return_statements = update_database_document(queue_document, updated_queue_document, 'queue',
                                                                  platform_mongodb)
        filter_plate_details_id = filter_plate_details['_id']
        flag_return, return_statements = move_database_document(filter_plate_details_id, 'consumables', 'wellplates',
                                                                platform_mongodb)
        platform_mongodb.close()
        return ['Success', additional_details['schedule_path_output'], success_statements]
    except Exception:
        return ['Error', traceback.format_exc()]


if __name__ == "__main__":
    queuename = 'tztz_scaffold_20230517_1_20230524'
    prep_return = filter_plate_sequence(queuename, '10')
    pprint(prep_return, width = 250)
