# -*- coding: utf-8 -*-
"""
Updated on 4/09/2022 --- MongoDB containing version of plate preparation

@author: Brent
"""

import copy
import pymongo
import traceback
from pprint import pprint
from copy import deepcopy
import datetime
import math
import sys
import natsort

if __name__ == "__main__":
    from supporting_functions import get_queue_information, convert_to_method_strings
    from platform_library_functions import update_database, update_database_document
    from useful_other_functions import liha_grouping, find_accessible_open_location, \
        check_for_problems, solvent_finder_origin, priority_sort, multi_transfer_grouping, \
        find_open_heater_shaker, previous_product_finder, find_labware, reagent_finder
    from worktable_cleaner import initial_worktable_prep, worktable_export
else:
    from Evoware_API.supporting_functions import get_queue_information, convert_to_method_strings
    from Evoware_API.platform_library_functions import update_database, update_database_document
    from Evoware_API.useful_other_functions import liha_grouping, find_accessible_open_location, \
        check_for_problems, solvent_finder_origin, priority_sort, multi_transfer_grouping, \
        find_open_heater_shaker, previous_product_finder, find_labware, reagent_finder
    from Evoware_API.worktable_cleaner import initial_worktable_prep, worktable_export

# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# These locations need to be updated with adjustments to the worktable
platform_locations = {'low_temperature_prep_locations': [[11, 1, 322], [11, 2, 322]],
                      'thermoshake_locations': [[11, 1, 322], [11, 2, 322], [61, 1, 324], [61, 2, 324]],
                      'teleshake_locations': [[25, 1, 323], [25, 2, 323]],

                      'standard_transfer_locations': [[28, 1, 329]],
                      'special_transfer_locations': [[22, 5, 84], [22, 6, 84], [28, 1, 329]],

                      'inert_box_block_hotel': [6, 1, 335],
                      'inert_box_shaker_location': [25, 4, 323],
                      'small_inert_box_block_hotel': [7, 1, 344]}

# Need to add the small inert_box location and hotel to this list


# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------


def reaction_plate_preparation(queue_name, operation_number):
    start_time = datetime.datetime.now()
    success_statements = []
    try:
        # ----------------------------------------------------------------------------
        # First we start with opening the queue document and getting operation info.
        queue_information = get_queue_information(queue_name, operation_number, 'prepare_wellplate')
        if queue_information[0] != 'Success':
            return ['Error', 'Problem: %s' % queue_information[1], success_statements]
        else:
            return_statement, full_library_mongo, queue_document, additional_details, statements = queue_information

        success_statements.extend(statements)
        initial_full_library = copy.deepcopy(full_library_mongo)
        additional_prep_details = additional_details['additional_prep_details']

        # Then build the current worktable dictionary for building output files, we are
        # also making a duplicate of the worktable, so it can be updated in this code
        return_statement, start_carriers, carriers_labware = initial_worktable_prep(
            additional_details['carriers_labware_filepath'], additional_details['empty_worktable_filepath'],
            full_library_mongo)

        if return_statement != 'Success':
            return ['Error', 'Problem: %s' % start_carriers, success_statements]
        updated_carriers = copy.deepcopy(start_carriers)
        success_statements.append('Loaded the carrier and labware library file!')
        # ----------------------------------------------------------------------------
        
        # Start by opening the reaction wellplate that we want to prepare
        reaction_plate_to_prepare = copy.deepcopy(queue_document['containers'][
                                                      additional_details['target_container']]['contents'])

        # We need to verify that all the chemicals are present and get their information
        problematic_wells = []
        for well_location in reaction_plate_to_prepare.keys():
            # Starting with reagents on the platform, we look for previous products (for multistep reactions)
            reagents = reaction_plate_to_prepare[well_location]['reagents']
            if 'final_product' in reaction_plate_to_prepare[well_location].keys():
                # TODO What to do when a previous product is needed from a different reaction tree?
                final_product_smiles = reaction_plate_to_prepare[well_location]['final_product'][2]
            else:
                final_product_smiles = ''
            return_statement, well_previous_products = previous_product_finder(full_library_mongo,
                                                                               final_product_smiles,
                                                                               reagents)
            if return_statement != 'Success':
                problematic_wells.append(well_previous_products)
            if well_previous_products:
                reaction_plate_to_prepare[well_location]['previous_product_info'] = deepcopy(well_previous_products)

            # Then check on the reagent stock solutions present
            return_statement, well_reagents = reagent_finder(full_library_mongo, reagents, {})
            if return_statement != 'Success':
                problematic_wells.append(well_reagents)
            reaction_plate_to_prepare[well_location]['reagents_info'] = deepcopy(well_reagents)

            # Finally check on the solvents present
            solvents = reaction_plate_to_prepare[well_location]['solvents']
            return_statement, well_solvents = solvent_finder_origin(full_library_mongo, solvents, {})
            if return_statement != 'Success':
                problematic_wells.append(well_solvents)
            reaction_plate_to_prepare[well_location]['solvents_info'] = deepcopy(well_solvents)
            reaction_plate_to_prepare[well_location]['remaining_volume'] = reaction_plate_to_prepare[
                well_location]['total_volume']
            
            # Now check to see if anything is missing, this happens when a chemical is not able to be found in the
            # library. Common causes are issues that arise when setting up final product trees
            if reagents:
                reagents_found = set()
                if type(well_previous_products) == dict:
                    reagents_found.update([item.split('>')[0] for item in well_previous_products.keys()])
                    # reagents_found.update(list(well_previous_products.keys()))
                if type(well_reagents) == dict:
                    reagents_found.update(list(well_reagents.keys()))
                missing_reagents = set([reagent[0] for reagent in reagents]).difference(reagents_found)
                if missing_reagents:
                    problematic_wells.append('Missing reagents %s in well %s' % (list(missing_reagents), well_location))
            solvents_found = set()
            if type(well_solvents) == dict:
                solvents_found.update(list(well_solvents.keys()))
            missing_solvents = set([solvent[0] for solvent in solvents]).difference(solvents_found)
            if missing_solvents:
                problematic_wells.append('Missing solvents %s in well %s' % (list(missing_solvents), well_location))
        if problematic_wells:
            return ['Error', 'Problem: %s' % problematic_wells, success_statements]

        # Once we know where all the chemicals that we need are located, we can
        # figure out how much to add and when to add the reagents
        reagents_needed = {}
        previous_products_needed = {}
        for well_location in reaction_plate_to_prepare.keys():
            reaction_scale = 1
            if 'previous_product_info' in reaction_plate_to_prepare[well_location].keys():
                # reaction_scaling = []
                for previous_product_key in reaction_plate_to_prepare[well_location]['previous_product_info'].keys():
                    product_info = deepcopy(reaction_plate_to_prepare[well_location]['previous_product_info'][
                                                previous_product_key])
                    initial_reagent_details = deepcopy(full_library_mongo[product_info['collection']][
                                                           product_info['container_id']]['contents'][
                                                           product_info['origin_well_location']])
                    if previous_product_key not in previous_products_needed.keys():
                        previous_products_needed[previous_product_key] = initial_reagent_details
                        reaction_scale = product_info['reaction_scale']
                        reaction_scale = 1
                        previous_products_needed[previous_product_key].update(
                            {'collection': product_info['collection'],
                             'amount_needed': product_info['amount_needed'] * reaction_scale,
                             'container_id': product_info['container_id'],
                             'container_name': product_info['container_name'],
                             'origin_well_location': product_info['origin_well_location'],
                             'transfer_sequence': product_info['transfer_sequence'],
                             'destinations': {}})
                    previous_sequences = previous_products_needed[previous_product_key]['destinations']
                    if product_info['transfer_sequence'] not in previous_sequences.keys():
                        previous_sequences[product_info['transfer_sequence']] = []
                    previous_sequences[product_info['transfer_sequence']].append([well_location,
                                                                                  product_info['amount_needed'] *
                                                                                  reaction_scale])
                    previous_volume = 0
                    previous_solvents = initial_reagent_details['solvents']
                    for previous_solvent in previous_solvents:
                        if previous_solvent[0] in ['dmf', 'dmso']:
                            previous_volume += previous_solvent[1]
                    reaction_plate_to_prepare[well_location]['remaining_volume'] -= previous_volume

            if 'reagents_info' in reaction_plate_to_prepare[well_location].keys():
                for reagent_key in reaction_plate_to_prepare[well_location]['reagents_info'].keys():
                    reagent_info = reaction_plate_to_prepare[well_location]['reagents_info'][reagent_key]
                    if reagent_key not in reagents_needed.keys():
                        initial_reagent_details = deepcopy(full_library_mongo[reagent_info['collection']][
                                                               reagent_info['container_id']]['contents'][
                                                               reagent_info['origin_well_location']])
                        reagents_needed[reagent_key] = initial_reagent_details
                        reagents_needed[reagent_key].update({'volume_needed': 0, 'destinations': {},
                                                             'container_id': reagent_info['container_id'],
                                                             'collection': reagent_info['collection']})
                    reagent_sequences = reagents_needed[reagent_key]['destinations']
                    t_seq = reagent_info['transfer_sequence']
                    if t_seq not in reagent_sequences.keys():
                        reagent_sequences[t_seq] = []
                    reagent_sequences[t_seq].append([well_location, reagent_info['amount_needed'] * reaction_scale,
                                                     reagent_info['volume_needed'] * reaction_scale])
                    reagents_needed[reagent_key]['volume_needed'] += reagent_info['volume_needed']
                    reaction_plate_to_prepare[well_location]['remaining_volume'] -= reagent_info['volume_needed']
        
        solvents_needed = {}
        for well_location in reaction_plate_to_prepare.keys():
            solvent_list = reaction_plate_to_prepare[well_location]['solvents']
            solvent_norm = sum([solvent[1] for solvent in solvent_list])
            for solvent in solvent_list:
                volume_needed = reaction_plate_to_prepare[well_location]['remaining_volume'] * (
                            solvent[1] / solvent_norm)
                if solvent[0] not in solvents_needed.keys():
                    solvents_needed[solvent[0]] = deepcopy(reaction_plate_to_prepare[well_location][
                                                               'solvents_info'][solvent[0]])
                    solvents_needed[solvent[0]].update({'destinations': {}, 'volume_needed': 0})
                if len(solvent) == 2:
                    solvent_sequence = 1
                else:
                    solvent_sequence = 3
                if solvent_sequence not in solvents_needed[solvent[0]]['destinations'].keys():
                    solvents_needed[solvent[0]]['destinations'][solvent_sequence] = []
                solvents_needed[solvent[0]]['destinations'][solvent_sequence].append([well_location,
                                                                                      volume_needed,
                                                                                      solvent_sequence])
                solvents_needed[solvent[0]]['volume_needed'] += volume_needed
                reaction_plate_to_prepare[well_location]['solvents_info'][solvent[0]]['volume_needed'] += volume_needed

        # Check the chemicals to see if there is enough on the platform
        problematic_chemicals = {'reagents': [], 'solvents': []}
        safety_factor = 1.1
        for solvent_key in solvents_needed.keys():
            required_volume = solvents_needed[solvent_key]['volume_needed'] * safety_factor
            original_volume = solvents_needed[solvent_key]['solvent_volume']
            if required_volume > original_volume:
                problematic_chemicals['solvents'].append('Not enough %s present (%s of %s)' % (solvent_key,
                                                                                               required_volume,
                                                                                               original_volume))
            
        for reagent_key in reagents_needed.keys():
            original_volume = reagents_needed[reagent_key]['volume_ul']
            required_volume = reagents_needed[reagent_key]['volume_needed'] * safety_factor
            if required_volume > original_volume:
                problematic_chemicals['reagents'].append('Not enough %s material present (%s v. %s)' % (reagent_key,
                                                                                                       required_volume,
                                                                                                       original_volume))
            elif original_volume - required_volume < 100:
                problematic_chemicals['reagents'].append('Not enough volume of %s present, final volume of %s' % (reagent_key,
                                                          original_volume-required_volume))
        if sum([len(problematic_chemicals[type_key]) for type_key in problematic_chemicals.keys()]) != 0:
            return ['Error', 'Problem: %s' % problematic_chemicals, success_statements]

        # With the solvents and reagents we can group them by originating location
        return_statement, previous_products_order = priority_sort(previous_products_needed, 'reagents')
        if return_statement == 'Error':
            return ['Error', 'Problem: %s' % previous_products_order, success_statements]
        return_statement, reagent_order = priority_sort(reagents_needed, 'reagents')
        if return_statement == 'Error':
            return ['Error', 'Problem: %s' % reagent_order, success_statements]
        return_statement, solvent_order = priority_sort(solvents_needed, 'solvents')
        if return_statement == 'Error':
            return ['Error', 'Problem: %s' % reagent_order, success_statements]
        success_statements.append('Found all of the reagents and solvents on the platform')
        
        # This gives us all the chemicals that we need to distribute into wells, so now
        # we can build up the schedule of tasks needed to distribute the chemicals, we
        # need to get the locations that we are placing the plates at and also find an
        # empty reaction plate that we can use for the method
        simplified_schedule, method_strings = [], []
        return_statement, reaction_plate_initial = find_labware(updated_carriers, full_library_mongo,
                                                                additional_details['wellplate_name'],
                                                                additional_details['target_well_plate_type'],
                                                                additional_details['target_container'])
        if return_statement != 'Success':
            return ['Error', 'Problem: %s' % reaction_plate_initial, success_statements]
        success_statements.append('Found an acceptable reaction plate %s at %s' %
                                  (reaction_plate_initial['container_name'], reaction_plate_initial['location'][1]))
        additional_prep_details = additional_details['additional_prep_details']

        # Check for special preparation routines
        if 'inert_atmosphere' in additional_prep_details.keys() and additional_prep_details['inert_atmosphere'] is True:
            reaction_well_plate_home = [25, 2, 323]
            if updated_carriers[reaction_well_plate_home[0]][reaction_well_plate_home[2]][
                'labware_types'][reaction_well_plate_home[1] - 1] != '':
                return ['Error', 'Problem: The inert reaction heater shaker is occupied']

        elif 'low_temperature' in additional_prep_details.keys() and additional_prep_details['low_temperature'] is True:
            labware_type = additional_details['target_well_plate_type']
            acceptable_locations = [[11, 1, 322], [11, 2, 322]]
            return_statement, reaction_well_plate_home = find_open_heater_shaker(updated_carriers, 5, labware_type,
                                                                                 platform_locations)
            if return_statement == 'Error':
                return ['Error', 'Problem: %s' % reaction_well_plate_home, success_statements]
            if reaction_well_plate_home not in acceptable_locations:
                return ['Error', 'Problem: %s' % reaction_well_plate_home, success_statements]

        else:
            allowed_prep_carriers = [32, 38, 44]
            if reaction_plate_initial['location'][1][0] in allowed_prep_carriers:
                reaction_well_plate_home = reaction_plate_initial['location'][1]
            else:
                return_statement, reaction_well_plate_home = find_accessible_open_location(updated_carriers)
                if return_statement != 'Success':
                    return ['Error', 'Problem: %s' % reaction_well_plate_home, success_statements]

        # Now double check for any problems
        problem_return = check_for_problems(updated_carriers, reaction_well_plate_home,
                                            additional_details['target_well_plate_type'])
        if problem_return != 'Acceptable':
            return ['Error', 'Problem: Site %s is not accessible' % reaction_well_plate_home, success_statements]
        success_statements.append('Found a suitable prep location at location: %s' % reaction_well_plate_home)

        if reaction_plate_initial['location'][1] != reaction_well_plate_home:
            simplified_schedule.append(['TRANSFER_LABWARE', reaction_plate_initial['container_name'],
                                        reaction_plate_initial['location'][1], reaction_well_plate_home,
                                        str(reaction_plate_initial['_id']), 'wellplates'])
            updated_carriers[reaction_well_plate_home[0]][reaction_well_plate_home[2]]['labware_labels'][
                reaction_well_plate_home[1] - 1] = 'location_%s_%s_%s' % (str(reaction_well_plate_home[0]),
                                                                          str(reaction_well_plate_home[1]),
                                                                          str(reaction_well_plate_home[2]))
            updated_carriers[reaction_well_plate_home[0]][reaction_well_plate_home[2]]['labware_types'][
                reaction_well_plate_home[1] - 1] = reaction_plate_initial['labware_type']

        # For the special preparation conditions there are additional methods needed
        if 'inert_atmosphere' in additional_prep_details.keys() and additional_prep_details['inert_atmosphere'] is True:
            if additional_details['target_well_plate_type'] in ['96 Well Microplate', '96 Well MediumWell']:
               inert_box_type = 'Small Inert Box Block'
               inert_box_block_hotel = platform_locations['small_inert_box_block_hotel']
            else:
               inert_box_type = 'Inert Box Block'
               inert_box_block_hotel = platform_locations['inert_box_block_hotel']
            for consumable_id in full_library_mongo['consumables'].keys():
                if 'labware_type' not in full_library_mongo['consumables'][consumable_id].keys():
                    continue
                if full_library_mongo['consumables'][consumable_id]['location'][0] != 'liquid_handler':
                    continue
                # TODO: Check to see what type of labware is coming to the heater-shaker
                if full_library_mongo['consumables'][consumable_id]['labware_type'] == inert_box_type:
                    relevant_document_id = consumable_id
                    break
                #if full_library_mongo['consumables'][consumable_id]['labware_type'] == 'Inert Box Block':
            else:
                return ['Error', 'Problem: Inert reaction block not found in consumables', success_statements]

            inert_box_shaker_location = platform_locations['inert_box_shaker_location']
            if full_library_mongo['consumables'][relevant_document_id]['location'][1] != inert_box_block_hotel:
                return ['Error', 'Problem: Inert reaction block not in %s, it is at %s' % (inert_box_block_hotel,
                                           full_library_mongo['consumables'][relevant_document_id]['location'][1]),
                                           success_statements]
            simplified_schedule.append(['SPECIAL_TRANSFER_LABWARE', 'Inert Box Block',
                                        inert_box_block_hotel, inert_box_shaker_location, '', ''])
            full_library_mongo['consumables'][relevant_document_id]['location'][1] = inert_box_shaker_location

        if ('inert_atmosphere' in additional_prep_details.keys() and additional_prep_details[
            'inert_atmosphere'] is True) or \
                ('low_temperature' in additional_prep_details.keys() and additional_prep_details[
                    'low_temperature'] is True):
            if 'temperature_setpoint' in additional_prep_details.keys():
                try:
                    simplified_schedule.append(['SHAKER_TEMPERATURE', reaction_well_plate_home, 'on',
                                                float(additional_prep_details['temperature_setpoint'])])
                except Exception:
                    simplified_schedule.append(['SHAKER_TEMPERATURE', reaction_well_plate_home, 'on', 20])
            if 'preparation_shaking' in additional_prep_details.keys() and additional_prep_details[
                'preparation_shaking'] is True:
                pass  # This is currently not enabled because of paradox plate problems
                # try:
                #    simplified_schedule.append(['START_SHAKER', reaction_well_plate_home,
                #                                int(additional_prep_details['shaking_rpms'])])
                # except:
                #    simplified_schedule.append(['START_SHAKER', reaction_well_plate_home, 500])
            if 'waiting_prep_time' in additional_prep_details.keys():
                try:
                    simplified_schedule.append(['START_TIMER', 1])
                    simplified_schedule.append(
                        ['WAIT_FOR_TIMER', 1, float(additional_prep_details['waiting_prep_time'])])
                except:
                    pass
            if 'inert_purge_time' in additional_prep_details.keys():
                try:
                    if int(additional_prep_details['inert_purge_time']) != 0:
                        simplified_schedule.append(['START_TIMER', 1])
                        simplified_schedule.append(
                            ['WAIT_FOR_TIMER', 1, int(additional_prep_details['inert_purge_time'])])
                except:
                    simplified_schedule.append(['START_TIMER', 1])
                    simplified_schedule.append(['WAIT_FOR_TIMER', 1, 5 * 60])
                    success_statements.append(
                        'Added 300s inert purge time but %s is not parsable' % additional_prep_details[
                            'inert_purge_time'])
            success_statements.append('Added the special prep method strings to schedule')

        # Rinse the LiHa tips since they will be used to prepare the well-plate
        if additional_prep_details['tip_prep_rinse']:
            simplified_schedule.append(['TIP_RINSE', 10, 10])
        else:
            simplified_schedule.append(['TIP_RINSE', 5, 5])
        simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 150])
        simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 150])
        if 'washing_frequency' in additional_prep_details.keys() and additional_prep_details['washing_frequency'] is not None:
            if additional_prep_details['washing_frequency'] == 'nonpolar_stepdown':
                if 'nonpolar_solvents' in additional_prep_details.keys():
                    if type(additional_prep_details['nonpolar_solvents']) == list:
                        nonpolar_solvents = additional_prep_details['nonpolar_solvents']
                    elif type(additional_prep_details['nonpolar_solvents']) == str:
                        nonpolar_solvents = [additional_prep_details['nonpolar_solvents']]
                    else:
                        nonpolar_solvents = []
                    for nonpolar_solvent in nonpolar_solvents:
                        simplified_schedule.append(['SOLVENT_RINSE', nonpolar_solvent, 150])

        # Now the aspirate and dispense functions are added to the simplified schedule
        relevant_collections = ['solvents', 'reagents']
        for solvent_sequence in solvent_order.keys():
            solvent_index = 0
            for target_container in solvent_order[solvent_sequence].keys():
                for relevant_collection in relevant_collections:
                    if target_container in full_library_mongo[relevant_collection].keys():
                        break
                else:
                    return ['Error', 'Problem: Unable to find solvent plate: %s' % target_container, success_statements]
                solvent_definition = full_library_mongo[relevant_collection][target_container]
                current_plate_location = solvent_definition['location'][1]
                if current_plate_location[0] == 44 and current_plate_location[1] == 3 and solvent_index == 0:
                    simplified_schedule.append(['MOVE_LiHa'])
                    simplified_schedule.append(['ROMA_MOVE', 'Tip Rack Waste Chute_Wide_1', [62, 1, 331], 0, {'back': 1}])
                for target_solvent in solvent_order[solvent_sequence][target_container].keys():
                    if solvent_definition['labware_type'] in ['8 Position Vial Carrier', '12 Position Vial Carrier',
                                                              '12 Position Inert Box']:
                        solvent_origin_well = solvents_needed[target_solvent]['origin_well_location']
                        
                        for target_destination in solvent_order[solvent_sequence][target_container][target_solvent][
                            'destinations']:
                            if 'liquid_class' in additional_prep_details:
                                simplified_schedule.append(
                                    ['LiHa_ASPIRATE_DISPENSE', solvent_definition['location'][1] + [[solvent_origin_well]],
                                     reaction_well_plate_home, [1], [target_destination[1]], [target_destination[0]],
                                    additional_prep_details['liquid_class']])
                            else:
                                simplified_schedule.append(
                                    ['LiHa_ASPIRATE_DISPENSE', solvent_definition['location'][1] + [[solvent_origin_well]],
                                     reaction_well_plate_home, [1], [target_destination[1]], [target_destination[0]]])
                    else:
                        destinations = solvent_order[solvent_sequence][target_container][target_solvent]['destinations']
                        return_statement, liha_groups = liha_grouping(destinations, carriers_labware,
                                                                      solvent_definition['labware_type'],
                                                                      reaction_plate_initial['labware_type'])
                        if return_statement == 'Error':
                            return ['Error', 'Problem: Something went wrong with LiHa grouping', success_statements]
                        for group_num, liha_group in enumerate(liha_groups.keys()):
                            group = liha_groups[liha_group]
                            wells, volumes, tips = map(list, zip(*group))
                            simplified_schedule.append(['LiHa_ASPIRATE_DISPENSE', solvent_definition['location'][1],
                                                        reaction_well_plate_home, tips, volumes, wells])
                            
                            
                            if 'washing_frequency' in additional_prep_details.keys() and additional_prep_details[
                                'washing_frequency'] is not None:
                                simplified_schedule.append(['MOVE_LiHa'])
                                simplified_schedule.append(['TIP_RINSE', 3, 3])
                                simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 150])
                                simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 150])
                                
                                if additional_prep_details['washing_frequency'] == 'nonpolar_stepdown':
                                    if 'nonpolar_solvents' in additional_prep_details.keys():
                                        if type(additional_prep_details['nonpolar_solvents']) == list:
                                            nonpolar_solvents = additional_prep_details['nonpolar_solvents']
                                        elif type(additional_prep_details['nonpolar_solvents']) == str:
                                            nonpolar_solvents = [additional_prep_details['nonpolar_solvents']]
                                        else:
                                            nonpolar_solvents = []
                                        for nonpolar_solvent in nonpolar_solvents:
                                            simplified_schedule.append(['SOLVENT_RINSE', nonpolar_solvent, 150])
                    solvent_index += 1
                    simplified_schedule.append(['MOVE_LiHa'])
                    simplified_schedule.append(['TIP_RINSE', 5, 5])
                    simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 150])
                    simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 150])
                    if 'washing_frequency' in additional_prep_details.keys() and additional_prep_details[
                        'washing_frequency'] == 'nonpolar_stepdown':
                        if type(additional_prep_details['nonpolar_solvents']) == list:
                            nonpolar_solvents = additional_prep_details['nonpolar_solvents']
                        elif type(additional_prep_details['nonpolar_solvents']) == str:
                            nonpolar_solvents = [additional_prep_details['nonpolar_solvents']]
                        else:
                            nonpolar_solvents = []
                        for nonpolar_solvent in nonpolar_solvents:
                            simplified_schedule.append(['SOLVENT_RINSE', nonpolar_solvent, 150])

                    
        success_statements.append('Added solvents to the wellplate being prepared')
        
        # Check if there was an optional statement to wait after adding solvents
        if 'solvent_waiting_time' in additional_prep_details.keys():
            try:
                simplified_schedule.append(['START_TIMER', 1])
                simplified_schedule.append(
                    ['WAIT_FOR_TIMER', 1, float(additional_prep_details['solvent_waiting_time'])])
                success_statements.append(
                    'Solvent waiting time of %s added' % additional_prep_details['solvent_waiting_time'])
            except:
                success_statements.append(
                    'Skipping waiting time, problem with: %s' % additional_prep_details['solvent_waiting_time'])

        # We need to keep track of the well-plate type moved into each bed position
        # otherwise there are problems with Evoware if two different labware types
        # are moved to the same position during the script
        location_occupancies = []
        liha_tip_number = 1
        acceptable_bed_locations = [32, 38, 44]
        
        # TODO Check on how previous products are being brough in now...
        previous_products_solvents = {}
        for previous_products_sequence in previous_products_order.keys():
            for target_container in previous_products_order[previous_products_sequence].keys():
                previous_product_definition = deepcopy(full_library_mongo['wellplates'][target_container])
                current_plate_location = previous_product_definition['location'][1]
                success_statements.append('Currently a previous product plate "%s" is located at location: %s' % (
                    previous_product_definition['container_name'], current_plate_location))
                if current_plate_location[0] not in acceptable_bed_locations:
                    exceptions = [location[0] for location in location_occupancies if
                                  location[1] != previous_product_definition['labware_type']]
                    return_statement, previous_product_plate_home = find_accessible_open_location(updated_carriers,
                                                                                                  exceptions,
                                                                                                  previous_product_definition[
                                                                                                      'labware_type'])
                    if return_statement != 'Success':
                        return ['Error', 'Problem: No positions available for reagent tray of type %s' %
                                previous_product_definition, success_statements]
                    if current_plate_location != previous_product_plate_home:
                        simplified_schedule.append(['TRANSFER_LABWARE', previous_product_definition['container_name'],
                                                    previous_product_definition['location'][1],
                                                    previous_product_plate_home,
                                                    str(previous_product_definition['_id']), 'wellplates'])
                        current_plate_location = previous_product_plate_home
                location_occupancies.append([current_plate_location, previous_product_definition['labware_type']])
                for target_product in previous_products_order[previous_products_sequence][target_container].keys():
                    product_origin_well = previous_products_needed[target_product]['plate_well']
                    if len(previous_products_order[previous_products_sequence][target_container][target_product][
                               'destinations']) > 1:
                        return ['Error', 'Problem: Previous product splitting into multiple wells not allowed',
                                success_statements]
                    for target_destination in \
                    previous_products_order[previous_products_sequence][target_container][target_product][
                        'destinations']:
                        previous_definition = previous_products_needed[target_product]
                        previous_volume = 0
                        for previous_solvent in previous_definition['solvents']:
                            if previous_solvent[0] in ['dmf', 'dmso']:
                                previous_volume += previous_solvent[1]
                        transfer_volume = 200
                        transfer_solvent = 'chloroform'
                        if previous_volume < transfer_volume:
                            additional_volume = transfer_volume - previous_volume
                            solvent_dict_stub = [[transfer_solvent,
                                                  additional_volume]]  # {'destination': [[target_destination[0], additional_volume, 1]]}
                            return_statements, solvent_information = solvent_finder_origin(full_library_mongo,
                                                                                           solvent_dict_stub, '')
                            if return_statements == 'Error':
                                return ['Error', 'Problem: %s' % solvent_information, success_statements]
                            solvent_definition = \
                            full_library_mongo[solvent_information[transfer_solvent]['collection']][
                                solvent_information[transfer_solvent]['container_id']]
                            solvent_origin_well = solvent_information[transfer_solvent]['origin_well_location']
                            simplified_schedule.append(
                                ['LiHa_ASPIRATE_DISPENSE', solvent_definition['location'][1] + [[solvent_origin_well]],
                                 current_plate_location, [liha_tip_number], [additional_volume], [product_origin_well]])
                            full_library_mongo[solvent_information[transfer_solvent]['collection']][
                                solvent_information[transfer_solvent]['container_id']]['contents'][solvent_origin_well][
                                'volume_ul'] -= additional_volume
                            if solvent_information[transfer_solvent]['solvent_volume'] - additional_volume < 5000:
                                return ['Error',
                                        'Problem: Not enough "%s" solvent present on the liquid handler' % transfer_solvent,
                                        success_statements]
                            if transfer_solvent not in previous_products_solvents.keys():
                                previous_products_solvents[transfer_solvent] = {'destinations': []}
                            previous_products_solvents[transfer_solvent]['destinations'].append(target_destination)
                            simplified_schedule.append(['MIX_LiHa', current_plate_location + [[product_origin_well]],
                                                        [liha_tip_number], [product_origin_well], [transfer_volume]])

                        product_transfer_volume = 450
                        if product_transfer_volume > 300:
                            number_of_transfers = int(product_transfer_volume // 200)
                        else:
                            number_of_transfers = 1
                        for transfer in range(0, number_of_transfers):
                            simplified_schedule.append(
                                ['LiHa_ASPIRATE_DISPENSE', current_plate_location + [[product_origin_well]],
                                 reaction_well_plate_home, [liha_tip_number],
                                 [product_transfer_volume / number_of_transfers], [target_destination[0]]])
                        if liha_tip_number == 8:
                            simplified_schedule.append(['TIP_RINSE', 5, 5])
                            simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 150])
                            simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 150])
                            liha_tip_number = 1
                        else:
                            liha_tip_number += 1
                simplified_schedule.append(['MOVE_LiHa'])
                if previous_product_definition['location'][1] != current_plate_location:
                    simplified_schedule.append(['TRANSFER_LABWARE', previous_product_definition['container_name'],
                                                current_plate_location, previous_product_definition['location'][1],
                                                str(previous_product_definition['_id']), 'wellplates'])
        if previous_products_order:
            simplified_schedule.append(['TIP_RINSE', 5, 5])
            simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 150])
            simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 150])
            simplified_schedule.append(['TIP_RINSE', 4, 4])
        
        reagent_sequence_raw = list(reagent_order.keys())
        reagent_sequence_raw = natsort.natsorted(reagent_sequence_raw)
        for reagent_sequence_index in range(0, len(list(reagent_order.keys()))):
            reagent_sequence = reagent_sequence_raw[reagent_sequence_index]
            
            for target_container in reagent_order[reagent_sequence].keys():
                reagent_tray_definition = deepcopy(full_library_mongo['reagents'][target_container])
                current_plate_location = reagent_tray_definition['location'][1]
                success_statements.append('Currently a reagent tray "%s" is located at location: %s' % (
                    reagent_tray_definition['container_name'], current_plate_location))
                if current_plate_location[0] not in acceptable_bed_locations:
                    exceptions = [location[0] for location in location_occupancies if
                                  location[1] != reagent_tray_definition['labware_type']]
                    return_statement, reagent_well_plate_home = find_accessible_open_location(updated_carriers,
                                                                                              exceptions,
                                                                                              reagent_tray_definition[
                                                                                                  'labware_type'])
                    if return_statement != 'Success':
                        return ['Error', 'Problem: No positions available for reagent tray of type %s' %
                                reagent_tray_definition, success_statements]
                    success_statements.append(
                        'Found a suitable prep location at location: %s' % reagent_well_plate_home)
                    if current_plate_location != reagent_well_plate_home:
                        simplified_schedule.append(['TRANSFER_LABWARE', reagent_tray_definition['container_name'],
                                                    reagent_tray_definition['location'][1], reagent_well_plate_home,
                                                    str(reagent_tray_definition['_id']), 'reagents'])
                        current_plate_location = reagent_well_plate_home
                if current_plate_location[0] == 44 and current_plate_location[1] == 3:
                    simplified_schedule.append(['MOVE_LiHa'])
                    simplified_schedule.append(['ROMA_MOVE', 'Tip Rack Waste Chute_Wide_1', [62, 1, 331], 0,  {'back': 1}])
                location_occupancies.append([current_plate_location, reagent_tray_definition['labware_type']])
                for target_reagent in reagent_order[reagent_sequence][target_container].keys():
                    reagent_origin_well = reagents_needed[target_reagent]['plate_well']
                    if 'single_transfer' in additional_prep_details.keys() and additional_prep_details['single_transfer'] == True:
                        for target_destination in reagent_order[reagent_sequence][target_container][target_reagent]['destinations']:
                            if 'liquid_class' in additional_prep_details:
                                simplified_schedule.append(
                                ['LiHa_ASPIRATE_DISPENSE', current_plate_location + [[reagent_origin_well]],
                                 reaction_well_plate_home, [liha_tip_number], [target_destination[2]],
                                 [target_destination[0]], additional_prep_details['liquid_class']])
                            else:
                                simplified_schedule.append(
                                    ['LiHa_ASPIRATE_DISPENSE', current_plate_location + [[reagent_origin_well]],
                                     reaction_well_plate_home, [liha_tip_number], [target_destination[2]],
                                     [target_destination[0]]])
                    else:
                        max_multi_transfer_volume = 400
                        target_destinations = reagent_order[reagent_sequence][target_container][target_reagent][
                            'destinations']
                        return_statement, reagent_groups = multi_transfer_grouping(target_destinations,
                                                                                   max_multi_transfer_volume)
                        if return_statement != 'Success':
                            return ['Error', 'Problem: %s' % reagent_groups, success_statements]

                        for reagent_group in reagent_groups.keys():
                            destinations = reagent_groups[reagent_group]
                            total_volume = sum([math.floor(10 * dest[2]) / 10 for dest in destinations])
                            simplified_schedule.append(
                                ['LiHa_ASPIRATE', current_plate_location + [[reagent_origin_well]],
                                 [liha_tip_number], [total_volume]])
                            for destination in destinations:
                                simplified_schedule.append(['LiHa_DISPENSE', reaction_well_plate_home,
                                                            [liha_tip_number], [math.floor(10 * destination[2]) / 10],
                                                            [destination[0]]])
                    if liha_tip_number == 8:
                        simplified_schedule.append(['TIP_RINSE', 5, 5])
                        simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 150])
                        simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 150])
                        if 'washing_frequency' in additional_prep_details.keys() and additional_prep_details[
                            'washing_frequency'] == 'nonpolar_stepdown':
                            if type(additional_prep_details['nonpolar_solvents']) == list:
                                nonpolar_solvents = additional_prep_details['nonpolar_solvents']
                            elif type(additional_prep_details['nonpolar_solvents']) == str:
                                nonpolar_solvents = [additional_prep_details['nonpolar_solvents']]
                            else:
                                nonpolar_solvents = []
                            for nonpolar_solvent in nonpolar_solvents:
                                simplified_schedule.append(['SOLVENT_RINSE', nonpolar_solvent, 150])
                        liha_tip_number = 1
                    else:
                        liha_tip_number += 1
                simplified_schedule.append(['MOVE_LiHa'])

                if reagent_tray_definition['location'][1] != current_plate_location:
                    simplified_schedule.append(['TRANSFER_LABWARE', reagent_tray_definition['container_name'],
                                                current_plate_location, reagent_tray_definition['location'][1],
                                                str(reagent_tray_definition['_id']), 'reagents'])
                    success_statements.append('Returned "%s" back to: %s' % (reagent_tray_definition['container_name'],
                                                                             reagent_tray_definition['location'][1]))
            if liha_tip_number != 1:
                simplified_schedule.append(['TIP_RINSE', 5, 5])
                simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 150])
                simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 150])
                if 'washing_frequency' in additional_prep_details.keys() and additional_prep_details[
                    'washing_frequency'] == 'nonpolar_stepdown':
                    if type(additional_prep_details['nonpolar_solvents']) == list:
                        nonpolar_solvents = additional_prep_details['nonpolar_solvents']
                    elif type(additional_prep_details['nonpolar_solvents']) == str:
                        nonpolar_solvents = [additional_prep_details['nonpolar_solvents']]
                    else:
                        nonpolar_solvents = []
                    for nonpolar_solvent in nonpolar_solvents:
                        simplified_schedule.append(['SOLVENT_RINSE', nonpolar_solvent, 150])
                liha_tip_number = 1
            if len(list(reagent_order.keys())) > 1 and 'reagent_sequence_waiting' in additional_prep_details.keys() and \
                    additional_prep_details['reagent_sequence_waiting'] == True:
                if 'reagent_waiting_time' in additional_prep_details.keys():
                    try:
                        reagent_waiting_time = float(additional_prep_details['reagent_waiting_time'])
                    except:
                        reagent_waiting_time = 60
                else:
                    reagent_waiting_time = 60
                simplified_schedule.append(['START_TIMER', 1])
                simplified_schedule.append(['WAIT_FOR_TIMER', 1, reagent_waiting_time])
        success_statements.append('Added all of the reagents to the reaction wellplate')

        # With the reagents added, we need to check if the reaction wellplate needs to be returned
        if ('inert_atmosphere' in additional_prep_details.keys() and additional_prep_details['inert_atmosphere'] == True) or (
                'low_temperature' in additional_prep_details.keys() and additional_prep_details['low_temperature'] == True):
            db_reaction_plate = full_library_mongo['wellplates'][str(reaction_plate_initial['_id'])]
            db_reaction_plate['location'][1] = reaction_well_plate_home
        else:
            if reaction_well_plate_home != reaction_plate_initial['location'][1]:
                simplified_schedule.append(['TRANSFER_LABWARE', reaction_plate_initial['container_name'],
                                            reaction_well_plate_home, reaction_plate_initial['location'][1],
                                            str(reaction_plate_initial['_id']), 'wellplates'])
        
        # Now we need to rebuild the worktable so we can write a file to run
        return_statement, filled_worktable_export = worktable_export(additional_details['carriers_labware_filepath'],
                                                                     additional_details['empty_worktable_filepath'],
                                                                     updated_carriers, full_library_mongo,
                                                                     simplified_schedule)
        if return_statement == 'Error':
            return ['Error', 'Problem: %s' % filled_worktable_export, success_statements]

        # We need to convert the simplified schedule into Tecan method strings
        return_statement, method_strings, updated_carriers = convert_to_method_strings(simplified_schedule,
                                                                                       start_carriers,
                                                                                       full_library_mongo,
                                                                                       carriers_labware)
        if return_statement == 'Error':
            return ['Error', 'Problem: %s' % method_strings, success_statements]

        with open(additional_details['schedule_path_output'], 'w', encoding='latin1') as outfile:
            for line in filled_worktable_export:
                outfile.write(line + '\n')
            for line in method_strings:
                outfile.write(line + '\n')

        # Update the new reaction plate with the contents added to the wellplate
        current_time = datetime.datetime.now()
        for well_key in reaction_plate_to_prepare.keys():
            reaction_plate_to_prepare[well_key]['date_updated'] = current_time.strftime('%m/%d/%Y %H:%M:%S')
            reaction_plate_to_prepare[well_key]['plate_well'] = well_key
            reaction_plate_to_prepare[well_key]['reagents'] = []
            well_solvents = {}
            if 'previous_product_info' in reaction_plate_to_prepare[well_key].keys():
                for previous_product in reaction_plate_to_prepare[well_key]['previous_product_info'].keys():
                    product_details = reaction_plate_to_prepare[well_key]['previous_product_info'][previous_product]
                    for previous_solvent in product_details['previous_solvent']:
                        if previous_solvent[0] not in well_solvents.keys():
                            well_solvents[previous_solvent[0]] = 0
                        well_solvents[previous_solvent[0]] += previous_solvent[1]
                    # TODO update the below when the scale is not 100%
                    reaction_plate_to_prepare[well_key]['reagents'].append([product_details['target_product'],
                                                                            product_details['previous_mols'],
                                                                            product_details['target_product']])
            for reagent in reaction_plate_to_prepare[well_key]['reagents_info'].keys():
                reagent_details = reaction_plate_to_prepare[well_key]['reagents_info'][reagent]
                if 'solvent' not in reagent_details.keys():
                    pprint(reagent_details)
                if reagent_details['solvent'] not in well_solvents.keys():
                    well_solvents[reagent_details['solvent']] = 0
                well_solvents[reagent_details['solvent']] += reagent_details['volume_needed']
                reaction_plate_to_prepare[well_key]['reagents'].append([reagent_details['chemical_name'],
                                                                        reagent_details['amount_needed'],
                                                                        reagent_details['chemical_smiles']])
            reaction_plate_to_prepare[well_key]['solvents'] = []
            for solvent in reaction_plate_to_prepare[well_key]['solvents_info'].keys():
                solvent_details = reaction_plate_to_prepare[well_key]['solvents_info'][solvent]
                if solvent not in well_solvents.keys():
                    well_solvents[solvent] = 0
                well_solvents[solvent] += solvent_details['volume_needed']
            for well_solvent in well_solvents.keys():
                if well_solvent not in solvents_needed.keys():
                    solvent_dict_stub = [[well_solvent, 1]]
                    # Potentially can add this solvent (if not found) to the list of solvents... need to check
                    return_statement, solvent_details = solvent_finder_origin(full_library_mongo, solvent_dict_stub, '')
                    if return_statement != 'Success':
                        return ['Error', 'Problem: Not able to find details for %s' % well_solvent, success_statements]
                    reaction_plate_to_prepare[well_key]['solvents'].append(
                        [solvent_details[well_solvent]['chemical_name'],
                         well_solvents[well_solvent],
                         solvent_details[well_solvent]['chemical_smiles']])
                else:
                    reaction_plate_to_prepare[well_key]['solvents'].append(
                        [solvents_needed[well_solvent]['chemical_name'],
                         well_solvents[well_solvent],
                         solvents_needed[well_solvent]['chemical_smiles']])
            reaction_plate_to_prepare[well_key].pop('reagents_info', None)
            reaction_plate_to_prepare[well_key].pop('previous_product_info', None)
            reaction_plate_to_prepare[well_key].pop('solvents_info', None)
            reaction_plate_to_prepare[well_key].pop('remaining_volume', None)

        if full_library_mongo['wellplates'][str(reaction_plate_initial['_id'])]['contents'] == 'Empty':
            new_wellplate_contents = reaction_plate_to_prepare
        else:
            new_wellplate_contents = full_library_mongo['wellplates'][str(reaction_plate_initial['_id'])]['contents']
            for well_key in reaction_plate_to_prepare.keys():
                if well_key not in new_wellplate_contents.keys():
                    new_wellplate_contents[well_key] = reaction_plate_to_prepare[well_key]
                else:
                    keys_to_merge = ['reagents', 'solvents']
                    for key_to_merge in keys_to_merge:
                        if key_to_merge in reaction_plate_to_prepare[well_key].keys() and key_to_merge in \
                                new_wellplate_contents[well_key].keys():
                            new_wellplate_contents[well_key][key_to_merge].extend(
                                reaction_plate_to_prepare[well_key][key_to_merge])
                        if key_to_merge in reaction_plate_to_prepare[well_key].keys() and key_to_merge not in \
                                new_wellplate_contents[well_key].keys():
                            new_wellplate_contents[well_key][key_to_merge] = reaction_plate_to_prepare[well_key][
                                key_to_merge]
                    keys_to_overwrite = ['total_volume']
                    for key_to_overwrite in keys_to_overwrite:
                        if key_to_overwrite in reaction_plate_to_prepare[well_key].keys():
                            new_wellplate_contents[well_key][key_to_overwrite] = reaction_plate_to_prepare[well_key][
                                key_to_overwrite]
            
        db_reaction_plate = full_library_mongo['wellplates'][str(reaction_plate_initial['_id'])]
        db_reaction_plate['contents'] = new_wellplate_contents
        db_reaction_plate['container_category'] = additional_details['target_container']
        db_reaction_plate['date_updated'] = current_time.strftime('%m/%d/%Y')
        success_statements.append('The completed well-plate is located at: %s' % \
                                  reaction_plate_initial['location'][1])
        
        # pprint(updated_carriers)
        # Update the amounts of all the solvents and reagents used in the plate prep
        for solvent_needed in solvents_needed.keys():
            solvent_details = solvents_needed[solvent_needed]
            full_library_mongo[solvent_details['collection']][solvent_details['container_id']]['contents'][
                solvent_details['origin_well_location']]['volume_ul'] -= solvent_details['volume_needed']
        for reagent_needed in reagents_needed.keys():
            reagent_details = reagents_needed[reagent_needed]
            full_library_mongo[reagent_details['collection']][reagent_details['container_id']]['contents'][
                reagent_details['plate_well']]['volume_ul'] -= deepcopy(reagent_details['volume_needed'])

        if __name__ == "__main__":
            pass
            for reagent_needed in reagents_needed.keys():
                print(reagent_needed, reagents_needed[reagent_needed]['container_name'],
                      reagents_needed[reagent_needed]['plate_well'])
            # pprint(reagents_needed, width=125)
            return ['Testing!', 'Testing!',success_statements]

        platform_mongodb = pymongo.MongoClient(additional_details['platform_mongodb_address'],
                                               additional_details['platform_mongodb_port'])
        for library_key in full_library_mongo.keys():
            previous = initial_full_library[library_key]
            updated = full_library_mongo[library_key]
            return_statements = update_database(previous, updated, library_key, platform_mongodb)

        updated_queue_document = copy.deepcopy(queue_document)
        if __name__ == "__main__":
            updated_queue_document['operations'][operation_number]['completed'] = 'yes'
        reaction_plate_name = reaction_plate_initial['container_name']
        target_container = updated_queue_document['operations'][operation_number]['container']
        updated_queue_document['containers'][target_container]['container_name'] = reaction_plate_name
        flag_return, return_statements = update_database_document(queue_document, updated_queue_document, 'queue',
                                                                  platform_mongodb)
        platform_mongodb.close()
        success_statements.append('Updated the platform library to be up-to-date following this run')
        return ['Success', additional_details['schedule_path_output'], success_statements]
    except Exception:
        return ['Error', 'Problem: %s' % traceback.format_exc(), success_statements]


if __name__ == "__main__":
    queuename = 'photochem_discovery_20230519_redo_20230620_1'
    prep_return = reaction_plate_preparation(queuename, '10')
    pprint(prep_return, width=150)
