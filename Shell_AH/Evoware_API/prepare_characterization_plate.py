# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 10:04:31 2022

@author: Brent
"""

import copy
import pymongo
import traceback
import math
from pprint import pprint
from mergedeep import merge
from copy import deepcopy
import sys
import numpy as np

if __name__ == "__main__":
    from supporting_functions import get_queue_information, convert_to_method_strings
    from platform_library_functions import update_database, update_database_document
    from useful_other_functions import solvent_finder, priority_sort, liha_grouping, solvent_finder_origin, \
        find_accessible_open_location, check_for_problems, chemical_characterization_lookup, find_labware, reagent_finder
    from worktable_cleaner import initial_worktable_prep, worktable_export
else:
    from Evoware_API.supporting_functions import get_queue_information, convert_to_method_strings
    from Evoware_API.platform_library_functions import update_database, update_database_document
    from Evoware_API.useful_other_functions import solvent_finder, priority_sort, liha_grouping, solvent_finder_origin, \
        find_accessible_open_location, check_for_problems, chemical_characterization_lookup, find_labware, reagent_finder
    from Evoware_API.worktable_cleaner import initial_worktable_prep, worktable_export


def characterization_plate_preparation(queue_name, operation_number):
    success_statements = []
    try:
        # ----------------------------------------------------------------------------
        # First we start with opening the queue document and getting operation info.
        queue_information = get_queue_information(queue_name, operation_number, 'prepare_characterization_plate')
        if queue_information[0] != 'Success':
            return ['Error', 'Problem: %s' % queue_information[1], success_statements]
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
            return ['Error', 'Problem: %s' % start_carriers, success_statements]
        updated_carriers = copy.deepcopy(start_carriers)
        success_statements.append('Loaded the carrier and labware library file!')
        # ----------------------------------------------------------------------------

        # Start off by sorting out the plate that we need to prepare
        well_plate_limit = 96
        contents_to_prepare = additional_details['queue_containers'][additional_details['target_container']]['contents']
        if 'targets_to_find' in additional_details['queue_containers'][additional_details['target_container']].keys():
            additional_wells_to_prepare = additional_details['queue_containers'][additional_details['target_container']]['targets_to_find']
        else:
            additional_wells_to_prepare = {}
        allowed_wells = [chr(ord('A') + ind//12) + str(ind%12 + 1) for ind in range(0, 96)]

        additional_wells_prepared = []
        missing_chemicals = []
        for target_index in additional_wells_to_prepare.keys():
            for allowed_well in allowed_wells:
                if allowed_well not in contents_to_prepare.keys():
                    next_well = allowed_well
                    break
            else:
                print('There are no more spaces for the other chemicals')
                continue
            reagent_needed = additional_wells_to_prepare[target_index]['target_molecules'][0]
            reagent = reagent_needed[0]
            if len(reagent_needed) == 4:
                destination = [next_well, *reagent_needed[1:], 1]
            else:
                destination = [next_well, *reagent_needed[1:]]
            reagent_details = {'destination': destination, 'container_category': reagent_needed[3]}
            return_statement, reagent_information = chemical_characterization_lookup(full_library_mongo, reagent,
                                                                                     reagent_details)
            if return_statement == 'Missing Chemical':
                missing_chemicals.append(target_index)
                continue
            elif return_statement == 'Error':
                return ['Error', 'Problem: %s' % reagent_information, success_statements]
            if next_well in allowed_wells and next_well not in contents_to_prepare.keys():
                contents_to_prepare[next_well] = additional_wells_to_prepare[target_index]
                additional_wells_prepared.append(target_index)

        reagents_needed = {}
        additional_dilutions_needed = {}
        problematic_chemicals = []
        for well_location in contents_to_prepare.keys():
            contents_to_prepare[well_location]['used_volume'] = 0
            target_molecules = contents_to_prepare[well_location]['target_molecules']
            for target_molecule in target_molecules:
                reagent_details = {'destination': target_molecule, 'container_category': target_molecule[3]}
                return_statement, reagent_information = chemical_characterization_lookup(full_library_mongo,
                                                                                         target_molecule[0],
                                                                                         reagent_details)
                if return_statement == 'Missing Chemical':
                    problematic_chemicals.append('Missing required chemical: %s' % target_molecule[0])
                elif return_statement == 'Error':
                    return ['Error', reagent_information, success_statements]
                if 'reagent_information' not in contents_to_prepare[well_location].keys():
                    contents_to_prepare[well_location]['reagent_information'] = {}
                labware_details = full_library_mongo[reagent_information[0]][reagent_information[1]]
                reagent_well_details = labware_details[reagent_information[2]][reagent_information[3]]

                prep_type = target_molecule[2]
                transfers_needed = []
                if target_molecule[3] == 'reagent_tray':
                    current_stock_concentration = reagent_well_details['concentration_molar']
                    current_stock_volume = reagent_well_details['volume_ul']
                else:
                    well_products = reagent_well_details['target_product']
                    relevant_product = []
                    for well_product in well_products:
                        if well_product[0] == target_molecule[0]:
                            relevant_product = well_product
                            break
                    if not relevant_product:
                        problematic_chemicals.append('Problem finding %s in well %s' % (target_molecule[0],
                                                                                        well_location))
                    current_stock_volume = sum([solvent[1] for solvent in reagent_well_details['solvents']])
                    current_stock_volume = 0
                    current_stock_volume = reagent_well_details['total_volume']
                    current_stock_concentration = relevant_product[1] / (current_stock_volume * 1E-6)

                potential_well_volume = contents_to_prepare[well_location]['total_volume']
                if prep_type.lower() == 'od':
                    if 'properties' not in reagent_well_details.keys():
                        problematic_chemicals.append('Problem: %s and %s, no properties key' % (prep_type,
                                                                                                target_molecule[0]))
                        continue
                    try:
                        extinct_coef = 10 ** reagent_well_details['properties'][target_molecule[0]]['log_e']['CC#N'][0]
                    except Exception:
                        extinct_coef = 1E4
                        
                    target_optical_density = float(target_molecule[1])
                    pathlength = (1/10) * potential_well_volume / (math.pi * (7/2)**2)
                    target_concentration = target_optical_density / (extinct_coef * pathlength)
                    dilution_needed = target_concentration / current_stock_concentration
                    potential_transfer_volume = potential_well_volume * dilution_needed
                elif prep_type.lower() == 'concentration':
                    target_concentration = float(target_molecule[1])
                    dilution_needed = target_concentration / current_stock_concentration
                    potential_transfer_volume = potential_well_volume * dilution_needed
                elif prep_type.lower() == 'amount':
                    target_mols = float(target_molecule[1])
                    target_concentration = target_mols / (potential_well_volume * 1E-6)
                    dilution_needed = target_concentration / current_stock_concentration
                    potential_transfer_volume = potential_well_volume * dilution_needed
                else:
                    problematic_chemicals.append('Problem with %s, prep type %s not implemented' % (target_molecule[0],
                                                                                                    prep_type))
                if dilution_needed > 1:
                    problematic_chemicals.append('Problem: %s, dilution needed is bad %s' % (target_molecule[0],
                                                                                             dilution_needed))
                    continue
                if potential_transfer_volume < 5:
                    dilution_sequence = []
                    current_concentration = current_stock_concentration
                    for dilution in range(0, 8):
                        if dilution == 0:
                            dilution_factor = 5 / potential_well_volume
                            current_concentration = dilution_factor * current_concentration
                            # dilution_sequence.append((1 - dilution_factor) * potential_well_volume)
                        else:
                            if target_concentration > 0.1 * current_concentration:
                                dilution_factor = target_concentration / current_concentration
                                dilution_sequence.append((1 - dilution_factor) * potential_well_volume)
                                break
                            else:
                                current_concentration = 0.1 * current_concentration
                                dilution_sequence.append(0.1 * potential_well_volume)
                    else:
                        problematic_chemicals.append('Problem with %s, dilution is too much %s' % (target_molecule[0],
                                                                                                   potential_transfer_volume))
                        continue
                    if well_location not in additional_dilutions_needed.keys():
                        additional_dilutions_needed[well_location] = {'initial_volume': 5,
                                                                      'dilution_sequence': dilution_sequence,
                                                                      'target_chemical': target_molecule[0],
                                                                      'chemical_type': 'reagent',
                                                                      'total_volume': potential_well_volume}
                    else:
                        problematic_chemicals.append('Well location %s needs multiple components to be \
                        diluted too much' % well_location)
                        continue
                    potential_transfer_volume = 5
                # potential_transfer_volume = 200
                transfers_needed.append([well_location, reagent_information[0], reagent_information[1],
                                         reagent_information[-1], potential_transfer_volume])
                contents_to_prepare[well_location]['used_volume'] += potential_transfer_volume
                transfer_sequence = target_molecule[-1]
                if target_molecule[0] not in reagents_needed.keys():
                    reagents_needed[target_molecule[0]] = deepcopy(reagent_well_details)
                    reagents_needed[target_molecule[0]].update({'volume_needed': 0,
                                                                'collection': reagent_information[0],
                                                                'container_id': reagent_information[1],
                                                                'container_name': labware_details['container_name'],
                                                                'destinations': {}})
                if transfer_sequence not in reagents_needed[target_molecule[0]]['destinations'].keys():
                    reagents_needed[target_molecule[0]]['destinations'][transfer_sequence] = []
                reagents_needed[target_molecule[0]]['destinations'][transfer_sequence].append([well_location,
                                                                                               potential_transfer_volume,
                                                                                               transfer_sequence])
                reagents_needed[target_molecule[0]]['volume_needed'] += potential_transfer_volume
        
        print(reagents_needed)
        print(additional_dilutions_needed)
        solvents_needed = {}
        for well_location in contents_to_prepare.keys():
            used_volume = contents_to_prepare[well_location]['used_volume']
            total_volume = contents_to_prepare[well_location]['total_volume']
            solvent_scaling = sum([solvent_needed[1] for solvent_needed in
                                   contents_to_prepare[well_location]['solvents']])
            for solvent_needed in contents_to_prepare[well_location]['solvents']:
                transfer_fraction = solvent_needed[1] / solvent_scaling
                transfer_volume = (total_volume - used_volume) * transfer_fraction
                if len(solvent_needed) == 2:
                    transfer_sequence = 1
                else:
                    transfer_sequence = solvent_needed[2]
                return_statement, solvent_details = solvent_finder_origin(full_library_mongo, [[solvent_needed[0],
                                                                                                transfer_volume]], '')
                if return_statement != "Success":
                    return ['Error', 'Problem: %s' % solvent_details, success_statements]
                for solvent_return in solvent_details:
                    if solvent_return not in solvents_needed.keys():
                        solvents_needed[solvent_return] = solvent_details[solvent_return]
                        solvents_needed[solvent_return].update({'destinations': {}})
                    if transfer_sequence not in solvents_needed[solvent_return]['destinations'].keys():
                        solvents_needed[solvent_return]['destinations'][transfer_sequence] = []
                    solvents_needed[solvent_return]['destinations'][transfer_sequence].append([well_location,
                                                                                               transfer_volume,
                                                                                               transfer_sequence])
                    solvents_needed[solvent_return]['volume_needed'] += transfer_volume
        if problematic_chemicals:
            return ['Error', 'Problem: %s' % problematic_chemicals, success_statements]

        # Now we need to figure out the transfer sequence for each of the groups
        return_statement, solvents_order = priority_sort(solvents_needed, 'solvents')
        if return_statement == 'Error':
            return ['Error', 'Problem: %s' % solvents_order, success_statements]

        return_statement, reagents_order = priority_sort(reagents_needed, 'previous_products')
        if return_statement == 'Error':
            return ['Error', reagents_order, success_statements]

        # This gives us all the chemicals that we need to distribute into wells, so now
        # we can build up the schedule of tasks needed to distribute the chemicals, we
        # need to get the locations that we are placing the plates at and also find an
        # empty reaction plate that we can use for the method
        simplified_schedule, method_strings = [], []
        return_statement, wellplate_initial = find_labware(updated_carriers, full_library_mongo,
                                                           additional_details['wellplate_name'],
                                                           additional_details['target_well_plate_type'],
                                                           additional_details['target_container'])
        if return_statement != 'Success':
            return ['Error', 'Problem: %s' % wellplate_initial, success_statements]
        success_statements.append('Found an acceptable reaction plate %s at %s' % (wellplate_initial['container_name'],
                                                                                   wellplate_initial['location'][1]))
        allowed_prep_carriers = [32, 38, 44]
        if wellplate_initial['location'][1][0] in allowed_prep_carriers:
            wellplate_home = wellplate_initial['location'][1]
        else:
            return_statement, wellplate_home = find_accessible_open_location(updated_carriers)
            if return_statement == 'Error':
                return ['Error', 'Problem: %s' % wellplate_home, success_statements]
        success_statements.append('Found a suitable prep location at location: %s' % wellplate_home)

        if wellplate_initial['location'][1] != wellplate_home:
            simplified_schedule.append(['TRANSFER_LABWARE', wellplate_initial['container_name'],
                                        wellplate_initial['location'][1], wellplate_home, str(wellplate_initial['_id']),
                                        'wellplates'])
            updated_carriers[wellplate_home[0]][wellplate_home[2]]['labware_labels'][wellplate_home[1] - 1] = \
                'location_%s_%s_%s' % (str(wellplate_home[0]), str(wellplate_home[1]), str(wellplate_home[2]))
            updated_carriers[wellplate_home[0]][wellplate_home[2]]['labware_types'][wellplate_home[1] - 1] = \
                wellplate_initial['labware_type']

        # Now the aspirate and dispense functions are added to the simplified schedule
        if additional_prep_details['tip_prep_rinse']:
            simplified_schedule.append(['TIP_RINSE', 10, 10])
        else:
            simplified_schedule.append(['TIP_RINSE', 5, 5])
        simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 150])
        simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 150])

        for solvent_seq in solvents_order.keys():
            solvent_index = 0
            for target_container in solvents_order[solvent_seq].keys():
                relevant_collections = ['solvents', 'reagents']
                for collection in relevant_collections:
                    if target_container in full_library_mongo[collection].keys():
                        relevant_collection = collection
                        break
                else:
                    return ['Error', 'Unable to find solvent plate: %s' % target_container, success_statements]
                solvent_definition = full_library_mongo[relevant_collection][target_container]
                solvent_origin = solvent_definition['location']
                origin_labware_type = carriers_labware['labware'][solvent_definition['labware_type']]
                destination_labware_type = carriers_labware['labware'][wellplate_initial['labware_type']]
                for target_solvent in solvents_order[solvent_seq][target_container].keys():
                    if solvent_definition['labware_type'] in ['8 Position Vial Carrier', '12 Position Vial Carrier']:
                        solvent_origin_well = solvents_needed[target_solvent]['initial_solvent_details']['plate_well']
                        for target_destination in solvents_order[solvent_seq][target_container][target_solvent]['destinations']:
                            tips = [1]
                            volumes = [target_destination[1]]
                            wells = [target_destination[0]]
                            simplified_schedule.append(['LiHa_ASPIRATE_DISPENSE', solvent_origin[1] + [[solvent_origin_well]],
                                                        wellplate_home, tips, volumes, wells])
                    else:
                        destinations = solvents_order[solvent_seq][target_container][target_solvent]['destinations']
                        return_statement, liha_groups = liha_grouping(destinations, carriers_labware,
                                                                      solvent_definition['labware_type'],
                                                                      wellplate_initial['labware_type'])
                        if return_statement == 'Error':
                            return ['Error', 'Something went wrong with LiHa grouping', success_statements]
                        for group_no, liha_group in enumerate(liha_groups.keys()):
                            group = liha_groups[liha_group]
                            wells, volumes, tips = map(list, zip(*group))
                            simplified_schedule.append(
                                ['LiHa_ASPIRATE_DISPENSE', solvent_origin[1], wellplate_home, tips, volumes, wells])
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

        if 'solvent_waiting_time' in additional_prep_details.keys():
            try:
                solvent_waiting_time = additional_prep_details['solvent_waiting_time']
                simplified_schedule.append(['START_TIMER', 1])
                simplified_schedule.append(['WAIT_FOR_TIMER', 1, float(solvent_waiting_time)])
            except:
                success_statements.append('Skipping waiting time, problem with: %s' %
                                          additional_prep_details['solvent_waiting_time'])

        acceptable_bed_locations = [32, 38, 44]
        location_occupancies = []
        liha_tip_number = 1
        
        if reagents_order:
            for reagent_seq in reagents_order.keys():
                for reagent_container in reagents_order[reagent_seq].keys():
                    relevant_collections = ['solvents', 'reagents', 'wellplates']
                    for collection in relevant_collections:
                        if reagent_container in full_library_mongo[collection].keys():
                            relevant_collection = collection
                            break
                    else:
                        return ['Error', 'Problem: Unable to find reagents plate: %s' % reagent_container,
                                success_statements]
                    wellplate_definition = full_library_mongo[relevant_collection][reagent_container]
                    current_plate_location = wellplate_definition['location'][1]
                    success_statements.append('Currently the reagent tray "%s" is located at location: %s' % (
                        wellplate_definition['container_name'], current_plate_location))
                    if wellplate_definition['location'][1][0] not in acceptable_bed_locations:
                        exceptions = [location[0] for location in location_occupancies if
                                      location[1] != wellplate_definition['labware_type']]
                        return_statement, new_wellplate_home = find_accessible_open_location(updated_carriers,
                                                                                             exceptions,
                                                                                             wellplate_definition[
                                                                                                 'labware_type'])
                        if return_statement == 'Error':
                            return ['Error',
                                    'Reagents well-plate %s not accessible' % wellplate_definition['container_name'],
                                    success_statements]

                        if wellplate_definition['location'] != new_wellplate_home:
                            simplified_schedule.append(['TRANSFER_LABWARE', wellplate_definition['container_name'],
                                                        wellplate_definition['location'][1], new_wellplate_home,
                                                        str(wellplate_definition['_id']), relevant_collection])
                            current_plate_location = new_wellplate_home
                    location_occupancies.append([current_plate_location, wellplate_definition['labware_type']])
                    for target_product in reagents_order[reagent_seq][reagent_container].keys():
                        product_origin_well = reagents_needed[target_product]['plate_well']
                        for target_destination in reagents_order[reagent_seq][reagent_container][target_product]['destinations']:
                            if target_destination[0] in additional_dilutions_needed.keys():
                                if additional_dilutions_needed[target_destination[0]]['target_chemical'] == target_product:
                                    dilution_sequence = additional_dilutions_needed[target_destination[0]]['dilution_sequence']
                                    number_of_dilutions = len(dilution_sequence)
                                    tips = [liha_tip_number]
                                    volumes = [target_destination[1]]
                                    wells = [target_destination[0]]
                                    simplified_schedule.append(['LiHa_ASPIRATE_DISPENSE', current_plate_location + [[product_origin_well]],
                                                                wellplate_home, tips, volumes, wells])

                                    simplified_schedule.append(['MIX_LIHA', wellplate_home, tips, wells, [200], 10])
                                    if len(contents_to_prepare[target_destination[0]]['solvents']) > 1:
                                        return ['Error', 'Too many solvents needed for dilution in %s' % target_destination[0]]
                                    relevant_solvent = contents_to_prepare[target_destination[0]]['solvents'][0][0]
                                    solvent_origin_well = solvents_needed[relevant_solvent]['origin_well_location']
                                    solvent_definition = full_library_mongo[solvents_needed[relevant_solvent]['collection']][solvents_needed[relevant_solvent]['container_id']]
                                    solvent_origin = solvent_definition['location']

                                    for seq_num, dilution_seq in enumerate(dilution_sequence):
                                        cleanup_volumes = [dilution_seq]
                                        cleanup_wells = [chr(ord('A') + liha_tip_number) + '1']
                                        simplified_schedule.append(['LiHa_ASPIRATE_DISPENSE', wellplate_home + [[target_destination[0]]],
                                                                    [1, 2, 30], tips, cleanup_volumes, cleanup_wells])

                                        if seq_num == 0 and liha_tip_number + number_of_dilutions + 1 >= 8:
                                            simplified_schedule.append(['TIP_RINSE', 5, 5])
                                            simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
                                            simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
                                            liha_tip_number = 1
                                        else:
                                            liha_tip_number += 1
                                        tips = [liha_tip_number]
                                        volumes = [dilution_seq]
                                        wells = [target_destination[0]]
                                        simplified_schedule.append(['LiHa_ASPIRATE_DISPENSE', solvent_origin[1] + [[solvent_origin_well]],
                                                            wellplate_home, tips, volumes, wells])
                                        solvents_needed[relevant_solvent]['volume_needed'] += dilution_seq
                                        simplified_schedule.append(['MIX_LIHA', wellplate_home, tips, wells, [200], 10])

                                    solvent_existing_volume = solvents_needed[relevant_solvent]['solvent_volume']
                                    if solvents_needed[relevant_solvent]['volume_needed'] >= solvent_existing_volume:
                                        return ['Error', 'Not enough %s present on the platform' % relevant_solvent]
                            else:
                                tips = [liha_tip_number]
                                volumes = [target_destination[1]]
                                wells = [target_destination[0]]
                                simplified_schedule.append(['LiHa_ASPIRATE_DISPENSE', current_plate_location + [[product_origin_well]],
                                                        wellplate_home, tips, volumes, wells])
                                simplified_schedule.append(['MIX_LIHA', wellplate_home, tips, wells, [200], 5])
                            if liha_tip_number == 8:
                                simplified_schedule.append(['TIP_RINSE', 5, 5])
                                simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
                                simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
                                liha_tip_number = 1
                            else:
                                liha_tip_number += 1
                    simplified_schedule.append(['MOVE_LiHa'])
                    if wellplate_definition['location'][1] != current_plate_location:
                        problem_statement = check_for_problems(updated_carriers, current_plate_location,
                                                               wellplate_definition['container_name'])
                        if problem_statement != 'Acceptable':
                            return ['Error', 'Cannot access the reagent tray %s for back transfer to its home' %
                                    wellplate_definition['container_name'], success_statements]
                        simplified_schedule.append(['TRANSFER_LABWARE', wellplate_definition['container_name'],
                                                    current_plate_location, wellplate_definition['location'][1],
                                                    str(wellplate_definition['_id']), relevant_collection])
                    if len(list(reagents_order.keys())) != 1:
                        if 'addition_sequence_waiting' in additional_prep_details.keys() and additional_prep_details['addition_sequence_waiting'] == True:
                            if 'reagent_waiting_time' in additional_prep_details.keys():
                                try:
                                    reagent_waiting_time = float(additional_prep_details['reagent_waiting_time'])
                                except:
                                    reagent_waiting_time = 60
                            else:
                                reagent_waiting_time = 60
                            simplified_schedule.append(['START_TIMER', 1])
                            simplified_schedule.append(['WAIT_FOR_TIMER', 1, reagent_waiting_time])

                simplified_schedule.append(['TIP_RINSE', 5, 5])
                simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
                simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
                simplified_schedule.append(['TIP_RINSE', 4, 4])

        if wellplate_home != wellplate_initial['location'][1]:
            problem_statement = check_for_problems(updated_carriers, wellplate_home, wellplate_initial['container_name'])
            if problem_statement != 'Acceptable':
                return ['Error', 'Cannot access the reaction plate %s for back transfer to its home' % wellplate_initial['container_name'],
                        success_statements]
            simplified_schedule.append(['TRANSFER_LABWARE', wellplate_initial['container_name'],
                                        wellplate_home, wellplate_initial['location'][1],
                                        str(wellplate_initial['_id']), 'wellplates'])

        # Now we need to rebuild the worktable, so we can write a file to run
        return_statement, filled_worktable_export = worktable_export(additional_details['carriers_labware_filepath'],
                                                                     additional_details['empty_worktable_filepath'],
                                                                     updated_carriers,
                                                                     full_library_mongo,
                                                                     simplified_schedule)
        if return_statement == 'Error':
            return ['Error', filled_worktable_export, success_statements]
        # pprint(simplified_schedule, width=200)
        # We need to convert the simplified schedule into Tecan method strings
        return_statement, method_strings, updated_carriers = convert_to_method_strings(simplified_schedule, updated_carriers, full_library_mongo, carriers_labware)
        if return_statement == 'Error':
            return ['Error', method_strings, success_statements]

        with open(additional_details['schedule_path_output'], 'w', encoding='latin1') as outfile:
            for line in filled_worktable_export:
                outfile.write(line + '\n')
            for line in method_strings:
                outfile.write(line + '\n')

        if full_library_mongo['wellplates'][str(wellplate_initial['_id'])]['contents'] == 'Empty':
            new_wellplate_contents = {}
            for key in contents_to_prepare.keys():
                new_wellplate_contents[key] = {}
        else:
            new_wellplate_contents = full_library_mongo['wellplates'][str(wellplate_initial['_id'])]['contents']

        for reagent in reagents_needed.keys():
            reagent_information = reagents_needed[reagent]
            if reagent_information['collection'] == 'wellplates':
                reagent_name = reagent
                reagent_smiles = reagent
                previous_well_details = full_library_mongo[reagent_information['collection']][
                    reagent_information['container_id']]['contents'][
                    reagent_information['plate_well']]
                previous_products = previous_well_details['target_product']
                transfer_fraction = reagent_information['volume_needed'] / reagent_information['total_volume']
                relevant_product = []
                for previous_product in previous_products:
                    if previous_product[0] == reagent:
                        relevant_product = deepcopy(previous_product)
                    previous_product[1] = previous_product[1] * (1 - transfer_fraction)
                if not relevant_product:
                    return ['Error', 'Problem: not able to find %s' % reagent, success_statements]
            else:
                previous_well_details = full_library_mongo[reagent_information['collection']][
                    reagent_information['container_id']]['contents'][reagent_information['plate_well']]
                reagent_name = previous_well_details['chemical_name']
                reagent_smiles = previous_well_details['chemical_smiles']
                transfer_fraction = reagent_information['volume_needed'] / reagent_information['volume_ul']

            for transfer_sequence in reagent_information['destinations'].keys():
                for destination in reagent_information['destinations'][transfer_sequence]:
                    if reagent_information['collection'] == 'wellplates':
                        reagent_moles = relevant_product[1] * destination[1] / previous_well_details['total_volume']
                        solvents = reagent_information['solvents']
                        reagent_total_volume = reagent_information['total_volume']
                    else:
                        reagent_amount = reagent_information['volume_ul'] * reagent_information['concentration_molar'] * 1E-6
                        reagent_moles = reagent_amount * destination[1] / reagent_information['volume_ul']
                        solvents = [[reagent_information['solvent'], reagent_information['volume_ul'],
                                     reagent_information['solvent']]]
                        reagent_total_volume = reagent_information['volume_ul']

                    updated_solvents = {}
                    if destination[0] in additional_dilutions_needed.keys():
                        dilution_details = additional_dilutions_needed[destination[0]]
                        additional_dilution = np.prod([1 - (volume / dilution_details['total_volume']) for
                                                       volume in dilution_details['dilution_sequence'][1:]])
                    else:
                        additional_dilution = 1

                    for solvent in solvents:
                        composition_fraction = solvent[1] / reagent_total_volume
                        solvent_volume = composition_fraction * destination[1] * additional_dilution
                        if solvent[0] not in updated_solvents.keys():
                            updated_solvents[solvent[0]] = {'volume': solvent_volume, 'name': solvent[0],
                                                            'smiles': solvent[2]}
                        else:
                            updated_solvents[solvent[0]]['volume'] += solvent_volume
                    used_volume = sum([updated_solvents[solvent]['volume'] for solvent in updated_solvents.keys()])
                    remaining_volume = contents_to_prepare[destination[0]]['total_volume'] - used_volume
                    solvent_norm = sum([solvent_needed[1] for solvent_needed in contents_to_prepare[destination[0]]['solvents']])
                    for solvent_needed in contents_to_prepare[destination[0]]['solvents']:
                        volume_fraction = solvent_needed[1] / solvent_norm
                        if solvent_needed[0] not in updated_solvents.keys():
                            updated_solvents[solvent_needed[0]] = {
                                'volume': remaining_volume * volume_fraction,
                                'name': solvents_needed[solvent_needed[0]]['chemical_name'],
                                'smiles': solvents_needed[solvent_needed[0]]['chemical_smiles']}
                        else:
                            updated_solvents[solvent_needed[0]]['volume'] += (remaining_volume * volume_fraction)
                    solvents_list = []
                    for solvent in updated_solvents.keys():
                        solvents_list.append([updated_solvents[solvent]['name'],
                                              updated_solvents[solvent]['volume'],
                                              updated_solvents[solvent]['smiles']])
                    if destination[0] not in new_wellplate_contents.keys():
                        new_wellplate_contents[destination[0]] = {}
                    if 'reagents' not in new_wellplate_contents[destination[0]].keys():
                        new_wellplate_contents[destination[0]]['reagents'] = []
                    if 'solvents' not in new_wellplate_contents[destination[0]].keys():
                        new_wellplate_contents[destination[0]]['solvents'] = solvents_list
                    if 'target_molecules' not in new_wellplate_contents[destination[0]].keys():
                        new_wellplate_contents[destination[0]]['target_molecules'] = []
                    new_wellplate_contents[destination[0]]['target_molecules'].append(reagent_smiles)
                    new_wellplate_contents[destination[0]]['reagents'].append([reagent_name, reagent_moles,
                                                                               reagent_smiles])
            reagent_details = full_library_mongo[reagent_information['collection']][reagent_information[
                'container_id']]['contents'][reagent_information['plate_well']]
            if reagent_information['collection'] == 'wellplates':
                updated_well_solvents = []
                for solvent in reagent_details['solvents']:
                    updated_well_solvents.append([solvent[0], solvent[1] * (1 - transfer_fraction), solvent[2]])
                reagent_details['solvents'] = updated_well_solvents
                reagent_details['total_volume'] *= (1 - transfer_fraction)
                updated_target_products = []
                for target_product in reagent_details['target_product']:
                    updated_target_products.append([target_product[0], target_product[1] * (1 - transfer_fraction)])
                reagent_details['product_remaining'] = updated_target_products
            else:
                reagent_details['volume_ul'] -= reagent_information['volume_needed']

        for well_key in new_wellplate_contents.keys():
            new_wellplate_contents[well_key]['total_volume'] = contents_to_prepare[well_key]['total_volume']
        full_library_mongo['wellplates'][str(wellplate_initial['_id'])]['contents'] = new_wellplate_contents
        full_library_mongo['wellplates'][str(wellplate_initial['_id'])]['container_category'] = additional_details['target_container']
        success_statements.append('The completed characterization well-plate is located at: %s' % wellplate_initial['location'][1])

        chemicals_left_to_find = list(set(list(additional_wells_to_prepare.keys())) - set(additional_wells_prepared) - set(missing_chemicals))
        remaining_targets = {}
        for chemical in chemicals_left_to_find:
            remaining_targets[chemical] = additional_wells_to_prepare[chemical]
        if remaining_targets:
            container_index = 2
            while True:
                if 'characterization_plate_%s' % container_index in queue_document['containers'].keys():
                    container_index += 1
                else:
                    new_container_name = 'characterization_plate_%s' % container_index
                    break
                if container_index > 100:
                    return ['Error', 'Something went wrong with container index']

            pprint(remaining_targets)
            queue_document['containers'][new_container_name] = {'container_name': None,
                                            'plate_type': queue_document['containers'][additional_details['target_container']]['plate_type'],
                                            'targets_to_find': remaining_targets}
            operations_to_copy = additional_prep_details['copy_next_n_operations']
            print(operations_to_copy)
            copied_operations = {}
            for operation_index in range(0, operations_to_copy + 1):
                operation_key = str(int(operation_number) + operation_index)
                new_operation_key = str(int(operation_number) + operation_index + operations_to_copy + 1)
                copy_operation = deepcopy(queue_document['operations'][operation_key])
                copy_operation['container'] = new_container_name
                copied_operations[new_operation_key] = copy_operation
            pprint(queue_document['operations'].keys())

            updated_operations = {}
            for operation_index in range(1, len(list(queue_document['operations'].keys())) + 1):
                if operation_index < int(operation_number) + operations_to_copy + 1:
                    updated_operations[str(operation_index)] = queue_document['operations'][str(operation_index)]
                elif operation_index == int(operation_number) + operations_to_copy + 1:
                    for new_operation in copied_operations.keys():
                        updated_operations[new_operation] = copied_operations[new_operation]
                    new_operation_index = operation_index + len(list(copied_operations.keys()))
                    updated_operations[str(new_operation_index)] = queue_document['operations'][str(operation_index)]
                else:
                    new_operation_index = operation_index + len(list(copied_operations.keys()))
                    updated_operations[str(new_operation_index)] = queue_document['operations'][str(operation_index)]

            queue_document['operations'] = updated_operations

        if __name__ == "__main__":
            pass
            return ['Testing!!', success_statements]

        platform_mongodb = pymongo.MongoClient(additional_details['platform_mongodb_address'],
                                               additional_details['platform_mongodb_port'])
        for library_key in full_library_mongo.keys():
            previous = initial_full_library[library_key]
            updated = full_library_mongo[library_key]
            return_statements = update_database(previous, updated, library_key, platform_mongodb)

        updated_queue_document = copy.deepcopy(queue_document)
        if __name__ == "__main__":
            updated_queue_document['operations'][operation_number]['completed'] = 'yes'
        wellplate_name = wellplate_initial['container_name']
        target_container = updated_queue_document['operations'][operation_number]['container']
        updated_queue_document['containers'][target_container]['container_name'] = wellplate_name
        flag_return, return_statements = update_database_document(queue_document, updated_queue_document, 'queue', platform_mongodb)
        platform_mongodb.close()
        success_statements.append('Updated the platform library to be up-to-date following this run')
        return ['Success', additional_details['schedule_path_output'], success_statements]
    except Exception:
        return ['Error', traceback.format_exc(), success_statements]


if __name__ == "__main__":
    queuename = 'test_spark_concentration_thing'
    prep_return = characterization_plate_preparation(queuename, '1')
    pprint(prep_return, width=250)
