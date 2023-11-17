# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 09:52:53 2022

@author: Brent
"""

import copy
import pymongo
import traceback
from pprint import pprint
from mergedeep import merge
from copy import deepcopy
from scipy import optimize
import sys

if __name__ == "__main__":
    from supporting_functions import get_queue_information, convert_to_method_strings
    from platform_library_functions import update_database, update_database_document
    from useful_other_functions import solvent_finder, priority_sort, find_accessible_open_location, \
        check_for_problems, previous_product_well_origin, find_labware, reagent_finder, chemical_characterization_lookup
    from worktable_cleaner import initial_worktable_prep, worktable_export
else:
    from Evoware_API.supporting_functions import get_queue_information, convert_to_method_strings
    from Evoware_API.platform_library_functions import update_database, update_database_document
    from Evoware_API.useful_other_functions import solvent_finder, priority_sort, find_accessible_open_location, \
        check_for_problems, previous_product_well_origin, find_labware, reagent_finder, chemical_characterization_lookup
    from Evoware_API.worktable_cleaner import initial_worktable_prep, worktable_export


def cleanup_characterization_plate(queue_name, operation_number):
    success_statements = []
    try:
        # ----------------------------------------------------------------------------
        # First we start with opening the queue document and getting operation info.
        queue_information = get_queue_information(queue_name, operation_number, 'cleanup_characterization_plate')
        if queue_information[0] != 'Success':
            return ['Error', queue_information[1], success_statements]
        else:
            return_statement, full_library_mongo, queue_document, additional_details, statements = queue_information

        success_statements.extend(statements)
        initial_full_library = copy.deepcopy(full_library_mongo)
        additional_prep_details = additional_details['additional_prep_details']

        # Then build the current worktable dictionary for building output files, we are
        # also making a duplicate of the worktable, so it can be updated in this code
        return_statement, start_carriers, carriers_labware = initial_worktable_prep(additional_details['carriers_labware_filepath'],
                                                                                    additional_details['empty_worktable_filepath'],
                                                                                    full_library_mongo)
        if return_statement != 'Success':
            return ['Error', start_carriers, success_statements]
        updated_carriers = copy.deepcopy(start_carriers)
        success_statements.append('Loaded the carrier and labware library file!')
        # ----------------------------------------------------------------------------
        
        if 'target_cleanup' not in additional_prep_details.keys():
            return ['Error', 'No target_cleanup key appears in the operation details', success_statements]
        target_cleanup = additional_prep_details['target_cleanup']

        wellplate_to_clean = additional_details['wellplate_name']
        current_container = additional_details['queue_containers'][additional_details['target_container']]['contents']
        wells_to_clean = {}
        for well_key in current_container.keys():
            current_well = current_container[well_key]
            if 'additional_cleanup' in current_well.keys():
                if target_cleanup in current_well['additional_cleanup'].keys():
                    wells_to_clean[well_key] = current_well

        if len(wells_to_clean.keys()) == 0:
            return ['Error', 'No wells were selected to be cleaned', success_statements]

        return_statement, wellplate_initial = find_labware(updated_carriers, full_library_mongo, wellplate_to_clean,
                                                           additional_details['target_well_plate_type'],
                                                           additional_details['target_container'])
        if return_statement != 'Success':
            return ['Error', wellplate_initial, success_statements]
        success_statements.append('Found an acceptable reaction plate %s at %s' % (wellplate_initial['container_name'],
                                                                                   wellplate_initial['location'][1]))

        problematic_wells = []
        additional_dilutions_needed = {}
        additional_reagents_needed = {}
        reagent_sequence = {}
        for well_to_clean in wells_to_clean.keys():
            cleanup_methods = wells_to_clean[well_to_clean]['additional_cleanup'][target_cleanup]
            known_cleanups = []
            if 'dilution_factor' in cleanup_methods.keys():
                known_cleanups.append('dilution_factor')
                solvents = current_container[well_to_clean]['solvents']
                # solvents = wellplate_initial['contents'][well_to_clean]['solvents']
                if len(solvents) != 1:
                    problematic_wells.append('Problem with well %s, bad number of solvents %s' % (well_to_clean, len(solvents)))
                    continue
                relevant_solvent = solvents[0][0]
                dilution_factor = cleanup_methods['dilution_factor']
                if dilution_factor > 1:
                    problematic_wells.append('Concentrating is currently implemented as concentration factor, well %s' % well_to_clean)
                potential_well_volume = wellplate_initial['contents'][well_to_clean]['total_volume']

                current_factor = 1
                dilution_sequence = []
                for dilution in range(0, 100):
                    if dilution_factor >= 0.5 * current_factor:
                        corr_dilution_factor = dilution_factor / current_factor
                        dilution_sequence.append((1 - corr_dilution_factor) * potential_well_volume)
                        break
                    else:
                        current_factor = 0.5 * current_factor
                        dilution_sequence.append(0.5 * potential_well_volume)
                else:
                    problematic_wells.append('Too many dilutions needed for %s, %s' % (well_to_clean, len(dilution_sequence)))
                    continue
                if well_to_clean in additional_dilutions_needed.keys():
                    problematic_wells.append('Too many dilutions specified for well %s' % well_to_clean)
                    continue
                additional_dilutions_needed[well_to_clean] = {'solvents': relevant_solvent, 'dilution_sequence': dilution_sequence,
                                                              'dilution_factor': dilution_factor}
            elif 'concentration_factor' in cleanup_methods.keys():
                known_cleanups.append('concentration_factor')
                solvents = current_container[well_to_clean]['solvents']
                # solvents = wellplate_initial['contents'][well_to_clean]['solvents']
                if len(solvents) != 1:
                    problematic_wells.append('Problem with well %s, bad number of solvents %s' % (well_to_clean, len(solvents)))
                    continue
                relevant_solvent = solvents[0][0]
                concentration_factor = cleanup_methods['concentration_factor']
                if concentration_factor < 1:
                    problematic_wells.append('Dilution is currently implemented as dilution factor, well %s' % well_to_clean)
                # Find the current concentration of the sample in the characterization plate
                current_well_volume = wellplate_initial['contents'][well_to_clean]['total_volume']
                target_reagents = wells_to_clean[well_to_clean]['target_molecules']
                if len(target_reagents) > 1:
                    problematic_wells.append('Too many reagents to concentrate %s, %s' % (well_to_clean, len(target_reagents)))
                    continue
                current_reagents = wellplate_initial['contents'][well_to_clean]['reagents']
                if len(current_reagents) > 1:
                    problematic_wells.append('Too many reagents in well %s, %s' % (well_to_clean, len(current_reagents)))
                    continue
                if target_reagents[0][0] not in current_reagents[0]:
                    problematic_wells.append('Target reagent %s not in wellplate' % (target_reagents[0][0]))
                    continue
                target_reagent = target_reagents[0][0]
                target_origin = target_reagents[0][3]
                if well_to_clean not in additional_reagents_needed.keys():
                    additional_reagents_needed[well_to_clean] = {'container_category': target_origin,
                                                                 'destinations': [], 'reagent_needed': target_reagent,
                                                                 'solvents': relevant_solvent, 'concentration_factor': concentration_factor}
                else:
                    problematic_wells.append('Well %s has multiple operations needed for same clean operation' % well_to_clean)
                    continue
                if target_origin == 'reagent_tray':
                    reagent_details = {'destination': well_to_clean, 'container_category': target_origin}
                    print(target_reagent, reagent_details)
                    return_statement, previous_product_list = chemical_characterization_lookup(full_library_mongo,
                                                                                               target_reagent,
                                                                                               reagent_details)
                    if return_statement != 'Success':
                        problematic_wells.append(previous_product_list)
                        continue
                    additional_reagents_needed[well_to_clean]['initial_well_details'] = deepcopy(full_library_mongo[
                         previous_product_list[0]][previous_product_list[1]][previous_product_list[2]][
                         previous_product_list[3]])
                    additional_reagents_needed[well_to_clean]['initial_well_details']['total_volume'] = additional_reagents_needed[well_to_clean]['initial_well_details']['volume_ul']
                    additional_reagents_needed[well_to_clean]['initial_well_details']['solvents'] = additional_reagents_needed[well_to_clean]['initial_well_details']['solvent']
                    
                    additional_reagents_needed[well_to_clean]['container_name'] = deepcopy(full_library_mongo[
                         previous_product_list[0]][previous_product_list[1]]['container_name'])
                    additional_reagents_needed[well_to_clean]['container_id'] = previous_product_list[1]
                    additional_reagents_needed[well_to_clean]['collection'] = previous_product_list[0]
                    additional_reagents_needed[well_to_clean]['origin_well_location'] = previous_product_list[3]
                    # pprint(additional_reagents_needed[well_to_clean])
                    #return_statement, reagent_information = reagent_finder(full_library_mongo, [target_reagents[0]], additional_reagents_needed[well_to_clean])
                    #if return_statement != 'Success':
                    #    problematic_wells.append(reagent_information)
                    #    continue
                    #
                    #additional_reagents_needed[well_to_clean]['initial_well_details'] = reagent_information[target_reagent]
                    #additional_reagents_needed[well_to_clean]['initial_well_details']['total_volume'] = reagent_information[target_reagent]['volume_ul']
                    #additional_reagents_needed[well_to_clean]['initial_well_details']['solvents'] = reagent_information[target_reagent]['solvent']
                    
                    #merge(additional_reagents_needed[well_to_clean], reagent_information[target_reagent])
                else:
                    reagent_details = {'destination': well_to_clean, 'container_category': target_origin}
                    return_statement, previous_product_list = chemical_characterization_lookup(full_library_mongo,
                                                                                               target_reagent,
                                                                                               reagent_details)
                    if return_statement != 'Success':
                        problematic_wells.append(previous_product_list)
                        continue
                    additional_reagents_needed[well_to_clean]['initial_well_details'] = deepcopy(full_library_mongo[
                         previous_product_list[0]][previous_product_list[1]][previous_product_list[2]][
                         previous_product_list[3]])
                    additional_reagents_needed[well_to_clean]['container_name'] = deepcopy(full_library_mongo[
                         previous_product_list[0]][previous_product_list[1]]['container_name'])
                    additional_reagents_needed[well_to_clean]['container_id'] = previous_product_list[1]
                    additional_reagents_needed[well_to_clean]['collection'] = previous_product_list[0]
                    additional_reagents_needed[well_to_clean]['origin_well_location'] = previous_product_list[3]
                current_reagent_amount = current_reagents[0][1]
                current_reagent_concentration = current_reagent_amount / (current_well_volume * 1E-6)
                print(current_reagent_amount, current_well_volume)
                target_reagent_concentration = current_reagent_concentration * concentration_factor
                print(well_to_clean, current_reagent_concentration, target_reagent_concentration, concentration_factor)
                # Find the concentration of the stock solution in the original wellplate
                if additional_reagents_needed[well_to_clean]['container_category'] == 'reaction_plate':
                    previous_products = additional_reagents_needed[well_to_clean]['initial_well_details']['target_product']
                    relevant_product = []
                    for previous_product in previous_products:
                        if target_reagent in previous_product:
                            relevant_product = previous_product
                    if not relevant_product:
                        problematic_wells.append('Previous product %s, not found' % target_reagent)
                        continue
                    previous_mols = previous_product[1]
                    previous_volume = additional_reagents_needed[well_to_clean]['initial_well_details']['total_volume']
                    stock_solution_concentration = (previous_mols / (previous_volume * 1E-6))
                else:
                    stock_solution_concentration = additional_reagents_needed[well_to_clean]['initial_well_details'][
                        'concentration_molar']
                    print('STOCK:', stock_solution_concentration)
                if target_reagent_concentration >= stock_solution_concentration:
                    target_reagent_concentration = stock_solution_concentration
                    #problematic_wells.append('Well %s target concentration %s incompatible with stock concentration %s' %(
                    #                            well_to_clean, target_reagent_concentration, stock_solution_concentration))
                conc_func = lambda x: abs((1 / current_well_volume) * ((current_well_volume - x) * current_reagent_concentration +
                                                                   x * stock_solution_concentration) - target_reagent_concentration)
                cons = ({'type': 'ineq', 'fun': lambda x: x - 1},
                        {'type': 'ineq', 'fun': lambda x: current_well_volume - x})
                best_volume_opt = optimize.minimize(conc_func, [1.0], constraints = cons, method = "cobyla")
                volume_to_add = list(best_volume_opt['x'])[0]
                new_concentration = (1 / current_well_volume) * ((current_well_volume - volume_to_add) * current_reagent_concentration +
                                                                   volume_to_add * stock_solution_concentration)
                dilution_sequence = []
                current_concentration = new_concentration
                for dilution in range(0, 8):
                    if target_reagent_concentration > 0.5 * current_concentration:
                        dilution_factor = target_reagent_concentration / current_concentration
                        dilution_volume = (1 - dilution_factor) * current_well_volume
                        if dilution_volume > 1:
                            dilution_sequence.append((1 - dilution_factor) * current_well_volume)
                        break
                    else:
                        current_concentration = 0.5 * current_concentration
                        dilution_sequence.append(0.5 * current_well_volume)
                else:
                    problematic_wells.append('Problem with %s, dilution is too much %s' % (target_reagent, len(dilution_sequence)))
                print(dilution_sequence)
                print(current_concentration, target_reagent_concentration)
                additional_reagents_needed[well_to_clean]['concentrating_volume_to_add'] = volume_to_add
                additional_reagents_needed[well_to_clean]['dilution_sequence'] = dilution_sequence
                if target_reagent not in reagent_sequence.keys():
                    reagent_sequence[target_reagent] = {'container_name': additional_reagents_needed[well_to_clean]['container_name'],
                                                        'container_id': additional_reagents_needed[well_to_clean]['container_id'],
                                                        'collection': additional_reagents_needed[well_to_clean]['collection'],
                                                        'destinations': {'1': []}}
                reagent_sequence[target_reagent]['destinations']['1'].append([well_to_clean, volume_to_add, 1])
            unknown_cleanups = list(set(list(cleanup_methods.keys())).difference(set(known_cleanups)))
            if len(unknown_cleanups) != 0:
                problematic_wells.append('Problem with well %s, a specified cleanup method %s is not implemented' % (well_to_clean, unknown_cleanups))
        simplified_schedule, method_strings = [], []
        # First we need to make sure that we have the wellplate in an accessible location
        allowed_prep_carriers = [32, 38, 44]
        if wellplate_initial['location'][1][0] in allowed_prep_carriers:
            wellplate_home = wellplate_initial['location'][1]
        else:
            return_statement, wellplate_home = find_accessible_open_location(updated_carriers)
            if return_statement == 'Error':
                    return ['Error', wellplate_home, success_statements]
        success_statements.append('Found a suitable prep location at location: %s' % wellplate_home)

        if wellplate_initial['location'][1] != wellplate_home:
            simplified_schedule.append(['TRANSFER_LABWARE', wellplate_initial['container_name'],
                                        wellplate_initial['location'][1], wellplate_home,
                                        str(wellplate_initial['_id']), 'wellplates'])
            updated_carriers[wellplate_home[0]][wellplate_home[2]]['labware_labels'][
                wellplate_home[1] - 1] = 'location_%s_%s_%s' % (str(wellplate_home[0]),
                                                    str(wellplate_home[1]), str(wellplate_home[2]))
            updated_carriers[wellplate_home[0]][wellplate_home[2]]['labware_types'][
                wellplate_home[1] - 1] = wellplate_initial['labware_type']

        # Now the aspirate and dispense functions are added to the simplified schedule
        simplified_schedule.append(['TIP_RINSE', 10, 10])
        simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
        simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])

        liha_tip_number = 1
        acceptable_bed_locations = [32, 38, 44]
        location_occupancies = []
        solvent_data = {}
        if additional_dilutions_needed:
            for target_destination in additional_dilutions_needed.keys():
                dilution_sequence = additional_dilutions_needed[target_destination]['dilution_sequence']
                relevant_solvent = additional_dilutions_needed[target_destination]['solvents']
                if relevant_solvent not in solvent_data.keys():
                    solvent_stub = {'destination': [[target_destination, 100, 1]]}
                    return_statement, solvent_information_return = solvent_finder(full_library_mongo, relevant_solvent, solvent_stub)
                    if return_statement != 'Success':
                        problematic_wells.append(solvent_information_return)
                        continue
                    solvent_data[relevant_solvent] = solvent_information_return
                solvent_information = full_library_mongo[solvent_data[relevant_solvent]['collection']][solvent_data[relevant_solvent]['container_id']]
                solvent_location = solvent_information['location'][1]
                solvent_origin_well = solvent_data[relevant_solvent]['well_location']

                if liha_tip_number + len(dilution_sequence) > 8:
                    simplified_schedule.append(['TIP_RINSE', 5, 5])
                    simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
                    simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
                    liha_tip_number = 1
                    
                for seq_num, dilution_seq in enumerate(dilution_sequence):
                    tips = [liha_tip_number]
                    volumes = [dilution_seq]
                    wells = [target_destination]
                    cleanup_wells = [chr(ord('A') + liha_tip_number - 1) + '1']
                    simplified_schedule.append(['LiHa_ASPIRATE_DISPENSE', wellplate_home + [wells],
                                                        [1, 2, 30], tips, volumes, cleanup_wells])
                    if liha_tip_number >= 8:
                        simplified_schedule.append(['TIP_RINSE', 5, 5])
                        simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
                        simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
                        liha_tip_number = 1
                    else:
                        liha_tip_number += 1
                    
                    tips = [liha_tip_number]
                    simplified_schedule.append(['LiHa_ASPIRATE_DISPENSE', solvent_location + [[solvent_origin_well]],
                                                        wellplate_home, tips, volumes, wells])
                    simplified_schedule.append(['MIX_LIHA', wellplate_home, tips, wells, [200], 10])
                    if liha_tip_number >= 8:
                        simplified_schedule.append(['TIP_RINSE', 5, 5])
                        simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
                        simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
                        liha_tip_number = 1
                    else:
                        liha_tip_number += 1
                    
                    solvent_data[relevant_solvent]['final_solvent_details']['volume_ul'] -= dilution_seq
                    if solvent_data[relevant_solvent]['final_solvent_details']['volume_ul'] <= 0:
                        return ['Error', 'Not enough %s present on the platform' % relevant_solvent]
                
                if liha_tip_number > 8:
                    simplified_schedule.append(['TIP_RINSE', 5, 5])
                    simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
                    simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
                    liha_tip_number = 1
                else:
                    liha_tip_number += 1

            simplified_schedule.append(['TIP_RINSE', 5, 5])
            simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
            simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
            simplified_schedule.append(['TIP_RINSE', 4, 4])
            liha_tip_number = 1

            simplified_schedule.append(['MOVE_LiHa'])
        if additional_reagents_needed:
            # pprint(reagent_sequence)
            return_statement, reagents_order = priority_sort(reagent_sequence, 'previous_products')
            if return_statement == 'Error':
                return ['Error', reagents_order, success_statements]
            for reagent_seq in sorted(list(reagents_order.keys())):
                for target_container in reagents_order[reagent_seq].keys():
                    for collection in full_library_mongo.keys():
                        if target_container in full_library_mongo[collection].keys():
                            wellplate_definition = full_library_mongo[collection][target_container]
                            break
                    else:
                        return ['Error', 'Problem: Unable to find %s' % target_container, success_statements]
                    current_plate_location = wellplate_definition['location'][1]
                    success_statements.append('Currently the reagent tray "%s" is located at location: %s' % (
                                            wellplate_definition['container_name'], current_plate_location))
                    if wellplate_definition['location'][1][0] not in acceptable_bed_locations:
                        exceptions = [location[0] for location in location_occupancies if location[1] != wellplate_definition['labware_type']]
                        return_statement, new_wellplate_home = find_accessible_open_location(updated_carriers, exceptions,
                                                                                                  wellplate_definition['labware_type'])
                        if return_statement == 'Error':
                            return ['Error', 'Reagents well-plate %s not accessible' % wellplate_definition['container_name'],
                                    success_statements]

                        if wellplate_definition['location'] != new_wellplate_home:
                            simplified_schedule.append(['TRANSFER_LABWARE', wellplate_definition['container_name'],
                                                        wellplate_definition['location'][1], new_wellplate_home,
                                                        str(wellplate_definition['_id']), collection])
                            current_plate_location = new_wellplate_home
                    location_occupancies.append([current_plate_location, wellplate_definition['labware_type']])
                    for target_molecule in reagents_order[reagent_seq][target_container].keys():
                        for target_destination in reagents_order[reagent_seq][target_container][target_molecule]['destinations']:
                            target_well = target_destination[0]
                            relevant_solvent = additional_reagents_needed[target_well]['solvents']
                            # print(solvent_data)
                            if relevant_solvent not in solvent_data.keys():
                                solvent_stub = {'destination': [[target_well, 0, 1]]}
                                return_statement, solvent_information_return = solvent_finder(full_library_mongo, relevant_solvent, solvent_stub)
                                if return_statement != 'Success':
                                    problematic_wells.append(solvent_information_return)
                                    continue
                                solvent_data[relevant_solvent] = solvent_information_return
                            solvent_information = full_library_mongo[solvent_data[relevant_solvent]['collection']][solvent_data[relevant_solvent]['container_id']]
                            solvent_location = solvent_information['location'][1]
                            solvent_origin_well = solvent_data[relevant_solvent]['well_location']
                            first_step = additional_reagents_needed[target_well]['concentrating_volume_to_add']
                            dilution_sequence = additional_reagents_needed[target_well]['dilution_sequence']
                            # pprint(additional_reagents_needed)
                            if liha_tip_number + len(dilution_sequence) >= 8:
                                simplified_schedule.append(['TIP_RINSE', 5, 5])
                                simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
                                simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
                                liha_tip_number = 1
                            tips = [liha_tip_number]
                            volumes = [first_step]
                            wells = [target_well]
                            cleanup_wells = [chr(ord('A') + liha_tip_number - 1) + '1']
                            simplified_schedule.append(['LiHa_ASPIRATE_DISPENSE', wellplate_home + [wells],
                                                                [1, 2, 30], tips, volumes, cleanup_wells])
                            liha_tip_number += 1
                            tips = [liha_tip_number]
                            if liha_tip_number >= 8:
                                simplified_schedule.append(['TIP_RINSE', 5, 5])
                                simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
                                simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
                                liha_tip_number = 1
                            if 'well_location' in additional_reagents_needed[target_well].keys():
                                origin_well = [additional_reagents_needed[target_well]['well_location']]
                            else:
                                origin_well = [additional_reagents_needed[target_well]['origin_well_location']]
                            simplified_schedule.append(['LiHa_ASPIRATE_DISPENSE', current_plate_location + [origin_well],
                                                                wellplate_home + [wells], tips, volumes, wells])
                            simplified_schedule.append(['MIX_LIHA', wellplate_home, tips, wells, [200], 10])
                            liha_tip_number += 1
                            if liha_tip_number >= 8:
                                simplified_schedule.append(['TIP_RINSE', 5, 5])
                                simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
                                simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
                                liha_tip_number = 1
                            tips = [liha_tip_number]
                            if 'final_well_details' in additional_reagents_needed[target_well].keys():
                                additional_reagents_needed[target_well]['final_well_details']['volume_ul'] -= first_step
                            else:
                                previous_volume = additional_reagents_needed[target_well]['initial_well_details']['total_volume']
                                fraction_remaining = 1 - first_step / previous_volume
                                previous_volume -= first_step
                                previous_solvents = additional_reagents_needed[target_well]['initial_well_details']['solvents']
                                # pprint(additional_reagents_needed[target_well]['initial_well_details'])
                                if len(solvents) != 1:
                                    updated_solvents = [[item[0], item[1] * fraction_remaining] for item in previous_solvents]
                                    additional_reagents_needed[target_well]['initial_well_details']['solvents'] = updated_solvents
                                    previous_moles = additional_reagents_needed[target_well]['initial_well_details']['target_product']
                                    updated_moles = [[item[0], item[1] * fraction_remaining] for item in previous_moles]
                                    additional_reagents_needed[target_well]['initial_well_details']['target_product'] = updated_moles

                            number_of_dilutions = len(dilution_sequence)
                            print(dilution_sequence)
                            for seq_num, dilution_seq in enumerate(dilution_sequence):
                                cleanup_volumes = [dilution_seq]
                                cleanup_wells = [chr(ord('A') + liha_tip_number) + '1']
                                simplified_schedule.append(['LiHa_ASPIRATE_DISPENSE', wellplate_home + [[target_destination[0]]],
                                                            [1, 2, 30], tips, cleanup_volumes, cleanup_wells])
                                liha_tip_number += 1
                                if seq_num == 0 and liha_tip_number + number_of_dilutions > 8:
                                    simplified_schedule.append(['TIP_RINSE', 5, 5])
                                    simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
                                    simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
                                    liha_tip_number = 1
                                else:
                                    liha_tip_number += 1
                                
                                tips = [liha_tip_number]
                                volumes = [dilution_seq]
                                wells = [target_destination[0]]
                                simplified_schedule.append(['LiHa_ASPIRATE_DISPENSE', solvent_location + [[solvent_origin_well]],
                                                    wellplate_home, tips, volumes, wells])
                                solvent_data[relevant_solvent]['final_solvent_details']['volume_ul'] -= dilution_seq
                                simplified_schedule.append(['MIX_LIHA', wellplate_home, tips, wells, [200], 10])
                                liha_tip_number += 1
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
                        return ['Error', 'Cannot access the reagent tray %s for back transfer to its home' % wellplate_definition['container_name'], success_statements]
                    simplified_schedule.append(['TRANSFER_LABWARE', wellplate_definition['container_name'],
                                                current_plate_location, wellplate_definition['location'][1],
                                                str(wellplate_definition['_id']), collection])
                simplified_schedule.append(['TIP_RINSE', 5, 5])
                simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
                simplified_schedule.append(['SOLVENT_RINSE', 'DMSO_WASH', 300])
                simplified_schedule.append(['TIP_RINSE', 4, 4])
                liha_tip_number = 1

        if len(problematic_wells) != 0:
            return ['Error', problematic_wells, success_statements]
        if wellplate_home != wellplate_initial['location'][1]:
            problem_statement = check_for_problems(updated_carriers, wellplate_home,
                                                   wellplate_initial['container_name'])
            if problem_statement != 'Acceptable':
                return ['Error', 'Cannot access the reaction plate %s for back transfer to its home' %
                        wellplate_initial['container_name'], success_statements]
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

        # We need to convert the simplified schedule into Tecan method strings
        return_statement, method_strings, updated_carriers = convert_to_method_strings(simplified_schedule, updated_carriers, full_library_mongo, carriers_labware)
        if return_statement == 'Error':
            return ['Error', method_strings, success_statements]

        with open(additional_details['schedule_path_output'], 'w', encoding='latin1') as outfile:
            for line in filled_worktable_export:
                outfile.write(line + '\n')
            for line in method_strings:
                outfile.write(line + '\n')

        final_wellplate_contents = deepcopy(wellplate_initial['contents'])
        if additional_dilutions_needed:
            for well_key in additional_dilutions_needed.keys():
                diluted_reagents = []
                dilution_factor = additional_dilutions_needed[well_key]['dilution_factor']
                for reagent in final_wellplate_contents[well_key]['reagents']:
                    diluted_reagents.append([reagent[0], reagent[1] * dilution_factor, reagent[2]])
                diluted_solvents = []
                for solvent in final_wellplate_contents[well_key]['solvents']:
                    if solvent[0] == additional_dilutions_needed[well_key]['solvents']:
                        dilution_volumes = sum(additional_dilutions_needed[well_key]['dilution_sequence'])
                        new_volume = solvent[1] * dilution_factor + dilution_volumes
                    else:
                        new_volume = solvent[1] * dilution_factor
                    diluted_solvents.append([solvent[0], new_volume, solvent[2]])
                final_wellplate_contents[well_key]['solvents'] = deepcopy(diluted_solvents)
                final_wellplate_contents[well_key]['reagents'] = deepcopy(diluted_reagents)
        if additional_reagents_needed:
            for well_key in additional_reagents_needed.keys():
                concentrated_reagents = []
                concentration_factor = additional_reagents_needed[well_key]['concentration_factor']
                for reagent in final_wellplate_contents[well_key]['reagents']:
                    concentrated_reagents.append([reagent[0], reagent[1] * concentration_factor, reagent[2]])
                final_wellplate_contents[well_key]['reagents'] = concentrated_reagents
        full_library_mongo['wellplates'][wellplate_initial['_id']]['contents'] = final_wellplate_contents
        for solvent in solvent_data.keys():
            solvent_info = solvent_data[solvent]
            full_library_mongo[solvent_info['collection']][solvent_info['container_id']]['contents'] = solvent_info['final_solvent_details']

        for well_key in additional_reagents_needed.keys():
            reagent_information = additional_reagents_needed[well_key]
            if 'initial_well_details' in reagent_information.keys():
                keys_to_overwrite = ['solvents', 'target_product', 'total_volume']
                
                for key in keys_to_overwrite:
                    if key not in reagent_information['initial_well_details']:
                        continue
                    full_library_mongo[reagent_information['collection']][reagent_information['container_id']][
                        'contents'][reagent_information['origin_well_location']][key] = additional_reagents_needed[
                        well_key]['initial_well_details'][key]

        success_statements.append('The completed characterization well-plate is located at: %s' % wellplate_initial['location'][1])
        if __name__ == "__main__":
            # pass
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
    prep_return = cleanup_characterization_plate(queuename, '8')
    pprint(prep_return, width=250)
