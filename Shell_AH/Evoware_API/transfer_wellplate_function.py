# -*- coding: utf-8 -*-
"""
Updated on 7/14/21 --- MongoDB containing version of plate transfer 

@author: Brent
"""

import copy
import traceback
import pymongo
from pprint import pprint
import sys

if __name__ == "__main__":
    from supporting_functions import get_queue_information, convert_to_method_strings
    from platform_library_functions import update_database, update_database_document
    from useful_other_functions import find_accessible_open_location, find_labware, \
        check_for_problems, find_open_heater_shaker, find_accessible_hotel_location, check_spark_availability
    from worktable_cleaner import initial_worktable_prep, worktable_export
else:
    from Evoware_API.supporting_functions import get_queue_information, convert_to_method_strings
    from Evoware_API.platform_library_functions import update_database, update_database_document
    from Evoware_API.useful_other_functions import find_accessible_open_location, find_labware, \
        check_for_problems, find_open_heater_shaker, find_accessible_hotel_location, check_spark_availability
    from Evoware_API.worktable_cleaner import initial_worktable_prep, worktable_export


# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# These locations need to be updated with adjustments to the worktable
platform_locations = {'low_temperature_prep_locations': [[11, 1, 322], [11, 2, 322]],
                      'thermoshake_locations': [[11, 1, 322], [11, 2, 322], [61, 1, 324], [61, 2, 324]],
                      'teleshake_locations': [[25, 1, 323], [25, 2, 323]],

                      'standard_transfer_locations': [[24, 1, 329]],
                      'special_transfer_locations': [[19, 5, 84], [19, 6, 84], [24, 1, 329]],

                      'inert_box_block_hotel': [6, 1, 335],
                      'inert_box_shaker_location': [25, 4, 323],
                      'small_inert_box_block_hotel': [7, 1, 344]}
# TODO Add the small box to the platform locations lookup here: [7, 1, 344]
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------


def transfer_wellplate(queue_name, operation_number):
    success_statements = []
    try:
        # ----------------------------------------------------------------------------
        # First we start with opening the queue document and getting operation info.
        queue_information = get_queue_information(queue_name, operation_number, 'transfer_wellplate')
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
        simplified_schedule = []
        method_strings = []
        necessary_transfers = []
        if additional_prep_details['target_destination'] == 'transfer_prep':
            if additional_details['target_well_plate_type'] == '96 Well Microplate':
                relevant_locations = platform_locations['special_transfer_locations']
            else:
                relevant_locations = platform_locations['standard_transfer_locations']
            # By looking at each of the locations, we can see if we need to move
            # wellplate off of the transfer hotel positions
            occupied_locations = []
            for relevant_location in relevant_locations:
                carrier_location = updated_carriers[relevant_location[0]][relevant_location[2]]
                if carrier_location['labware_types'][relevant_location[1] - 1] != '':
                    labware_label = carrier_location['labware_labels'][relevant_location[1] - 1]
                    occupied_locations.append([relevant_location, labware_label])
            if len(occupied_locations) == 0:
                return ['Success', 'Nothing to move', success_statements]
            if additional_details['target_well_plate_type'] == '96 Well Microplate' and len(occupied_locations) == 1:
                return ['Success', 'Nothing to move', success_statements]

            # Now check to see if there is a place to move those wellplates to
            bed_exceptions = []
            for occupied_location in occupied_locations:
                # We start by finding where the wellplate we want to move is located
                # TODO: This section of code does not work as intended, it should be moving
                # wellplates that occupy useful locations, it does not look up correctly
                # when there is a plate that needs to be moved
                return_statement, labware_details = find_labware(updated_carriers, full_library_mongo,
                                                                 additional_details['wellplate_name'],
                                                                 additional_details['target_well_plate_type'])
                if return_statement != 'Success':
                    return ['Error', labware_details, success_statements]
                initial_wellplate_definition = copy.deepcopy(labware_details)
                initial_location = initial_wellplate_definition['location'][1]

                # Then look for an open location to move the wellplate to
                return_statement, target_location = find_accessible_open_location(updated_carriers, exceptions=bed_exceptions)
                if return_statement == 'Error':
                    return ['Error', target_location, success_statements]
                simplified_schedule.append(['TRANSFER_LABWARE', occupied_location[1], occupied_location[0],
                                           target_location, str(labware_details['_id']), labware_details['collection']])
                bed_exceptions.append(target_location)
                necessary_transfers.append(occupied_location + [target_location])

                # We will update the carrier dictionary for this change
                updated_carriers[target_location[0]][target_location[2]]['labware_labels'][target_location[1] - 1] = \
                    'location_%s_%s_%s' % (str(target_location[0]), str(target_location[1]), str(target_location[2]))
                updated_carriers[target_location[0]][target_location[2]]['labware_types'][target_location[1] - 1] = \
                    labware_details['labware_type']
                full_library_mongo[labware_details['collection']][str(labware_details['_id'])]['location'][1] = target_location

                # If it is a 96 Well Microplate transfer sequence only one wellplate needs to be moved
                if additional_details['target_well_plate_type'] == '96 Well Microplate':
                    break

        elif additional_prep_details['target_destination'] == 'transfer_cleanup':
            # Check if this is brining in new MCA tips
            if additional_details['target_well_plate_type'] == 'DiTi SBS Waste':
                # First find the empty tip boxes on the liquid handler
                empty_tip_boxes = []
                nested_diti_key = None
                for consumable_key in full_library_mongo['consumables'].keys():
                    if full_library_mongo['consumables'][consumable_key]['container_name'] == 'DiTi_200uL_Nested_328':
                        tip_locations = full_library_mongo['consumables'][consumable_key]['tip_locations']
                        nested_diti_key = consumable_key
                        if '15' not in tip_locations:
                            empty_tip_boxes.append(25)
                        if '16' not in tip_locations:
                            empty_tip_boxes.append(26)
                        break
                if nested_diti_key is None:
                    return ['Error', 'Problem: Nested tips location not found in database', success_statements]
                # Then find out if those two base locations are occupied
                nested_location = []
                for grid in updated_carriers.keys():
                    for carrier in updated_carriers[grid].keys():
                        if carrier == 328:
                            nested_location = [grid, carrier]
                            break
                    else:
                        continue
                    break
                else:
                    return ['Error', 'Problem: Nested tips location not found', success_statements]
                occupied_locations = full_library_mongo['consumables'][nested_diti_key]['occupied_locations']
                transfers_needed = []
                for occupied_location in occupied_locations:
                    if occupied_location in empty_tip_boxes:
                        labware_type = updated_carriers[nested_location[0]][nested_location[1]]['labware_types'][occupied_location - 1]
                        transfers_needed.append([[nested_location[0], occupied_location, nested_location[1]], labware_type])
                success_statements.append('Found %s empty tip boxes to deal with' % len(transfers_needed))
                # Move those empty tip boxes to the waste chute
                for transfer_needed in transfers_needed:
                    simplified_schedule.append(['LARGE_WASTE_DUMP', transfer_needed])
                    full_library_mongo['consumables'][nested_diti_key]['occupied_locations'].remove(transfer_needed[0][1])
                success_statements.append('Needed to do %s transfers to get rid of tip boxes' % len(transfers_needed))
                # Find any new tips that are present on the platform
                new_tip_locations = []
                for grid in updated_carriers.keys():
                    for carrier in updated_carriers[grid].keys():
                        if carrier == 328:
                            continue
                        for site_index, labware_type in enumerate(updated_carriers[grid][carrier]['labware_types']):
                            if labware_type in ['DiTi SBS Waste']:
                                labware_label = updated_carriers[grid][carrier]['labware_labels'][site_index]
                                for consumable_key in full_library_mongo['consumables'].keys():
                                    if full_library_mongo['consumables'][consumable_key]['container_name'] == labware_label:
                                        mongo_key = consumable_key
                                        break
                                else:
                                    return ['Error', 'Problem: Could not find the %s wellplate in the database' % labware_label, success_statements]
                                new_tip_locations.append([[grid, site_index+1, carrier], labware_label, mongo_key, labware_type])
                # Move the new tips to the nested tip carrier
                appropriate_sites = [25, 26]
                completed_transfers = []
                incomplete_transfers = []
                for new_tip_location in new_tip_locations:
                    for appropriate_site in appropriate_sites:
                        if appropriate_site not in full_library_mongo['consumables'][nested_diti_key]['occupied_locations']:
                            target_site = [nested_location[0], appropriate_site, nested_location[1]]
                            break
                    else:
                        incomplete_transfers.append([new_tip_location[0], new_tip_location[1], new_tip_location[2]])
                        continue

                    # Here we will first check to make sure that the new tip container is accessible
                    problem_return = check_for_problems(updated_carriers, new_tip_location[0], new_tip_location[3])
                    if problem_return != 'Acceptable':
                        return ['Error', 'Problem: Unable to access container %s to move' % new_tip_location[1], success_statements]

                    # Then we quickly check to make sure that we can access the target destination, we need to
                    # do this because we have not checked to make sure that the target destination is accessible
                    problem_return = check_for_problems(updated_carriers, target_site, new_tip_location[3])
                    if problem_return != 'Acceptable':
                        return ['Error', 'Problem: Unable to access site %s to move a container' % target_site, success_statements]
                    simplified_schedule.append(['TRANSFER_LABWARE', new_tip_location[1], new_tip_location[0],
                                               target_site, new_tip_location[2],'consumables'])
                    completed_transfers.append([new_tip_location[0], new_tip_location[1], new_tip_location[2], target_site])
                    consumable_definition = full_library_mongo['consumables'][new_tip_location[2]]
                    full_library_mongo['consumables'][nested_diti_key]['occupied_locations'].append(appropriate_site)

                    # Consolidate the tips into the consumable document
                    existing_tips = [int(item) for item in full_library_mongo['consumables'][nested_diti_key]['tip_locations']]
                    for index in range(0, len(consumable_definition['tip_locations'])):
                        existing_tips.append(target_site[1] - 10 - 2*index)
                    existing_tips.sort()
                    full_library_mongo['consumables'][nested_diti_key]['tip_locations'] = list(map(lambda x: str(x), existing_tips))
                if len(completed_transfers) < len(new_tip_locations):
                    success_statements.append('Was not able to transfer all tips: %s of %s completed' % (len(completed_transfers),
                                                                                                         len(new_tip_locations)))
                else:
                    success_statements.append('Completed all %s new tip rack transfers' % len(completed_transfers))

            else:
                # This happens for all other labwares on the platform for transfer_cleanup
                return_statement, target_location = find_accessible_open_location(updated_carriers)
                if return_statement != 'Success':
                    return ['Error', target_location, success_statements]
                return_statement, labware_details = find_labware(updated_carriers, full_library_mongo,
                                                                 additional_details['wellplate_name'],
                                                                 additional_details['target_well_plate_type'])
                initial_wellplate_definition = copy.deepcopy(labware_details)
                initial_location = initial_wellplate_definition['location'][1]
                simplified_schedule.append(['TRANSFER_LABWARE', additional_details['wellplate_name'], initial_location,
                                            target_location, str(labware_details['_id']), labware_details['collection']])
                necessary_transfers.append([initial_location, additional_details['wellplate_name'], target_location])
                full_library_mongo[labware_details['collection']][str(labware_details['_id'])]['location'][1] = target_location
            if 'transfers' not in additional_prep_details.keys():
                additional_transfers = []
            else:
                additional_transfers = additional_prep_details['transfers']
            for additional_transfer in additional_transfers:
                if 'DiTi SBS Waste' in additional_transfer[1]:
                    continue
                return_statement, labware_details = find_labware(updated_carriers, full_library_mongo,
                                                                 additional_transfer[1], 'name_only_lookup')
                if return_statement != 'Success':
                    return ['Error', 'Problem: Did not find %s to back-transfer' % additional_transfer[1], success_statements]

                original_site = additional_transfer[0]
                problem_return = check_for_problems(updated_carriers, original_site, labware_details['labware_type'])
                if problem_return != 'Acceptable':
                    return ['Error', 'Problem: Unable to access %s to move %s back' % (original_site, additional_transfer[1]), success_statements]

                current_site = additional_transfer[2]
                problem_return = check_for_problems(updated_carriers, current_site, labware_details['labware_type'])
                if problem_return != 'Acceptable':
                    return ['Error', 'Problem: Unable to access %s to move %s back' % (current_site, additional_transfer[1]), success_statements]

                simplified_schedule.append(['TRANSFER_LABWARE', additional_transfer[1], current_site,
                                           original_site, str(labware_details['_id']), labware_details['collection']])
                full_library_mongo[labware_details['collection']][str(labware_details['_id'])]['location'][1] = original_site

        else:
            # These are the standard transfers that are possible at the moment,
            # start by finding where the wellplate we want to move is located
            return_statement, labware_details = find_labware(updated_carriers, full_library_mongo,
                                                             additional_details['wellplate_name'], 'name_only_lookup')
            if return_statement != 'Success':
                return ['Error', labware_details, success_statements]
            initial_wellplate_definition = copy.deepcopy(labware_details)
            initial_location = initial_wellplate_definition['location'][1]

            # Then check to see if there are any transfer issues for the location
            problem_return = check_for_problems(updated_carriers, initial_location,
                                                initial_wellplate_definition['labware_type'])
            if problem_return != 'Acceptable':
                return ['Error', 'Problem: Unable to access location %s for transfers' % initial_location, success_statements]

            # Depending on the target destination we need to find the appropriate
            # postions to transfer the wellplate to
            target_destination = additional_prep_details['target_destination']
            if target_destination == 'heater_shaker':
                if 'temperature' not in additional_prep_details.keys():
                    return ['Error', 'Problem: No temperature specified for %s, %s' % (queue_name, operation_number), success_statements]
                target_temperature = additional_prep_details['temperature']
                return_statement, target_location = find_open_heater_shaker(updated_carriers, target_temperature,
                                                                            labware_details['labware_type'],
                                                                            platform_locations)
            elif target_destination == 'transfer_hotel':
                return_statement, target_location = find_accessible_hotel_location(updated_carriers, 'transfer',
                                                                                   labware_details['labware_type'],
                                                                                   platform_locations)
            elif target_destination == 'storage_hotel':
                return_statement, target_location = find_accessible_hotel_location(updated_carriers, 'storage',
                                                                                   labware_details['labware_type'],
                                                                                   platform_locations)
            elif target_destination == 'bed_position':
                return_statement, target_location = find_accessible_open_location(updated_carriers)

            elif target_destination == 'spark':
                return_statement, target_location = check_spark_availability(updated_carriers)

            else:
                return ['Error', 'Problem: Destination not implemented in the transfer script', success_statements]

            # Now check to make sure that a location was actually returned
            if return_statement == 'Error':
                return ['Error', target_location, success_statements]
            else:
                success_statements.append('Found a valid destination to transfer to: %s' % target_location)

            # There is an optional statement to disassemble the inert reactor block
            if 'disassemble_inert' in additional_prep_details.keys() and additional_prep_details['disassemble_inert'] == True:
                # Start by looking for the inert box block in the database
                if initial_wellplate_definition['labware_type'] in ['96 Well Microplate', '96 Well MediumWell']:
                    inert_box_type = 'Small Inert Box Block'
                    inert_box_block_hotel = platform_locations['small_inert_box_block_hotel']
                else:
                    inert_box_type = 'Inert Box Block'
                    inert_box_block_hotel = platform_locations['inert_box_block_hotel']
                return_statement, inert_box_details = find_labware(updated_carriers, full_library_mongo, None, inert_box_type)
                
                # return_statement, inert_box_details = find_labware(updated_carriers, full_library_mongo, None, 'Inert Box Block')
                if return_statement != 'Success':
                    return ['Error', inert_box_details, success_statements]
                # The locations for the inert box are hard-defined due to site requirements
                inert_box_block_hotel = platform_locations['inert_box_block_hotel']
                inert_box_shaker_location = platform_locations['inert_box_shaker_location']
                # We will check to see if the inert box is ready to be transferred
                # however, we will return an error if it is not
                if inert_box_details['location'][1] != inert_box_shaker_location:
                    return ['Error', 'Problem: Inert reaction block not in %s, instead it is at %s' % (inert_box_shaker_location,
                                inert_box_details['location'][1])]
                simplified_schedule.append(['MCA_MOVE', [8, 2, 48]])
                simplified_schedule.append(['SPECIAL_TRANSFER_LABWARE', 'Inert Box Block',
                                                inert_box_shaker_location, inert_box_block_hotel, '', ''])
                full_library_mongo['consumables'][str(inert_box_details['_id'])]['location'][1] = inert_box_block_hotel

            # Now we add the operation to move the target well-plate and update
            # the database for this movement
            if initial_location != target_location:
                simplified_schedule.append(['MCA_MOVE', [8, 2, 48]])
                simplified_schedule.append(['TRANSFER_LABWARE', additional_details['wellplate_name'], initial_location,
                                            target_location, str(labware_details['_id']), labware_details['collection']])
                full_library_mongo[labware_details['collection']][labware_details['_id']]['location'][1] = target_location

        if len(simplified_schedule) == 0:
            return ['Success', 'Nothing to move', success_statements]

        # ----------------------------------------------------------------------------
        # Now we need to rebuild the worktable so we can write a file to run
        return_statement, filled_worktable_export = worktable_export(additional_details['carriers_labware_filepath'],
                                                                     additional_details['empty_worktable_filepath'],
                                                                     updated_carriers, full_library_mongo,
                                                                     simplified_schedule)
        if return_statement == 'Error':
            return ['Error', filled_worktable_export, success_statements]

        # We need to convert the simplified schedule into Tecan method strings
        return_statement, method_strings, updated_carriers = convert_to_method_strings(simplified_schedule, updated_carriers,
                                                                                       full_library_mongo, carriers_labware)
        if return_statement == 'Error':
            return ['Error', method_strings, success_statements]

        with open(additional_details['schedule_path_output'], 'w', encoding='latin1') as outfile:
            for line in filled_worktable_export:
                outfile.write(line + '\n')
            for line in method_strings:
                outfile.write(line + '\n')

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
        if additional_prep_details['target_destination'] == 'transfer_prep':
            cleanup_operation_number = str(int(operation_number) + 1)
            updated_queue_document['operations'][cleanup_operation_number]['details']['transfers'] = necessary_transfers
        flag_return, return_statements = update_database_document(queue_document, updated_queue_document, 'queue', platform_mongodb)
        platform_mongodb.close()
        success_statements.append('Updated the platform library to be up-to-date following this run')

    except Exception:
        return ['Error', 'Problem: %s' % traceback.format_exc(), success_statements]

    return ['Success', additional_details['schedule_path_output'], success_statements]


if __name__ == "__main__":
    queuename = 'second_sulfur_explore_extra_20221019'
    prep_return = transfer_wellplate(queuename, '14')
    pprint(prep_return, width=150)
