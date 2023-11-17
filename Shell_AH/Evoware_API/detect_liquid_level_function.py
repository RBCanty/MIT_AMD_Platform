# -*- coding: utf-8 -*-
"""
Updated on 7/14/21 --- MongoDB containing version of detect liquid level function

@author: Brent
"""

import copy
import traceback
import pymongo
from pprint import pprint
from copy import deepcopy
from datetime import datetime
import sys

if __name__ == "__main__":
    from supporting_functions import get_queue_information, convert_to_method_strings
    from platform_library_functions import update_database_document, query_document
    from useful_other_functions import find_labware, evaporation_estimation, solvent_finder_origin
    from worktable_cleaner import initial_worktable_prep, worktable_export
else:
    from Evoware_API.supporting_functions import get_queue_information, convert_to_method_strings
    from Evoware_API.platform_library_functions import update_database_document, query_document
    from Evoware_API.useful_other_functions import find_labware, evaporation_estimation, solvent_finder_origin
    from Evoware_API.worktable_cleaner import initial_worktable_prep, worktable_export


def detect_liquid_level(queue_name, operation_number):
    success_statements = []
    try:
        # ----------------------------------------------------------------------------
        # First we start with opening the queue document and getting operation info.

        queue_information = get_queue_information(queue_name, operation_number, 'detect_liquid_level')
        if queue_information[0] != 'Success':
            return ['Error', 'Problem: %s' % queue_information[1], success_statements]
        else:
            return_statement, full_library_mongo, queue_document, additional_details, statements = queue_information

        success_statements.extend(statements)
        additional_prep_details = additional_details['additional_prep_details']

        # Then build the current worktable dictionary for building output files, we are
        # also making a duplicate of the worktable, so it can be updated in this code
        return_statement, start_carriers, carriers_labware = initial_worktable_prep(
            additional_details['carriers_labware_filepath'],
            additional_details['empty_worktable_filepath'],
            full_library_mongo)
        if return_statement != 'Success':
            return ['Error', 'Problem: %s' % start_carriers, success_statements]
        updated_carriers = copy.deepcopy(start_carriers)
        success_statements.append('Loaded the carrier and labware library file!')
        # ----------------------------------------------------------------------------

        # Start by finding where the wellplate is located on the platform, if it is on
        # a heater-shaker we need to make sure it is not shaking, otherwise the plate
        # should be ready for detection
        return_statement, labware_details = find_labware(updated_carriers, full_library_mongo,
                                                         additional_details['wellplate_name'],
                                                         additional_details['target_well_plate_type'],
                                                         additional_details['target_container'])
        if return_statement != 'Success':
            return ['Error', 'Problem: %s' % labware_details, success_statements]
        initial_labware_details = deepcopy(labware_details)
        wellplate_location = labware_details['location'][1]
        labware_collection = labware_details['collection']
        labware_container_id = labware_details['_id']

        # Check to make sure that the target detection well has been specified
        detection_well = additional_prep_details['detection_well']
        wellplate_name = additional_details['wellplate_name']

        simplified_schedule = []
        current_time = datetime.now()
        date_updated = current_time.strftime('%m/%d/%Y %H:%M:%S')
        heater_shaker_carrier_ids = [322, 323, 324]
        if wellplate_location[2] in heater_shaker_carrier_ids:
            simplified_schedule.append(['STOP_SHAKER', wellplate_location, 'off', 25])
            success_statements.append('The well-plate at %s has been set to stop shaking' % wellplate_location)

        # Check to see if the volume_history key is in the details section
        simplified_schedule.append(['TIP_RINSE', 5, 5])
        if 'volume_history' not in additional_prep_details.keys():
            if detection_well in labware_details['contents'].keys():
                if labware_details['contents'][detection_well]['confirmation'] != 'reference_well':
                    return ['Error', 'Problem: Well %s not ready for reference (%s)' % (detection_well, wellplate_name),
                            success_statements]
            else:
                success_statements.append('Wellplate %s has well %s ready' % (additional_details['wellplate_name'],
                                                                              detection_well))
                # Look through all the wells to find the one to replicate
                volume_to_replicate = 0
                for well_key in labware_details['contents'].keys():
                    well_details = labware_details['contents'][well_key]
                    well_volume = 0
                    for solvent in well_details['solvents']:
                        if 'water' in solvent[0]:
                            well_volume += solvent[1]
                    if well_volume > volume_to_replicate:
                        volume_to_replicate = well_volume

                if volume_to_replicate < 50:
                    return ['Error', 'Problem: Volume to replicate too small %s' % volume_to_replicate,
                            success_statements]
                solvent_to_find = 'water'
                solvent_dict_stub = [[solvent_to_find, volume_to_replicate]]
                return_statements, solvent_information = solvent_finder_origin(full_library_mongo,
                                                                               solvent_dict_stub, '')
                if return_statements == 'Error':
                    return ['Error', 'Problem: %s' % solvent_information, success_statements]
                collection = solvent_information[solvent_to_find]['collection']
                container_id = solvent_information[solvent_to_find]['container_id']
                solvent_definition = full_library_mongo[collection][container_id]
                solvent_origin_well = solvent_information[solvent_to_find]['origin_well_location']
                liha_tip_number = 1

                number_of_transfers = volume_to_replicate//250
                for transfer_number in range(0, int(number_of_transfers)):
                    solvent_location = solvent_definition['location'][1]
                    simplified_schedule.append(['LiHa_ASPIRATE_DISPENSE', solvent_location + [[solvent_origin_well]],
                                                wellplate_location, [liha_tip_number],
                                                [volume_to_replicate/number_of_transfers], [detection_well]])

                relevant_well = full_library_mongo[collection][container_id]['contents'][solvent_origin_well]
                relevant_well['volume_ul'] -= volume_to_replicate
                if solvent_information[solvent_to_find]['solvent_volume'] - volume_to_replicate < 5000:
                    return ['Error', 'Problem: Not enough "%s" solvent present on the liquid handler' % solvent_to_find,
                            success_statements]

                new_well_contents = {'confirmation': 'reference_well', 'final_product': [], 'reagents': [],
                                     'solvents': solvent_dict_stub, 'target_product': [],
                                     'total_volume': volume_to_replicate, 'date_updated': date_updated}
                full_library_mongo[labware_collection][labware_container_id]['contents'][detection_well] = new_well_contents

        variable_name = 'DETECTED_VOLUME_1'
        simplified_schedule.append(['SET_VARIABLE', variable_name, 0, 0, 2])
        simplified_schedule.append(['LIHA_DETECT_LIQUID', [wellplate_location, [1], [detection_well]]])
        simplified_schedule.append(['EXPORT_VARIABLE', [variable_name], additional_details['liquid_detection_path']])
        simplified_schedule.append(['MOVE_LiHa'])

        if wellplate_location[2] in heater_shaker_carrier_ids:
            simplified_schedule.append(['START_SHAKER', wellplate_location, 600, 100])
            success_statements.append('The well-plate at %s has been set to stop shaking' % wellplate_location)

        # ----------------------------------------------------------------------------
        # Now we need to rebuild the worktable, so we can write a file to run
        return_statement, filled_worktable_export = worktable_export(additional_details['carriers_labware_filepath'],
                                                                     additional_details['empty_worktable_filepath'],
                                                                     updated_carriers, full_library_mongo,
                                                                     simplified_schedule)
        if return_statement == 'Error':
            return ['Error', 'Problem: %s' % filled_worktable_export, success_statements]

        # We need to convert the simplified schedule into Tecan method strings
        return_statement, method_strings, updated_carriers = convert_to_method_strings(simplified_schedule,
                                                                                       updated_carriers,
                                                                                       full_library_mongo,
                                                                                       carriers_labware)
        if return_statement == 'Error':
            return ['Error', 'Problem: %s' % method_strings, success_statements]

        # Now we need to write the worktable and method string to a file, this
        # requires a specific type of encoding that includes extended ASCII
        # characters, latin1 is pretty close to it
        with open(additional_details['schedule_path_output'], 'w', encoding='latin1') as outfile:
            for line in filled_worktable_export:
                outfile.write(line + '\n')
            for line in method_strings:
                outfile.write(line + '\n')

        pprint(method_strings, width=200)
        if __name__ == "__main__":
            pass
            return ['Testing!!', success_statements]

        platform_mongodb = pymongo.MongoClient(additional_details['platform_mongodb_address'],
                                               additional_details['platform_mongodb_port'])
        if __name__ == "__main__":
            updated_queue_document = copy.deepcopy(queue_document)
            updated_queue_document['operations'][operation_number]['completed'] = 'yes'
            flag_return, return_statements = update_database_document(queue_document, updated_queue_document,
                                                                      'queue', platform_mongodb)
            if flag_return == 'Error':
                platform_mongodb.close()
                return ['Error', 'Issue with the database queue document update: %s' % return_statements,
                        success_statements]
        flag_return, return_statements = update_database_document(initial_labware_details,
                                                           full_library_mongo[labware_collection][labware_container_id],
                                                           'wellplates', platform_mongodb)
        if flag_return == 'Error':
            platform_mongodb.close()
            return ['Error', 'Issue with updating %s in %s' % (labware_container_id, labware_collection),
                    success_statements]
        platform_mongodb.close()
        success_statements.append('Updated the platform library to be up-to-date following this run')
        detect_details = {'wellplate_name': additional_details['wellplate_name'],
                          'target_container': additional_details['target_container'],
                          'detection_well': detection_well}
        return ['Success', additional_details['schedule_path_output'], success_statements]

    except Exception:
        return ['Error', 'Problem: %s' % traceback.format_exc(), success_statements]


if __name__ == "__main__":
    queuename = 'sulfur_explore_scaffold_round2_4_20220804'
    prep_return = detect_liquid_level(queuename, '18')
    pprint(prep_return, width = 200)