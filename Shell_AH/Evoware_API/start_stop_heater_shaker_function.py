# -*- coding: utf-8 -*-
"""
Updated on 7/14/21 --- MongoDB containing version of heater shaker operations

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
    from useful_other_functions import find_labware, evaporation_estimation
    from worktable_cleaner import initial_worktable_prep, worktable_export
else:
    from Evoware_API.supporting_functions import get_queue_information, convert_to_method_strings
    from Evoware_API.platform_library_functions import update_database_document, query_document
    from Evoware_API.useful_other_functions import find_labware, evaporation_estimation
    from Evoware_API.worktable_cleaner import initial_worktable_prep, worktable_export


def start_stop_heater_shaker(queue_name, operation_number):
    success_statements = []
    try:
        # ----------------------------------------------------------------------------
        # First we start with opening the queue document and getting operation info.

        queue_information = get_queue_information(queue_name, operation_number, 'start_stop_heater_shaker')
        if queue_information[0] != 'Success':
            return ['Error', 'Problem: %s' % queue_information[1], success_statements]
        else:
            return_statement, full_library_mongo, queue_document, additional_details, statements = queue_information

        success_statements.extend(statements)
        additional_prep_details = additional_details['additional_prep_details']

        # Then build the current worktable dictionary for building output files, we are
        # also making a duplicate of the worktable so it can be updated in this code
        return_statement, start_carriers, carriers_labware = initial_worktable_prep(
            additional_details['carriers_labware_filepath'],
            additional_details['empty_worktable_filepath'],
            full_library_mongo)
        if return_statement != 'Success':
            return ['Error', 'Problem: %s' % start_carriers, success_statements]
        updated_carriers = copy.deepcopy(start_carriers)
        success_statements.append('Loaded the carrier and labware library file!')
        # ----------------------------------------------------------------------------

        # Start by finding where the wellplate is located on the platform, if
        # it is not on a heater shaker then return an Error, there has either
        # been database mismanagament or an incorrect operations sequence
        return_statement, labware_details = find_labware(updated_carriers, full_library_mongo,
                                                         additional_details['wellplate_name'],
                                                         additional_details['target_well_plate_type'])
        if return_statement != 'Success':
            return ['Error', 'Problem: %s' % labware_details, success_statements]
        wellplate_location = labware_details['location'][1]

        acceptable_heater_shaker_carrier_ids = [322, 323, 324]
        if wellplate_location[2] not in acceptable_heater_shaker_carrier_ids:
            return ['Error',
                    'Problem: The wellplate %s is not in a heater shaker' % additional_details['wellplate_name'],
                    success_statements]

        # The control options for the heater shakers are quite restricted, these
        # come from the details provided in the queue document, heating and shaking
        # are controlled by two separate functions
        simplified_schedule = []
        target_power = additional_prep_details['power']
        updated_database_document = {}
        if target_power == 'on':
            target_temperature = additional_prep_details['temperature']
            target_rpms = additional_prep_details['rpms']
            simplified_schedule.append(['SHAKER_TEMPERATURE', wellplate_location, target_power, target_temperature])
            success_statements.append(
                'The temperature has been set to %s at %s' % (target_temperature, wellplate_location))
            if target_rpms > 0:
                simplified_schedule.append(['START_SHAKER', wellplate_location, target_rpms])

            updated_database_document = deepcopy(labware_details)
            current_time = datetime.now()
            for well_key in updated_database_document['contents'].keys():
                updated_database_document['contents'][well_key]['date_updated'] = current_time.strftime(
                    '%m/%d/%Y %H:%M:%S')
                updated_database_document['contents'][well_key]['temperature_setpoint'] = target_temperature
            updated_database_document['date_updated'] = current_time.strftime('%m/%d/%Y')

        elif target_power == 'off':
            target_temperature = 25
            simplified_schedule.append(['SHAKER_TEMPERATURE', wellplate_location, target_power, target_temperature])
            success_statements.append('The temperature has been turned off for location %s' % wellplate_location)
            simplified_schedule.append(['STOP_SHAKER', wellplate_location, target_power, target_temperature])
            success_statements.append('The well-plate at %s has been set to stop shaking' % wellplate_location)

            updated_database_document = deepcopy(labware_details)
            previous_operation = str(int(operation_number)- 1)
            if queue_document['operations'][previous_operation]['operation'] == 'detect_liquid_level' \
                and updated_database_document['contents']['H12']['confirmation'] == 'reference_well':
                current_time = datetime.now()
                for well_key in updated_database_document['contents'].keys():
                    well_information = updated_database_document['contents'][well_key]
                    updated_database_document['contents'][well_key]['solvents'] = []
                    well_information['date_updated'] = current_time.strftime('%m/%d/%Y %H:%M:%S')
                updated_database_document['contents'].pop('H12', None)
            else:
                current_time = datetime.now()
                for well_key in updated_database_document['contents'].keys():
                    well_information = updated_database_document['contents'][well_key]
                    labware_type = updated_database_document['labware_type']
                    previous_time = datetime.strptime(well_information['date_updated'], '%m/%d/%Y %H:%M:%S')
                    elapsed_time = (current_time - previous_time).total_seconds()
                    temperature_setpoint = well_information['temperature_setpoint']
                    return_statement, solvent_return = evaporation_estimation(well_information['solvents'], labware_type,
                                                                              elapsed_time, temperature_setpoint)
                    if return_statement != 'Success':
                        return ['Error', 'Problem: %s' % solvent_return, success_statements]
                    updated_database_document['contents'][well_key]['solvents'] = solvent_return
                    well_information['date_updated'] = current_time.strftime('%m/%d/%Y %H:%M:%S')

        else:
            return ['Error', 'Problem: Input format not understood for heater-shaker operation', success_statements]

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

        if __name__ == "__main__":
            pass
            return ['Testing!!', success_statements]
        
        if updated_database_document:
            platform_mongodb = pymongo.MongoClient(additional_details['platform_mongodb_address'],
                                                   additional_details['platform_mongodb_port'])
            flag_return, return_statements = update_database_document(labware_details, updated_database_document,
                                                                      'wellplates', platform_mongodb)
            platform_mongodb.close()
            success_statements.append('Updated the platform library to be up-to-date following this run')
        
        if __name__ == "__main__":
            platform_mongodb = pymongo.MongoClient(additional_details['platform_mongodb_address'],
                                                   additional_details['platform_mongodb_port'])
            updated_queue_document = copy.deepcopy(queue_document)
            updated_queue_document['operations'][operation_number]['completed'] = 'yes'
            flag_return, return_statements = update_database_document(queue_document, updated_queue_document, 'queue',
                                                                      platform_mongodb)
            platform_mongodb.close()
        success_statements.append('Updated the platform library to be up-to-date following this run')
        return ['Success', additional_details['schedule_path_output'], success_statements]
    except Exception:
        return ['Error', 'Problem: %s' % traceback.format_exc(), success_statements]


if __name__ == "__main__":
    queuename = 'dryrun_queue_20221018_1'
    prep_return = start_stop_heater_shaker(queuename, '2')
    pprint(prep_return, width=200)
