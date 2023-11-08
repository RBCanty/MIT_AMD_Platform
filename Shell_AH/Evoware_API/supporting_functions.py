# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 15:23:44 2022

@author: Brent
"""

import traceback
import yaml
import pymongo
import datetime
import copy
import os

if __name__ == "__main__":
    from platform_library_functions import query_collection, query_document
    from useful_other_functions import solvent_rinse_location
    from tecan_functions import transfer_labware, wash_tips, liha_aspirate_dispense, special_lid_transfer_labware, \
        move_liha, special_transfer_labware, start_timer, wait_for_timer, stop_heater_shaker, large_waste_dump, \
        start_heater_shaker, set_heater_shaker_temperature, user_prompt, mca_move, tevac_string_return, \
        mix_liha, mca_droptips, mca_get_tips, mca_mix, mca_aspirate, mca_dispense, set_variable, \
        liha_detect_liquid, export_variable, liha_aspirate, liha_dispense, execute_vbs_script, roma_vector
else:
    try:
        from Evoware_API.platform_library_functions import query_collection, query_document
        from Evoware_API.useful_other_functions import solvent_rinse_location
        from Evoware_API.tecan_functions import transfer_labware, wash_tips, liha_aspirate_dispense, \
            special_lid_transfer_labware, move_liha, special_transfer_labware, start_timer, wait_for_timer, \
            stop_heater_shaker, large_waste_dump, start_heater_shaker, set_heater_shaker_temperature, user_prompt, \
            mca_move, tevac_string_return, mix_liha, mca_droptips, mca_get_tips, mca_mix, mca_aspirate, mca_dispense, \
            set_variable, liha_detect_liquid, export_variable, liha_aspirate, liha_dispense, execute_vbs_script, roma_vector
    except ModuleNotFoundError:
        from platform_library_functions import query_collection, query_document
        from useful_other_functions import solvent_rinse_location
        from tecan_functions import transfer_labware, wash_tips, liha_aspirate_dispense, special_lid_transfer_labware, \
            move_liha, special_transfer_labware, start_timer, wait_for_timer, stop_heater_shaker, large_waste_dump, \
            start_heater_shaker, set_heater_shaker_temperature, user_prompt, mca_move, tevac_string_return, \
            mix_liha, mca_droptips, mca_get_tips, mca_mix, mca_aspirate, mca_dispense, set_variable, \
            liha_detect_liquid, export_variable, liha_aspirate, liha_dispense, execute_vbs_script, roma_vector


def directory_check(path_to_check):
    if not os.path.exists(path_to_check):
        os.makedirs(path_to_check)
        return ['Success', 'Made a new directory: %s' % path_to_check]
    return ['Success', 'Directory already exists: %s' % path_to_check]


def get_queue_information(queue_name, operation_number, target_operation):
    success_statements = []
    # noinspection PyBroadException
    try:
        # Open the configuration file before starting anything
        config_filepath = r'C:\Users\admin\PycharmProjects\AMD_Control_Platform\venv\Shell_AH\Evoware_API\config.yaml'
        with open(config_filepath, 'r') as yaml_file:
            config = yaml.load(yaml_file, Loader=yaml.Loader)

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
            return ['Error', queue_return_statements, success_statements, '', '']
        else:
            success_statements.append('Found element %s in the queue collection' % queue_name)

        # Now we check to make sure that the plate is ready to be prepared
        plate_operations = queue_document['operations']
        queue_containers = queue_document['containers']

        # If it is the first operation in the queue then we just check it is not completed
        if operation_number == '1':
            if plate_operations[operation_number]['completed'] != 'no':
                success_statements.append('Operation %s is not ready to be executed, status %s' % (
                    operation_number, plate_operations[operation_number]['completed']))
                return ['Error', 'Queue completion error', success_statements, '', '']
        else:
            previous_operation = str(int(operation_number) - 1)
            previous_operation_status = plate_operations[previous_operation]['completed']
            if previous_operation_status == 'error':
                return ['Error', 'Operation %s has an error status, not ready for operation %s' %
                        (previous_operation, operation_number), success_statements, '', '']
            elif previous_operation_status == 'no':
                return ['Error', 'Previous operation %s not marked as complete in database, not ready for %s' %
                        (previous_operation, operation_number), success_statements, '', '']
            elif previous_operation_status == 'yes' and plate_operations[operation_number]['completed'] == 'no':
                success_statements.append('Operation %s is in valid position in queue document' % operation_number)
            else:
                return ['Error', 'The operation %s is already marked as complete' % operation_number,
                        success_statements, '', '']

        if plate_operations[operation_number]['operation'] != target_operation:
            return ['Error', 'Expected operation %s does not match queue document operation %s' %
                    (plate_operations[operation_number]['operation'], target_operation), success_statements, '', '']
        else:
            success_statements.append('Expected operation matches queue document operation %s' % target_operation)

        additional_prep_details = plate_operations[operation_number]['details']
        target_container = plate_operations[operation_number]['container']
        target_well_plate_type = queue_containers[target_container]['plate_type']
        success_statements.append('The %s operation is in a valid position and ready to be executed' % target_operation)

        # Here there is an optional statement from config to send to the testing folder
        current_time = datetime.datetime.now()
        if config['startup']['testing']:
            config_path_output = config['paths/files']['schedule_testing_output_path']
            schedule_path_output = config_path_output % current_time.strftime('%Y%m%d')
        else:
            config_path_output = config['paths/files']['schedule_output_path']
            folder_out_path = config_path_output % current_time.strftime('%Y%m%d')
            directory_check(folder_out_path)
            schedule_path_output = folder_out_path + r'\%s_%s_%s_%s.esc' % (
                queue_name, operation_number, target_operation, current_time.strftime('%Y%m%d%H%M%S'))

        # Then load all the documents from the platform library which includes empty
        # wellplates, reagents, solvents, and consumables. This could probably be made
        # more efficient later when it is needed but for now it will be fast enough
        full_library_mongo = {}
        important_collections = ['wellplates', 'solvents', 'reagents', 'consumables']
        for important_collection in important_collections:
            message = {'auth': {'port': '27017', 'user': username, 'password': password},
                       'request': {'request_type': 'query_collection', 'collection': important_collection}}
            document_return, extra_statements = query_collection(message['request']['collection'], platform_mongodb)
            if document_return == 'Error':
                return ['Error', extra_statements, success_statements, '', '']
            full_library_mongo[important_collection] = document_return
        platform_mongodb.close()

        return_information = {'carriers_labware_filepath': carriers_labware_filepath,
                              'empty_worktable_filepath': empty_worktable_filepath,
                              'target_container': target_container,
                              'target_well_plate_type': target_well_plate_type,
                              'queue_containers': queue_containers,
                              'schedule_path_output': schedule_path_output,
                              'wellplate_name': queue_containers[target_container]['container_name'],
                              'platform_mongodb_address': platform_mongodb_address,
                              'platform_mongodb_port': platform_mongodb_port,
                              'additional_prep_details': additional_prep_details,
                              'liquid_detection_path': config['paths/files']['liquid_detection_filepath']}
    except Exception:
        return ['Error', traceback.format_exc(), '', '', success_statements]
    return ['Success', full_library_mongo, queue_document, return_information, success_statements]


def convert_to_method_strings(simplified_schedule, updated_carriers, full_library_mongo, carriers_labware):
    try:
        method_strings = []
        bad_methods = []
        for item_index, item in enumerate(simplified_schedule):
            if item[0] == 'TRANSFER_LABWARE':
                return_statement, method_string = transfer_labware(updated_carriers,
                                                                   full_library_mongo[item[5]],
                                                                   item[0:5])
                if return_statement == 'Success':
                    # We need to quickly update the positions so that the other
                    # functions will know about worktable changes being made
                    updated_carriers[item[3][0]][item[3][2]]['labware_labels'][item[3][1] - 1] = copy.deepcopy(
                        updated_carriers[item[2][0]][item[2][2]]['labware_labels'][item[2][1] - 1])
                    updated_carriers[item[2][0]][item[2][2]]['labware_labels'][item[2][1] - 1] = ''
                    updated_carriers[item[3][0]][item[3][2]]['labware_types'][item[3][1] - 1] = copy.deepcopy(
                        updated_carriers[item[2][0]][item[2][2]]['labware_types'][item[2][1] - 1])
                    updated_carriers[item[2][0]][item[2][2]]['labware_types'][item[2][1] - 1] = ''
            elif item[0] == 'SHAKER_TEMPERATURE':
                return_statement, method_string = set_heater_shaker_temperature(item[1], item[2], item[3])
            elif item[0] == 'START_SHAKER':
                return_statement, method_string = start_heater_shaker(item[1], item[2])
            elif item[0] == 'STOP_SHAKER':
                return_statement, method_string = stop_heater_shaker(item[1])
            elif item[0] == 'SPECIAL_TRANSFER_LABWARE':
                return_statement, method_string = special_transfer_labware(updated_carriers, item[1], item)
            elif item[0] == 'START_TIMER':
                return_statement, method_string = start_timer(item[1])
            elif item[0] == 'WAIT_FOR_TIMER':
                return_statement, method_string = wait_for_timer(*item[1::])
            elif item[0] == 'EXECUTE_VBS_SCRIPT':
                return_statement, method_string = execute_vbs_script(*item[1::])
            elif item[0] == 'TIP_RINSE':
                return_statement, method_string = wash_tips()
            elif item[0] == 'MOVE_LiHa':
                return_statement, method_string = move_liha()
            elif item[0] == 'MIX_LiHa':
                return_statement, method_string = mix_liha(updated_carriers, carriers_labware,
                                                           item[1], item[2], item[3], item[4])
            elif item[0] == 'USER_PROMPT':
                return_statement, method_string = user_prompt(item[1], sound=3)
            elif item[0] == 'SET_VARIABLE':
                return_statement, method_string = set_variable(*item[1::])
            elif item[0] == 'EXPORT_VARIABLE':
                return_statement, method_string = export_variable(*item[1::])
            elif item[0] == 'MCA_ASPIRATE':
                return_statement, method_string = mca_aspirate(*item[1::])
            elif item[0] == 'MCA_DISPENSE':
                return_statement, method_string = mca_dispense(*item[1::])
            elif item[0] == 'MCA_DROP_TIPS':
                return_statement, method_string = mca_droptips(updated_carriers)
            elif item[0] == 'MCA_MOVE':
                if len(item) > 1:
                    destination = item[1]
                else:
                    destination = None
                return_statement, method_string = mca_move(location=destination)
            elif item[0] == 'ROMA_MOVE':
                retrun_statement, method_string = roma_vector(*item[1::])
            elif item[0] == 'MCA_GET_TIPS':
                return_statement, method_string = mca_get_tips(item[1])
            elif item[0] == 'MCA_MIX':
                return_statement, method_string = mca_mix(*item[1::])
            elif item[0] == 'LARGE_WASTE_DUMP':
                return_statement, method_string = large_waste_dump(*item[1::], updated_carriers)
            elif item[0] == 'TEVAC_FUNCTION':
                return_statement, method_string = tevac_string_return(*item[1::])
            elif item[0] == 'SPECIAL_LID_TRANSFER':
                return_statement, method_string = special_lid_transfer_labware(updated_carriers,
                                                                               full_library_mongo[item[7]],
                                                                               item)
            elif item[0] == 'LIHA_DETECT_LIQUID':
                return_statement, method_string = liha_detect_liquid(updated_carriers, carriers_labware, *item[1::])
            elif item[0] == 'MIX_LIHA':
                return_statement, method_string = mix_liha(updated_carriers, carriers_labware,
                                                           item[1], item[2], item[3], item[4], cycles=item[5])

            # There are a couple of special ones that return more than one method,
            # so we will add them to method strings and continue to the next item
            elif item[0] == 'LiHa_ASPIRATE':
                return_statement, aspirate_val = liha_aspirate(updated_carriers, carriers_labware, item[1:])
                if return_statement == 'Error':
                    bad_methods.append(aspirate_val)
                    continue
                method_strings.append(aspirate_val)
                continue
            elif item[0] == 'LiHa_DISPENSE':
                return_statement, dispense_val = liha_dispense(updated_carriers, carriers_labware, item[1:])
                if return_statement == 'Error':
                    bad_methods.append(dispense_val)
                    continue
                method_strings.append(dispense_val)
                continue
            elif item[0] == 'LiHa_ASPIRATE_DISPENSE':
                if len(item) == 7:
                    return_statement, aspirate_val, dispense_val = liha_aspirate_dispense(updated_carriers,
                                                                                          carriers_labware,
                                                                                          item[1:-1],
                                                                                          liquidclass=item[-1])
                else:
                    return_statement, aspirate_val, dispense_val = liha_aspirate_dispense(updated_carriers,
                                                                                          carriers_labware,
                                                                                          item[1:])
                if return_statement == 'Error':
                    bad_methods.append(aspirate_val)
                    continue
                method_strings.append(aspirate_val)
                method_strings.append(dispense_val)
                continue
            elif item[0] == 'SOLVENT_RINSE':
                method_string = ''
                return_statement, solvent_rinse_origin = solvent_rinse_location(full_library_mongo, item[1])
                if return_statement == 'Error':
                    bad_methods.append(solvent_rinse_origin)
                    continue
                return_statement, aspirate_val, dispense_val = liha_aspirate_dispense(updated_carriers,
                                                                                      carriers_labware,
                                                                                      [solvent_rinse_origin, [1, 2, 30],
                                                                                       [1, 2, 3, 4, 5, 6, 7, 8],
                                                                                       [int(item[2])] * 8,
                                                                                       ['A1', 'B1', 'C1', 'D1', 'E1',
                                                                                        'F1', 'G1', 'H1']])
                if return_statement == 'Error':
                    bad_methods.append(aspirate_val)
                    continue
                method_strings.append(aspirate_val)
                method_strings.append(dispense_val)
                continue

            elif item[0] == 'MCA_SOFT_MIX':
                method_strings.append('BeginLoop("%s","soft_mixing");' % item[5])
                plate_type = item[6]
                if plate_type == 'Paradox Thermal Plate':
                    liquid_class = 'Adjusted Water free dispense DiTi Paradox %s'
                elif plate_type == '96 Well DeepWell':
                    liquid_class = 'Adjusted Water free dispense DiTi DeepWell %s'
                else:
                    bad_methods.append('Input wellplate %s not acceptable for LLE' % plate_type)
                    continue
                if item[4] == 'top':
                    return_statement, method_string = mca_aspirate(item[1], item[3], liquid_class % 'Top')
                    if return_statement == 'Error':
                        bad_methods.append('Problem with MCA_ASPIRATE with %s' % item)
                        continue
                    method_strings.append(method_string)
                    return_statement, method_string = mca_dispense(item[2], item[3], liquid_class % 'Bottom')
                    if return_statement == 'Error':
                        bad_methods.append('Problem with MCA_DISPENSE with %s' % item)
                        continue
                    method_strings.append(method_string)
                elif item[4] == 'bottom':
                    return_statement, method_string = mca_aspirate(item[1], item[3], liquid_class % 'Bottom')
                    if return_statement == 'Error':
                        bad_methods.append('Problem with MCA_ASPIRATE with %s' % item)
                        continue
                    method_strings.append(method_string)
                    return_statement, method_string = mca_dispense(item[2], item[3], liquid_class % 'Top')
                    if return_statement == 'Error':
                        bad_methods.append('Problem with MCA_DISPENSE with %s' % item)
                        continue
                    method_strings.append(method_string)
                method_strings.append('EndLoop();')

            elif item[0] == 'MCA_LAYER':
                plate_type = item[5]
                if plate_type == 'Paradox Thermal Plate':
                    liquid_class = 'Adjusted Water free dispense DiTi Paradox %s'
                elif plate_type == '96 Well DeepWell':
                    liquid_class = 'Adjusted Water free dispense DiTi DeepWell %s'
                else:
                    bad_methods.append('Input wellplate %s not acceptable for LLE' % plate_type)
                    continue
                if item[4] == 'top':
                    return_statement, method_string = mca_aspirate(item[1], item[3], liquid_class % 'Top')
                    if return_statement == 'Error':
                        bad_methods.append('Problem with MCA_ASPIRATE with %s' % item)
                        continue
                    method_strings.append(method_string)
                    return_statement, method_string = mca_dispense(item[2], item[3])
                    if return_statement == 'Error':
                        bad_methods.append('Problem with MCA_DISPENSE with %s' % item)
                        continue
                    method_strings.append(method_string)
                elif item[4] == 'bottom':
                    return_statement, method_string = mca_aspirate(item[1], item[3], liquid_class % 'Bottom')
                    if return_statement == 'Error':
                        bad_methods.append('Problem with MCA_ASPIRATE with %s' % item)
                        continue
                    method_strings.append(method_string)
                    return_statement, method_string = mca_dispense(item[2], item[3])
                    if return_statement == 'Error':
                        bad_methods.append('Problem with MCA_DISPENSE with %s' % item)
                        continue
                    method_strings.append(method_string)
                method_strings.append('EndLoop();')

            else:
                bad_methods.append('Method %s not implemented' % item[0])
                return_statement = 'Error'
                continue

            # Now we can check the method_string return for errors
            if return_statement == 'Error':  # noqa
                bad_methods.append('item #%s : %s' % (item_index, method_string))  # noqa
                continue
            else:
                method_strings.append(method_string)

    except Exception:  # noqa
        return ['Error', traceback.format_exc(), '']

    if len(bad_methods) != 0:
        return ['Error', bad_methods, '']
    else:
        return ['Success', method_strings, updated_carriers]
