# -*- coding: utf-8 -*-
"""
Created on 09/27/2022

Code that converts the simple workflows in the queue operations

@author: Brent Koscher
"""

import traceback
from copy import deepcopy
from pprint import pprint


def operations_list_to_details(workflow_type, operations_list, workflow_details):
    try:
        queue_operations = {}
        for operation_index, operation in enumerate(operations_list):
            print(operation)
            operation_key = str(operation_index + 1)
            # -----------------------------------------------------------------------------------------------
            # --------------------------- Below are the liquid handler operations ---------------------------
            if operation[0] == 'prepare_wellplate':
                pprint(workflow_details)
                extra_prep_details = {}
                if 'additional_prep_details' in workflow_details.keys():
                    if 'washing_frequency' in workflow_details['additional_prep_details'].keys():
                        washing_frequency = workflow_details['additional_prep_details']['washing_frequency']
                        if washing_frequency == 'nonpolar_stepdown' and 'nonpolar_solvents' in workflow_details['additional_prep_details'].keys():
                            extra_prep_details.update({'nonpolar_solvents': workflow_details['additional_prep_details']['nonpolar_solvents']})
                    else:
                        washing_frequency = None
                    if 'single_transfer' in workflow_details.keys():
                        extra_prep_details['single_transfer'] = workflow_details['single_transfer']
                else:
                    washing_frequency = None
                prepare_details = {'empty_wells': True,
                                   'washing_frequency': washing_frequency,
                                   'tip_prep_rinse': True}
                prepare_details.update(extra_prep_details)
                if operation[2] == 'photochemistry':
                    # TODO What are the special preparations needed for photochemistry? None are obvious
                    pass
                operation_stub = {'operation': 'prepare_wellplate',
                                  'container': operation[1],
                                  'details': prepare_details,
                                  'completed': 'no',
                                  'time_est': 3600,
                                  'agent': 'Lh',
                                  'start_time': None,
                                  'end_time': None}
                if 'inert' in operation or 'inert_atmosphere' in operation:
                    operation_stub['details'].update({'inert_atmosphere': True,
                                                      'inert_purge_time': 60 * 10,
                                                      'solvent_waiting_time': 60 * 10})
                if 'liquid_class' in workflow_details:
                    operation_stub['details'].update({'liquid_class': workflow_details['liquid_class']})
                if 'transfer_type' in workflow_details.keys():
                    operation_stub['details'].update({workflow_details['transfer_type']: True})
                if 'single_transfer' in workflow_details.keys() and workflow_details['single_transfer']:
                    operation_stub['details'].update({'single_transfer': True})
            elif operation[0] == 'prepare_characterization_plate':
                operation_stub = {'agent': 'Lh',
                                  'operation': 'prepare_characterization_plate',
                                  'container': operation[1],
                                  'details': {'empty_wells': True,
                                              'tip_prep_rinse': True,
                                              'copy_next_n_operations': operation[2]},
                                  'completed': 'no',
                                  'time_est': 60 * 30,
                                  'start_time': None,
                                  'end_time': None}
            elif operation[0] == 'copy_wellplate':
                operation_stub = {'agent': 'Lh',
                                  'operation': 'copy_wellplate',
                                  'container': operation[1],
                                  # TODO Actually look at the inputs
                                  # Queues to update
                                  # Solvent to add
                                  # Total volume
                                  'details': {},
                                  'completed': 'no',
                                  'time_est': 60 * 30,
                                  'start_time': None,
                                  'end_time': None}
                if operation[3]:
                    operation_stub['details'] = deepcopy(operation[4]['aliquot_details'])
                    operation_stub['details'].update({'target_container': operation[2]})
                    if 'queues_to_update' in operation[4].keys():
                        operation_stub['details']['queues_to_update'] = operation[4]['queues_to_update']
                    iteration_number = operation[4]['aliquot_iteration']
                    print(iteration_number, len(operation[4]['aliquot_time_series']))
                    if iteration_number < len(operation[4]['aliquot_time_series']) - 1:
                        schedule_time = operation[4]['aliquot_time_series'][iteration_number + 1] - operation[4]['aliquot_time_series'][iteration_number]
                        #operation_stub['details'].update({'schedule_time': schedule_time, 'is_paired': 'yes'})

            elif operation[0] == 'transfer_wellplate':
                operation_stub = {'operation': 'transfer_wellplate',
                                  'container': operation[1],
                                  'details': {'target_destination': operation[2]},
                                  'completed': 'no',
                                  'time_est': 60,
                                  'agent': 'Lh',
                                  'start_time': None,
                                  'end_time': None}
                if 'is_paired' in operation:
                    operation_stub['details']['is_paired'] = 'yes'
                if operation[2] == 'heater_shaker' and workflow_type == 'photochemistry':
                    operation_stub['details']['temperature'] = 20
                elif operation[2] == 'heater_shaker':
                    operation_stub['details']['temperature'] = operation[4]
                if 'inert' in operation:
                    operation_stub['details']['disassemble_inert'] = True
                if operation[3] == 'inert_atmosphere':
                    operation_stub['details']['disassemble_inert'] = True
            elif operation[0] == 'start_stop_heater_shaker':
                if operation[2] == 'start':
                    if operation[3] == 'mix':
                        temperature_set = 20
                    else:
                        temperature_set = workflow_details['reaction_temperature']
                    operation_details = {'power': 'on',
                                         'rpms': 600,
                                         'temperature': temperature_set}
                    if 'aliquot_time_series' in workflow_details.keys():
                        if workflow_details['aliquot_time_series'][0] != 0:
                            schedule_time = workflow_details['aliquot_time_series'][0]
                            operation_details.update({'is_paired': 'yes',
                                                      'schedule_time': schedule_time})
                    else:
                        operation_details.update({'is_paired': 'yes',
                                                  'schedule_time': operation[4]})
                else:
                    operation_details = {'power': 'off',
                                         'rpms': 200,
                                         'temperature': 20}
                operation_stub = {'operation': 'start_stop_heater_shaker',
                                  'container': operation[1],
                                  'details': operation_details,
                                  'completed': 'no',
                                  'time_est': 60,
                                  'agent': 'Lh',
                                  'start_time': None,
                                  'end_time': None}
            elif operation[0] == 'filter_wellplate':
                operation_stub = {'operation': 'filter_wellplate',
                                  'container': operation[1],
                                  'details':
                                      {'target_container': operation[2],
                                       'initial_solvent_volume': 400,
                                       'initial_filter_solvent': 'chloroform',
                                       'mca_filter_volume': 250,
                                       'mca_mix': True,
                                       'tevac_filter_time': 180,
                                       'tevac_location': 'back'
                                       },
                                  'completed': 'no',
                                  'time_est': 1800,
                                  'agent': 'Lh',
                                  'start_time': None,
                                  'end_time': None}
                if operation[3] == 'hplc_plate':
                    operation_stub['details']['hplc_container'] = operation[3]
                    operation_stub['details']['hplc_container_solvent'] = operation[4]
                if workflow_type == 'photochemistry':
                    operation_stub['details']['initial_solvent_volume'] = 250
                    operation_stub['details']['initial_filter_solvent'] = 'dmso'


            # -----------------------------------------------------------------------------------------------
            # ---------------------------- Below are the robotic arm operations -----------------------------
            elif operation[0] == 'move_wellplate':
                operation_stub = {'operation': 'move_wellplate',
                                  'container': operation[1],
                                  'details': {'target_destination': operation[2],
                                              'is_paired': 'yes'},
                                  'completed': 'no',
                                  'time_est': 60,
                                  'agent': 'Ra',
                                  'start_time': None,
                                  'end_time': None}
                if 'is_paired' in operation:
                    operation_stub['details']['is_paired'] = 'yes'

            # -----------------------------------------------------------------------------------------------
            # ---------------------------- Below are the HPLC Module operations -----------------------------
            elif operation[0] == 'run_analytical_batch':
                operation_stub = {'operation': 'run_analytical_batch',
                                  'container': operation[1],
                                  'details': {'update_plates': operation[2],
                                              'queue_number': operation[3]},
                                  'completed': 'no',
                                  'time_est': 96 * 60 * 2,
                                  'agent': 'LC',
                                  'start_time': None,
                                  'end_time': None}
                if 'hplc_processing' in workflow_details.keys() and workflow_details['hplc_processing'] != '':
                    operation_stub['details']['hplc_processing_method'] = workflow_details['hplc_processing']
            elif operation[0] == 'run_uplc_batch':
                operation_stub = {'operation': 'run_uplc_batch',
                                  'container': operation[1],
                                  'details': {'update_plates': operation[2],
                                              'queue_number': operation[3]},
                                  'completed': 'no',
                                  'time_est': 96 * 60 * 2,
                                  'agent': 'LC',
                                  'start_time': None,
                                  'end_time': None}
                if 'hplc_processing' in workflow_details.keys() and workflow_details['hplc_processing'] != '':
                    operation_stub['details']['hplc_processing_method'] = workflow_details['hplc_processing']

            # -----------------------------------------------------------------------------------------------
            # ---------------------------- Below are the photo-reactor operations ----------------------------
            elif operation[0] == 'Ph.prepare_to_send_or_receive':
                operation_stub = {'operation': 'Ph.prepare_to_send_or_receive',
                                  'container': operation[1],
                                  'details': {},
                                  'completed': 'no',
                                  'time_est': 60,
                                  'agent': 'Ph',
                                  'start_time': None,
                                  'end_time': None}
                if 'is_paired' in operation:
                    operation_stub['details']['is_paired'] = 'yes'
            elif operation[0] == 'Ph.run':
                photo_settings = {'incubation_time': 0,
                                  'reaction_time': 0,
                                  'stabilization_time': 0,
                                  'reaction_temperature': 20,
                                  'led_power': 0}
                reaction_time = 60
                for detail_key in workflow_details.keys():
                    if detail_key in photo_settings.keys():
                        photo_settings[detail_key] = workflow_details[detail_key]
                        if '_time' in detail_key:
                            reaction_time += workflow_details[detail_key]
                if 'led_power_setting' in workflow_details.keys():
                    photo_settings['led_power'] = workflow_details['led_power_setting']
                operation_stub = {'operation': 'Ph.run',
                                  'container': operation[1],
                                  'details': {'photo_settings': photo_settings},
                                  'completed': 'no',
                                  'time_est': reaction_time,
                                  'agent': 'Ph',
                                  'start_time': None,
                                  'end_time': None}
            elif operation[0] == 'Ph.go_to_initial_state':
                operation_stub = {'operation': 'Ph.go_to_initial_state',
                                  'container': operation[1],
                                  'details': {},
                                  'completed': 'no',
                                  'time_est': 60,
                                  'agent': 'Ph',
                                  'start_time': None,
                                  'end_time': None}
                if 'is_paired' in operation:
                    operation_stub['details']['is_paired'] = 'yes'

            # -----------------------------------------------------------------------------------------------
            # --------------------------- Below are the thermal reactor operations --------------------------
            elif operation[0] == 'Th.prepare_to_send_receive':
                operation_stub = {'operation': 'Th.prepare_to_send_receive',
                                  'container': operation[1],
                                  'details': {},
                                  'completed': 'no',
                                  'time_est': 60,
                                  'agent': 'Th',
                                  'start_time': None,
                                  'end_time': None}
                if 'is_paired' in operation:
                    operation_stub['details']['is_paired'] = 'yes'
            elif operation[0] == 'Th.run':
                operation_stub = {'operation': 'Th.run',
                                  'container': operation[1],
                                  'details': {'heating_profile': {
                                      'venting_time': 300,
                                      'heat_temperature': operation[3],
                                      'heat_hold_time': operation[2],
                                      'evap_temperature': 160,  # Boiling point of solvent
                                      'evap_hold_time': -1,  # 60 * 60 * 2, # Some time estimation from solvent
                                      'safe_temperature': 40
                                  }
                                  },
                                  'completed': 'no',
                                  'time_est': 60 * 60 * 12,
                                  'agent': 'Th',
                                  'start_time': None,
                                  'end_time': None}
            elif operation[0] == 'Th.go_to_initial_state':
                operation_stub = {'operation': 'Th.go_to_initial_state',
                                  'container': operation[1],
                                  'details': {},
                                  'completed': 'no',
                                  'time_est': 60,
                                  'agent': 'Th',
                                  'start_time': None,
                                  'end_time': None}
                if 'is_paired' in operation:
                    operation_stub['details']['is_paired'] = 'yes'

            # -----------------------------------------------------------------------------------------------
            # -------------------------- Below are the storage carousel operations --------------------------
            elif operation[0] == 'Ss.stow':
                operation_stub = {'agent': 'Ss',
                                  'operation': 'Ss.stow',
                                  'container': operation[1],
                                  'details': {},
                                  'completed': 'no',
                                  'time_est': 60,
                                  'start_time': None,
                                  'end_time': None}

            # -----------------------------------------------------------------------------------------------
            # ------------------------------ Below are the MCN level operations -----------------------------
            elif operation[0] == 'complete_queue':
                operation_stub = {'operation': 'complete_queue',
                                  'container': '',
                                  'details': {},
                                  'completed': 'no',
                                  'time_est': 30,
                                  'agent': 'MC',
                                  'start_time': None,
                                  'end_time': None}

            elif operation[0] == 'soft_wait':
                operation_stub = {'operation': 'soft_wait',
                                  'container': '',
                                  'details': {},
                                  'completed': 'no',
                                  'time_est': 30,
                                  'agent': 'MC',
                                  'start_time': None,
                                  'end_time': None}
                if 'wait_time' in operation[1].keys() and operation[1]['wait_time'] > 30:
                    operation_stub['time_est'] = operation[1]['wait_time']

            # -----------------------------------------------------------------------------------------------
            # ----------------------------- Below are the plate reader operations ---------------------------
            elif operation[0] == 'Pr.prepare_to_send_or_receive':
                operation_stub = {'agent': 'Pr',
                                  'operation': 'Pr.prepare_to_send_or_receive',
                                  'container': operation[1],
                                  'details': {'direction': operation[2]},
                                  'completed': 'no',
                                  'time_est': 30,
                                  'start_time': None,
                                  'end_time': None}
            elif operation[0] == 'Pr.go_to_initial_state':
                operation_stub = {'agent': 'Pr',
                                  'operation': 'Pr.go_to_initial_state',
                                  'container': operation[1],
                                  'details': {},
                                  'completed': 'no',
                                  'time_est': 30,
                                  'start_time': None,
                                  'end_time': None}
                if 'is_paired' in operation:
                    operation_stub['details']['is_paired'] = 'yes'
            elif operation[0] == 'Pr.check_wellplate':
                operation_stub = {'agent': 'Pr',
                                  'operation': 'Pr.check_wellplate',
                                  'container': operation[1],
                                  'details': operation[2],
                                  'completed': 'no',
                                  'time_est': 5 * 60 * 1.05,
                                  'start_time': None,
                                  'end_time': None}
            elif operation[0] == 'Pr.run_platereader':
                operation_stub = {'agent': 'Pr',
                                  'operation': 'Pr.run_platereader',
                                  'container': operation[1],
                                  'details': {},
                                  'completed': 'no',
                                  'time_est': 60 * 30,
                                  'start_time': None,
                                  'end_time': None}
                if operation[2] == 'abs':
                    operation_stub['details'] = {'assay': 'peaks', 'mode': operation[2]}
                    operation_stub['time_est'] = 10 * 60
                elif operation[2] == 'photodeg_abs':
                    operation_stub['details'] = operation[3]
                    operation_stub['time_est'] = (operation_stub['details']['iterations'] * 60 * (
                            5 + operation_stub['details']['interval'])) * 1.05
                else:
                    return ['Error', 'Characterization of type %s not implemented' % operation[2], '']
            elif operation[0] == 'Pr.analyze_data':
                operation_stub = {'agent': 'Pr',
                                  'operation': 'Pr.analyze_data',
                                  'container': operation[1],
                                  'details': {},
                                  'completed': 'no',
                                  'time_est': 5 * 60 * 1.05,
                                  'start_time': None,
                                  'end_time': None}
                if operation[2] == 'kinetic_abs':
                    operation_stub['details'] = operation[3]
                elif operation[2] == 'static_abs':
                    operation_stub['details'] = operation[3]
            # Otherwise return an error
            else:
                print(operation)
                return ['Error', 'Key "%s" not present in operation table' % operation[0]]
            queue_operations[operation_key] = operation_stub
        return ['Success', queue_operations]
    except:
        return ['Error', traceback.format_exc()]
