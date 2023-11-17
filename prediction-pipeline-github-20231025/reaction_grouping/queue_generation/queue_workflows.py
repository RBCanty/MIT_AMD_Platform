# -*- coding: utf-8 -*-
"""
Created on 09/27/2022

Code that defines the different allowed workflows for queue generation

@author: Brent Koscher
"""

import sys
import traceback
from pprint import pprint
import re
from copy import deepcopy
import datetime
from platform_operations import operations_list_to_details


def queue_workflows(workflow_type, grouping_document, grouping_information):
    return_statement, reaction_plate_queues, reagents_needed = reaction_queues(workflow_type, grouping_document,
                                                                               grouping_information)
    if return_statement != 'Success':
        return ['Error', reaction_plate_queues, '']
    return ['Success', reaction_plate_queues, reagents_needed]


def find_queue_dependencies(queue_list, workflow_details):
    # Now we need to check the queue_documents for dependencies
    # TODO Change to check for the reaction dependencies between wellplates in queue
    dependencies = {}
    queue_documents_list = list(queue_list.keys())
    if 'aliquot_time_series' in workflow_details.keys():
        for queue_index, queue_document in enumerate(queue_documents_list):
            dependencies[queue_document] = []
    else:
        for queue_index, queue_document in enumerate(queue_documents_list):
            if queue_index == 0:
                dependencies[queue_document] = []
            else:
                dependencies[queue_document] = queue_documents_list[0:queue_index]
    return ['Success', dependencies]


def intermediate_purification(grouping_information, workflow_details, wellplate_index, number_of_groups):
    extra_containers_needed = []
    intermediate_purification_operations = []
    if 'hplc_analysis' in workflow_details and workflow_details['hplc_analysis']:
        if 'intermediate_purification' in workflow_details and workflow_details['intermediate_purification']:
            target_container = 'reaction_plate'
            extra_containers_needed.extend(['filtrate_plate', 'solid_filter_plate'])
            if 'separate_hplc_plate' in workflow_details and workflow_details['separate_hplc_plate']:
                intermediate_purification_operations.append(['filter_wellplate', target_container,
                                                             'filtrate_plate', 'hplc_plate', 'dmso'])
                target_container = 'hplc_plate'
                extra_containers_needed.append('hplc_plate')
            else:
                intermediate_purification_operations.append(['filter_wellplate', target_container,
                                                             'filtrate_plate', '', ''])
                target_container = 'filtrate_plate'
        else:
            target_container = 'reaction_plate'
        intermediate_purification_operations.append(['transfer_wellplate', target_container,
                                                     'transfer_hotel', '', 'is_paired'])
        intermediate_purification_operations.append(['move_wellplate', target_container, 'autosampler', 'is_paired'])
        if target_container == 'hplc_plate':
            containers_to_update = ['filtrate_plate', 'hplc_plate']
        else:
            containers_to_update = [target_container]
        if 'untargeted_hplc_analysis' in workflow_details.keys():
            extra_analysis = workflow_details['untargeted_hplc_analysis']
        else:
            extra_analysis = ''
        if 'hplc_processing' in workflow_details.keys() and workflow_details['hplc_processing'] == 'photochemistry':
            intermediate_purification_operations.append(['run_uplc_batch', target_container, containers_to_update,
                                                         [wellplate_index + 1, number_of_groups], extra_analysis])
        else:
            intermediate_purification_operations.append(['run_analytical_batch', target_container, containers_to_update,
                                                         [wellplate_index + 1, number_of_groups], extra_analysis])
        intermediate_purification_operations.append(['move_wellplate', target_container, 'lpx', 'is_paired'])
        intermediate_purification_operations.append(['Ss.stow', target_container])
    return ['Success', intermediate_purification_operations, extra_containers_needed]


def group_to_wellplate(wellplate_mapping, reaction_grouping, reaction_queue_name):
    try:
        wellplate_contents = {}
        reagents_needed = {}
        for well_key in wellplate_mapping:
            if wellplate_mapping[well_key] == '':
                continue
            mapped_reaction_key = wellplate_mapping[well_key]
            well_details = reaction_grouping[mapped_reaction_key]
            wellplate_contents[well_key] = {'confirmation': 'none',
                                            'total_volume': -1,
                                            'target_product': [],
                                            'reagents': [],
                                            'solvents': [],
                                            'reaction_smiles': '',
                                            'templates': []}
            fields_to_check = ['reactants', 'reagents', 'catalysts']

            well_reagents = []
            for field_to_check in fields_to_check:
                if field_to_check in well_details.keys() and well_details[field_to_check]:
                    well_reagents.extend(well_details[field_to_check])
                    wellplate_contents[well_key]['predicted_%s' % field_to_check] = [chemical[0] for chemical in
                                                                                     well_details[field_to_check]]
            if well_reagents:
                wellplate_contents[well_key]['reagents'] = well_reagents

            if reaction_queue_name not in reagents_needed.keys():
                reagents_needed[reaction_queue_name] = {'reagents': {}, 'solvents': {}}
            for well_reagent in well_reagents:
                if well_reagent[0] not in reagents_needed[reaction_queue_name]['reagents']:
                    reagents_needed[reaction_queue_name]['reagents'][well_reagent[0]] = 0
                reagents_needed[reaction_queue_name]['reagents'][well_reagent[0]] += well_reagent[1]

            if 'solvents' in well_details.keys() and well_details['solvents']:
                for solvent in well_details['solvents']:
                    if type(solvent) == str:
                        wellplate_contents[well_key]['solvents'].append([solvent, 1])
                    elif type(solvent) == list and len(solvent) == 1:
                        wellplate_contents[well_key]['solvents'].append([solvent[0], 1])
                    else:
                        wellplate_contents[well_key]['solvents'].append(solvent)
                wellplate_contents[well_key]['predicted_solvents'] = [chemical[0] for chemical in
                                                                      wellplate_contents[well_key]['solvents']]
            else:
                wellplate_contents[well_key]['predicted_solvents'] = []
            wellplate_contents[well_key]['total_volume'] = well_details['total_volume']

            if 'reaction_smiles' in well_details.keys():
                wellplate_contents[well_key]['reaction_smiles'] = well_details['reaction_smiles']
                target_product = well_details['reaction_smiles'].split('>>')[1]
                reactants_scale = min([chemical[1] for chemical in well_details['reactants']])
                wellplate_contents[well_key]['target_product'] = [[target_product, reactants_scale]]
                wellplate_contents[well_key]['product_remaining'] = reactants_scale
            else:
                wellplate_contents[well_key]['reaction_smiles'] = '.'.join([el[0] for el in well_details['reactants']])
                wellplate_contents[well_key]['target_product'] = None
                reactants_scale = min([chemical[1] for chemical in well_details['reactants']])
                wellplate_contents[well_key]['product_remaining'] = reactants_scale

            if 'templates' in well_details.keys():
                wellplate_contents[well_key]['templates'] = well_details['templates']

            if 'final_product' in well_details.keys():
                wellplate_contents[well_key]['final_product'] = well_details['final_product']
            else:
                wellplate_contents[well_key]['final_product'] = [0, 'no', '', 1]
    except:
        return ['Error', 'Issue with wellplate grouping: %s' % traceback.format_exc(), '']
    return ['Success', wellplate_contents, reagents_needed]


def get_wellplate_mapping(workflow_details, reaction_grouping):
    # Check the target wellplate type and see if it is an acceptable type of labware
    # TODO Replace some of this section with a lookup to the Evo ware Library file
    if 'wellplate_type' not in workflow_details.keys():
        return ['Error', 'Target wellplate type was not specified', '']
    allowed_wellplate_types = ['96 Well Microplate', 'Paradox Thermal Plate', '96 Well DeepWell', '96 Well MediumWell']
    if workflow_details['wellplate_type'] not in allowed_wellplate_types:
        return ['Error', 'Wellplate of type %s is not allowed %s' % (workflow_details['wellplate_type'],
                                                                     allowed_wellplate_types)]
    normal_wellplate_types = ['96 Well Microplate', 'Paradox Thermal Plate', '96 Well DeepWell', '96 Well MediumWell']

    # Figure out the allowed layout of the found wellplate type
    if workflow_details['wellplate_type'] in normal_wellplate_types:
        n_rows, n_cols, n_wells = 8, 12, 96
    else:
        return ['Error', 'Wellplate of type %s has not been given details' % workflow_details['wellplate_type']]
    wellplate_mapping = {}
    col, row = 1, 1
    for i in range(0, n_wells):
        wellplate_mapping[chr(ord('A') + row - 1) + str(col)] = ''
        if row == n_rows:
            row = 1
            col += 1
        else:
            row += 1

    # Check to see if there are specified wells already in the reaction grouping, giving priority
    # to pre-defined well keys and then filling with remaining reactions
    for reaction_key in reaction_grouping.keys():
        if reaction_key in wellplate_mapping.keys():
            wellplate_mapping[reaction_key] = reaction_key
    for reaction_key in reaction_grouping.keys():
        if reaction_key in wellplate_mapping.keys():
            continue
        for wellplate_key in wellplate_mapping.keys():
            if not wellplate_mapping[wellplate_key]:
                wellplate_mapping[wellplate_key] = reaction_key
                break
        else:
            return ['Error', 'Something is wrong with wellplate mapping with %s' % reaction_key]
    return ['Success', wellplate_mapping]


def reaction_queues(workflow_type, grouping_document, grouping_information):
    # First we will define the basic template for some predefined containers
    predefined_containers = {'reaction_plate': {'plate_type': '',
                                                'container_name': None,
                                                'contents': {}},
                             'filtrate_plate': {'plate_type': '96 Well Microplate',
                                                'container_name': None,
                                                'contents': {}},
                             'solid_filter_plate': {'plate_type': '96 Well Filtration Plate',
                                                    'container_name': None,
                                                    'contents': {}},
                             'extraction_plate': {'plate_type': '96 Well Microplate',
                                                  'container_name': None,
                                                  'contents': {}},
                             'characterization_plate': {'plate_type': '96 Well Microplate',
                                                        'container_name': None,
                                                        'contents': {}},
                             'aliquot_plate': {'plate_type': '96 Well Microplate',
                                               'container_name': None,
                                               'contents': {}},
                             'hplc_plate': {'plate_type': '96 Well PCR Plate',
                                            'container_name': None,
                                            'contents': {}}}

    # Now we need to give all the reactions a well to go into and using the workflow scripts we
    # need to generate the operations to run those reactions
    reagents_needed = {}
    queues_return = {}
    number_of_groups = len(list(grouping_document['wellplate_grouping'].keys()))
    current_time = datetime.datetime.now()
    for wellplate_index, wellplate_group in enumerate(grouping_document['wellplate_grouping'].keys()):
        group_information = grouping_document['wellplate_grouping'][wellplate_group]
        reaction_grouping = deepcopy(group_information['reaction_grouping'])
        workflow_details = group_information['workflow_details']
        reaction_queue_name = '%s_%s_%s' % (grouping_document['campaign_name'],
                                            current_time.strftime('%Y%m%d'),
                                            wellplate_index + 1)
        if reaction_queue_name in queues_return.keys():
            return ['Error', 'Reaction queue naming has an unexpected error', '']
        queues_return[reaction_queue_name] = {'operations': {},
                                                       'containers': {},
                                                       'queue_name': reaction_queue_name,
                                                       'status': 'idle',
                                                       'date_created': current_time.strftime('%m/%d/%Y')}
        if workflow_details['workflow_type'] == 'reaction':
            # Now add the reaction plate type and predefined stub to the queue document
            queue_document_containers = queues_return[reaction_queue_name]['containers']
            queue_document_containers['reaction_plate'] = deepcopy(predefined_containers['reaction_plate'])
            queue_document_containers['reaction_plate']['plate_type'] = workflow_details['wellplate_type']

            # Now we need to convert the group of reactions into a wellplate mapping
            return_statement, wellplate_mapping = get_wellplate_mapping(workflow_details, reaction_grouping)
            if return_statement != 'Success':
                return ['Error', wellplate_mapping, '']

            # With the well mapping in hand, we are ready to build the queue document
            return_statement, plate_grouping, plate_reagents = group_to_wellplate(wellplate_mapping,
                                                                                  reaction_grouping,
                                                                                  reaction_queue_name)
            if return_statement != 'Success':
                return ['Error', plate_grouping, '']
            queue_document_containers['reaction_plate']['contents'] = plate_grouping
            reagents_needed.update(plate_reagents)

        elif workflow_details['workflow_type'] == 'characterization':
            # Now add the reaction plate type and predefined stub to the queue document
            queue_document_containers = queues_return[reaction_queue_name]['containers']
            queue_document_containers['characterization_plate'] = deepcopy(predefined_containers['characterization_plate'])
            queue_document_containers['characterization_plate']['plate_type'] = workflow_details['wellplate_type']
            queue_document_containers['characterization_plate']['targets_to_find'] = reaction_grouping

        # Take the incoming workflow type and build the operations list
        # TODO Include the characterization steps into the operations as well
        operations_needed = []
        additional_plate_operations = []
        additional_queue_name = ''

        if workflow_type == 'photochemistry':
            operations_needed.append(['prepare_wellplate', 'reaction_plate', 'photochemistry'])
            operations_needed.append(['transfer_wellplate', 'reaction_plate', 'heater_shaker', ''])
            operations_needed.append(['start_stop_heater_shaker', 'reaction_plate', 'start', 'mix', 300])
            operations_needed.append(['start_stop_heater_shaker', 'reaction_plate', 'stop', '', ''])
            operations_needed.append(['transfer_wellplate', 'reaction_plate', 'transfer_hotel', '', 'is_paired'])
            operations_needed.append(['Ph.prepare_to_send_or_receive', 'reaction_plate', 'is_paired'])
            operations_needed.append(['move_wellplate', 'reaction_plate', 'photo_reactor', 'is_paired'])
            operations_needed.append(['Ph.run', 'reaction_plate'])
            operations_needed.append(['Ph.prepare_to_send_or_receive', 'reaction_plate', 'is_paired'])
            operations_needed.append(['move_wellplate', 'reaction_plate', 'liquid_handler', 'is_paired'])
            operations_needed.append(['Ph.go_to_initial_state', 'reaction_plate', 'is_paired'])
            operations_needed.append(['transfer_wellplate', 'reaction_plate', 'storage_hotel', ''])
        elif workflow_type == 'dual_catalyst_exploration':
            operations_needed.append(['prepare_wellplate', 'reaction_plate', 'dual_catalyst_exploration', 'inert'])
            operations_needed.append(['start_stop_heater_shaker', 'reaction_plate', 'start', '', ''])
            if 'aliquot_time_series' in workflow_details.keys() and workflow_details['aliquot_time_series']:
                additional_queue_name = '%s_%s_%s-%s' % (grouping_document['campaign_name'],
                                                         current_time.strftime('%Y%m%d'),
                                                         wellplate_index + 1,
                                                         'aliquot')
                if additional_queue_name not in queues_return.keys():
                    queues_return[additional_queue_name] = {'operations': {},
                                                                     'containers': {},
                                                                     'queue_name': additional_queue_name,
                                                                     'status': 'error',
                                                                     'date_created': current_time.strftime('%m/%d/%Y')}
                for iteration in range(0, len(workflow_details['aliquot_time_series'])):
                    aliquot_plate_name = 'aliquot_plate_%s' % (iteration+1)
                    operations_needed.append(['copy_wellplate', 'reaction_plate', aliquot_plate_name,
                                              'dual_catalyst_exploration', {'aliquot_iteration': iteration,
                                                                            'aliquot_time_series': workflow_details['aliquot_time_series'],
                                                                            'analysis_queue_name': additional_queue_name,
                                                                            'reaction_queue_name': reaction_queue_name,
                                                                            'queues_to_update': [additional_queue_name, reaction_queue_name],
                                                                            'aliquot_details': workflow_details['aliquot_details']}])
                    
                    # temporary fix while is-paired is being weird
                    operations_needed.append(['soft_wait', {'wait_time': workflow_details['aliquot_time_series'][iteration]}])
                    #temporary addition end
                    queues_return[additional_queue_name]['containers'][aliquot_plate_name] = predefined_containers['aliquot_plate']
                    queues_return[reaction_queue_name]['containers'][aliquot_plate_name] = predefined_containers['aliquot_plate']
                    if 'optional_characterization' in workflow_details.keys():
                        char_type = workflow_details['optional_characterization']['characterization_type']
                        if char_type == 'uvvis':
                            additional_plate_operations.append(['Pr.prepare_to_send_or_receive', aliquot_plate_name, 'receive'])
                            additional_plate_operations.append(['transfer_wellplate', aliquot_plate_name, 'spark', '', 'is_paired'])
                            additional_plate_operations.append(['Pr.run_platereader', aliquot_plate_name, 'abs'])
                            # additional_plate_operations.append(['Pr.analyze_data', aliquot_plate_name, 'static_abs', {'assay': 'peaks'}])
                            additional_plate_operations.append(['Pr.prepare_to_send_or_receive', aliquot_plate_name, 'send'])
                            additional_plate_operations.append(['transfer_wellplate', aliquot_plate_name, 'transfer_hotel', '', 'is_paired'])
                            additional_plate_operations.append(['Pr.go_to_initial_state', 'characterization_plate', 'is_paired'])
                    else:
                        additional_plate_operations.append(['transfer_wellplate', aliquot_plate_name, 'transfer_hotel', '', 'is_paired'])
                    additional_plate_operations.append(['move_wellplate', aliquot_plate_name, 'autosampler', 'is_paired'])
                    additional_plate_operations.append(['run_uplc_batch', aliquot_plate_name, [aliquot_plate_name],
                                                        [1, 1], workflow_details['hplc_analysis']])
                    
                    # temporary fix while lpx is being weird
                    additional_plate_operations.append(['move_wellplate', aliquot_plate_name, 'liquid_handler', 'is_paired'])
                    additional_plate_operations.append(['transfer_wellplate', aliquot_plate_name, 'storage_hotel', '', ''])
                    #additional_plate_operations.append(['move_wellplate', aliquot_plate_name, 'lpx', 'is_paired'])
                    #additional_plate_operations.append(['Ss.stow', aliquot_plate_name])

            operations_needed.append(['start_stop_heater_shaker', 'reaction_plate', 'stop', '', ''])
            #semi-permanent fix until we figure out how to take the small inert box off
            #operations_needed.append(['transfer_wellplate', 'reaction_plate', 'storage_hotel', 'inert', ''])
        elif workflow_type == 'standard' and workflow_details['workflow_type'] == 'reaction':
            if 'extra_conditions' in workflow_details.keys() and workflow_details['extra_conditions'] == 'inert_atmosphere':
                operations_needed.append(['prepare_wellplate', 'reaction_plate', workflow_details['extra_conditions']])
                operations_needed.append(['start_stop_heater_shaker', 'reaction_plate', 'start', '', workflow_details['reaction_time']])
                operations_needed.append(['start_stop_heater_shaker', 'reaction_plate', 'stop', '', ''])
                operations_needed.append(['transfer_wellplate', 'reaction_plate', 'storage_hotel', '', 'is_paired'])
            elif 'extra_conditions' in workflow_details.keys() and workflow_details['extra_conditions'] == 'low_temperature':
                operations_needed.append(['prepare_wellplate', 'reaction_plate', workflow_details['extra_conditions']])
                operations_needed.append(['start_stop_heater_shaker', 'reaction_plate', 'start', '', workflow_details['reaction_time']])
                operations_needed.append(['start_stop_heater_shaker', 'reaction_plate', 'stop', '', ''])
                operations_needed.append(['transfer_wellplate', 'reaction_plate', 'storage_hotel', '', 'is_paired'])
            elif workflow_details['reaction_temperature'] > 120:
                operations_needed.append(['prepare_wellplate', 'reaction_plate', ''])
                operations_needed.append(['transfer_wellplate', 'reaction_plate', 'transfer_hotel', ''])
                operations_needed.append(['Th.prepare_to_send_receive', 'reaction_plate', ''])
                operations_needed.append(['move_wellplate', 'reaction_plate', 'thermal_reactor', ''])
                operations_needed.append(['Th.run', 'reaction_plate', workflow_details['reaction_time'], workflow_details['reaction_temperature']])
                operations_needed.append(['Th.prepare_to_send_receive', 'reaction_plate', ''])
                operations_needed.append(['move_wellplate', 'reaction_plate', 'liquid_handler', ''])
                operations_needed.append(['transfer_wellplate', 'reaction_plate', 'storage_hotel', ''])
            else:
                operations_needed.append(['prepare_wellplate', 'reaction_plate', ''])
                operations_needed.append(['transfer_wellplate', 'reaction_plate', 'heater_shaker', '', workflow_details['reaction_temperature']])
                operations_needed.append(['start_stop_heater_shaker', 'reaction_plate', 'start', '', workflow_details['reaction_time']])
                operations_needed.append(['start_stop_heater_shaker', 'reaction_plate', 'stop', '', ''])
                operations_needed.append(['transfer_wellplate', 'reaction_plate', 'storage_hotel', '', 'is_paired'])

        elif workflow_type == 'standard' and workflow_details['workflow_type'] == 'characterization':
            characterizations_needed = workflow_details['characterizations']
            if 'uv_vis' in characterizations_needed.keys() and 'photodeg' not in characterizations_needed.keys():
                operations_needed.append(['prepare_characterization_plate', 'characterization_plate'])
                operations_needed.append(['Pr.prepare_to_send_or_receive', 'characterization_plate', 'receive'])
                operations_needed.append(['transfer_wellplate', 'characterization_plate', 'spark'])
                operations_needed.append(['Pr.check_wellplate', 'characterization_plate', {'assay': 'check', 'iter': 3, 'mode': 'abs'}])
                operations_needed.append(['Pr.run_platereader', 'characterization_plate', 'abs'])
                operations_needed.append(['Pr.analyze_data', 'characterization_plate', 'static_abs', {'assay': 'peaks'}])
            if 'photodeg' in characterizations_needed.keys():
                operations_needed.append(['prepare_characterization_plate', 'characterization_plate'])
                operations_needed.append(['Pr.prepare_to_send_or_receive', 'characterization_plate', 'receive', 'is_paired'])
                operations_needed.append(['transfer_wellplate', 'characterization_plate', 'spark', 'is_paired'])
                operations_needed.append(['Pr.check_wellplate', 'characterization_plate', {'assay': 'check', 'iter': 3, 'mode': 'abs'}])
                operations_needed.append(['Pr.run_platereader', 'characterization_plate', 'photodeg_abs', {'iterations': 15,
                                                                                                           'interval': 15,
                                                                                                           'is_photo': True,
                                                                                                           'temperature': 40,
                                                                                                           'photo_temp': 40,
                                                                                                           'wait_for_t': True,
                                                                                                           'assay': 'photodeg',
                                                                                                           'mode': 'abs'}])
                operations_needed.append(['Pr.analyze_data', 'characterization_plate', 'kinetic_abs', {'assay': 'photodeg', 'mode': 'abs'}])
                if 'uv_vis' in characterizations_needed.keys():
                    operations_needed.append(['Pr.analyze_data', 'characterization_plate', 'static_abs', {'from_assay': 'photodeg',
                                                                                                          'assay': 'peaks',
                                                                                                          'mode': 'abs'}])
            if operations_needed:
                operations_needed.append(['Pr.prepare_to_send_or_receive', 'characterization_plate', 'send', 'is_paired'])
                operations_needed.append(['transfer_wellplate', 'characterization_plate', 'transfer_hotel', 'is_paired'])
                operations_needed.append(['Pr.go_to_initial_state', 'characterization_plate', 'is_paired'])
                operations_needed.append(['move_wellplate', 'characterization_plate', 'lpx', 'is_paired'])
                operations_needed.append(['Ss.stow', 'characterization_plate'])
                operations_needed[0].append(len(operations_needed))
        else:
            return ['Error', 'Workflow of type %s not implemented or allowed' % workflow_type, '']

        # Check to see if intermediate purification is needed (for most workflows it is needed pre-HPLC)
        return_statement, additional_operations, additional_containers = intermediate_purification(grouping_information,
                                                                                                   workflow_details,
                                                                                                   wellplate_index,
                                                                                                   number_of_groups)
        if return_statement != 'Success':
            return ['Error', additional_operations]
        if additional_operations:
            operations_needed.extend(additional_operations)
        for container in additional_containers:
            queue_document_containers[container] = predefined_containers[container]
        operations_needed.append(['complete_queue'])

        return_statement, plate_operations = operations_list_to_details(workflow_type,
                                                                        operations_needed,
                                                                        workflow_details)
        if return_statement == 'Error':
            return ['Error', plate_operations, '']
        queues_return[reaction_queue_name]['operations'] = plate_operations
        if additional_plate_operations:
            if additional_plate_operations[-1] != ['complete_queue']:
                additional_plate_operations.append(['complete_queue'])
            return_statement, plate_operations = operations_list_to_details(workflow_type,
                                                                            additional_plate_operations,
                                                                            workflow_details)
            if return_statement == 'Error':
                return ['Error', plate_operations, '']
            if not additional_queue_name:
                return ['Error', 'No additional queue name specified']
            queues_return[additional_queue_name]['operations'] = plate_operations

        if additional_queue_name:
            pass
            # pprint(queues_return[additional_queue_name]['operations'])

    return_statement, queue_dependencies = find_queue_dependencies(queues_return, workflow_details)
    if return_statement == 'Error':
        return ['Error', queue_dependencies]
    for queue_document_name in queue_dependencies.keys():
        queues_return[queue_document_name]['dependency'] = list(queue_dependencies[queue_document_name])
    return ['Success', queues_return, reagents_needed]
