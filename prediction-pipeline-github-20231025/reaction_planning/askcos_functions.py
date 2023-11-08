# -*- coding: utf-8 -*-
"""
Created on 01/21/22

Condensing of the ASKCOS related functions into this file
 -- get_reaction_trees
 -- get_reaction_contexts
 -- compile_reaction_details

@author: brent
"""

import glob
import re
import json
import itertools
from copy import deepcopy
import yaml
import os
import time
from pprint import pprint
from rdkit import Chem

from reaction_planning.api_client import APIClient
from reaction_planning.context_correction.context_correction import contexts_handler
from functions.supporting_functions import print_progress_bar


def extract_values(obj, key):
    """Pull all values of specified key from nested JSON."""
    arr = []
    def extract(obj, arr, key):
        """Recursively search for values of key in JSON tree."""
        if isinstance(obj, dict):
            for k, v in obj.items():
                if k == 'tforms':
                    arr.append(v)
                elif isinstance(v, (dict, list)):
                    extract(v, arr, key)
                elif k == key:
                    arr.append(v)
        elif isinstance(obj, list):
            for item in obj:
                extract(item, arr, key)
        return arr
    results = extract(obj, arr, key)
    return results


def get_reaction_trees(job_status_func, incoming_dict, file_directory, job_name, cleanup=False):
    askcos_url = r'https://askcos.mit.edu:7000/'
    askcos_client = APIClient(askcos_url, verify=False)
    smiles_to_request = list(set(job_status_func['molecules_to_run']).difference(set(job_status_func['completed_molecules'])))
    size_of_groups = 4
    request_groups = [smiles_to_request[i*size_of_groups:(i+1)*size_of_groups] for
                      i in range((len(smiles_to_request)+size_of_groups-1)//size_of_groups)]
    print('\t\tNeed to look at %s/%s target SMILES in %s request groups' % (len(smiles_to_request),
                                                                            len(job_status_func['molecules_to_run']),
                                                                            len(request_groups)))
    askcos_data_directory = r'%s\askcos_data' % file_directory
    max_askcos_tries = 3
    for product_index, request_group in enumerate(request_groups):
        print('\t\tWorking on reaction request group %s of %s' % (product_index+1, len(request_groups)))
        current_ids = []
        for smiles in request_group:
            askcos_params = deepcopy(incoming_dict['reaction_planning']['askcos_params']['tree_builder'])
            askcos_params['smiles'] = smiles
            for iterations in range(0, max_askcos_tries):
                try:
                    result = askcos_client.post('tree-builder', data=askcos_params)
                    current_ids.append([smiles, result['task_id']])
                    break
                except:
                    time.sleep(5)

        try:
            task_results = [askcos_client.get_result(tid[1], timeout=600, interval=5) for tid in current_ids]
        except json.JSONDecodeError:
            print('\tJSON DECODE ERROR!!! Problem with it...')
            continue
        save_dataset = {}
        for current_index, task_result in enumerate(task_results):
            try:
                if task_result['state'] == 'SUCCESS':
                    returned_trees = []
                    if 'output' not in task_result.keys():
                        print('\tASKCOS failed to return something...')
                        continue
                    for output in task_result['output']:
                        output_smiles = extract_values(output, 'smiles')
                        output_reactions = [item for item in output_smiles if '>>' in item]
                        current_unique_reactions = set(job_status_func['unique_reactions'])
                        current_unique_reactions.update(output_reactions)
                        job_status_func['unique_reactions'] = list(current_unique_reactions)
                    save_dataset[current_ids[current_index][0]] = task_result
                    job_status_func['completed_molecules'].append(current_ids[current_index][0])
                    if len(returned_trees) == 0:
                        job_status_func['no_trees'].append(current_ids[current_index][0])
                elif cleanup:
                    job_status_func['completed_molecules'].append(current_ids[current_index][0])
                    job_status_func['failed_molecules'].append(current_ids[current_index][0])
            except Exception as e:
                print(e)
                continue
        with open(file_directory + '\\%s_status.yaml' % job_name, 'w') as outfile:
            yaml.dump(job_status_func, outfile)
        file_index = 0
        while os.path.exists(askcos_data_directory + '\\%s_trees_%s.json' % (job_name, file_index)):
            file_index += 1
        file_path = askcos_data_directory + '\\%s_trees_%s.json' % (job_name, file_index)
        with open(file_path, 'w') as outfile:
            json.dump(save_dataset, outfile)
        return_statement, number_of_contexts, job_status_func = get_reaction_contexts(job_status_func,
                                                                                      incoming_dict,
                                                                                      file_directory,
                                                                                      job_name)
        if return_statement != 'Error':
            print('\t\t\tNeeded to get %s reaction contexts from ASKCOS' % number_of_contexts)
    return ['Success', len(smiles_to_request), job_status_func]


def get_reaction_contexts(job_status_func, incoming_dict, file_directory, job_name):
    askcos_url = r'https://askcos.mit.edu:7000/'
    askcos_client = APIClient(askcos_url, verify=False)
    completed_contexts = set(job_status_func['completed_contexts'])
    contexts_to_request = list(set(job_status_func['unique_reactions']).difference(completed_contexts))
    planning_details = incoming_dict['reaction_planning']
    if len(contexts_to_request) > 0:
        size_of_groups = 4
        context_groups = [contexts_to_request[i*size_of_groups:(i+1)*size_of_groups] for 
                          i in range((len(contexts_to_request)+size_of_groups-1)//size_of_groups)]
        print('\t\t\tNeed to request contexts for %s reactions (%s groups)' % (len(contexts_to_request),
                                                                             len(context_groups)))
        save_dataset = {}
        for context_index, context_group in enumerate(context_groups):
            print('\t\t\tWorking on context group %s of %s' % (context_index + 1, len(context_groups)))
            current_ids = []
            for reaction_smiles in context_group:
                askcos_params = deepcopy(planning_details['askcos_params']['context_builder'])
                askcos_params['reactants'] = reaction_smiles.split('>>')[0]
                askcos_params['products'] = reaction_smiles.split('>>')[1]
                result = askcos_client.post('context', data=askcos_params)
                current_ids.append([reaction_smiles, result['task_id']])
            try:
                task_results = [askcos_client.get_result(tid[1], timeout=600, interval=1) for tid in current_ids]
            except json.JSONDecodeError:
                print('\tJSON DECODE ERROR!!! Problem with it...')
                continue
            for current_index, task_result in enumerate(task_results):
                #try:
                if 'context_correction' in planning_details.keys() and planning_details['context_correction']:
                    cleaned_contexts = contexts_handler(current_ids[current_index][0], task_result)
                    # print(cleaned_contexts)
                else:
                    cleaned_contexts = task_result
                save_dataset[current_ids[current_index][0]] = cleaned_contexts
                job_status_func['completed_contexts'].append(current_ids[current_index][0])
                #except:
                #    continue

            askcos_data_directory = r'%s\askcos_data' % file_directory
            with open(file_directory + '\\%s_status.yaml' % job_name, 'w') as outfile:
                yaml.dump(job_status_func, outfile)
            file_index = 0
            while os.path.exists(askcos_data_directory + '\\%s_contexts_%s.json' % (job_name, file_index)):
                file_index += 1
            file_path = askcos_data_directory + '\\%s_contexts_%s.json' % (job_name, file_index)
            with open(file_path, 'w') as outfile:
                json.dump(save_dataset, outfile)
    else:
        print('\t\tNo additional context requests needed!')
    return ['Success', len(contexts_to_request), job_status_func]
   
   
def compile_reaction_details(details_filepath, incoming_parameters):
    compiled_reaction_details = {}
    reaction_contexts = {}

    template_context_corrections_path = r'.\reaction_planning\context_correction\template_corrections.yaml'
    with open(template_context_corrections_path, 'r') as yaml_file:
        template_context_corrections = yaml.load(yaml_file, Loader=yaml.Loader)
    context_files = glob.glob(details_filepath + '\\*_contexts_*.json')
    for file_index, output_file in enumerate(context_files):
        print_progress_bar(file_index + 1, len(context_files), prefix='\t\tProcessing:', suffix='Complete', length=50)
        if not re.search(r'_contexts_\d+.json', output_file):
            continue
        with open(output_file, 'r') as jsonfile:
            incoming_reaction_contexts = json.load(jsonfile)
            reaction_contexts = {**reaction_contexts, **incoming_reaction_contexts}
    tree_files = glob.glob(details_filepath + '\\*_trees_*.json')
    for file_index, output_file in enumerate(tree_files):
        print(output_file)
        print_progress_bar(file_index + 1, len(tree_files), prefix='\t\tProcessing:', suffix='Complete', length=50)
        if not re.search(r'_trees_\d+.json', output_file):
            continue
        with open(output_file, 'r') as jsonfile:
            incoming_reaction_trees = json.load(jsonfile)
        for product_key in incoming_reaction_trees.keys():
            if len(incoming_reaction_trees[product_key]['output']) == 0:
                print('\t', product_key)
                continue
            for tree_index, reaction_tree in enumerate(incoming_reaction_trees[product_key]['output']):
                if tree_index >= incoming_parameters['grouping_information']['tree_rank_limit']:
                    break
                raw_smiles_set = extract_values(reaction_tree, 'smiles')
                smiles_set = []
                template_set = []
                for smiles_index, raw_smiles in enumerate(raw_smiles_set):
                    if '>>' in raw_smiles:
                        smiles_set.append(raw_smiles)
                        template_set.append(raw_smiles_set[smiles_index + 1])
                smiles_set = smiles_set[::-1]
                template_set = template_set[::-1]
                if len(smiles_set) > incoming_parameters['grouping_information']['tree_depth_limit']:
                    continue
                if len(smiles_set) == 0:
                    continue
                if product_key not in compiled_reaction_details.keys():
                    compiled_reaction_details[product_key] = {}
                tree_key_name = '%s_Tree_%s' % (product_key, tree_index + 1)
                # Here we need to determine the reaction dependence of things
                number_of_reactions = len(smiles_set)
                reaction_pairs = itertools.combinations(list(range(0, number_of_reactions)), 2)
                dependencies = {}
                for reaction_pair in reaction_pairs:
                    first_reaction = smiles_set[reaction_pair[0]].split('>>')
                    first_reaction_products = set(first_reaction[1].split('.'))
                    second_reaction = smiles_set[reaction_pair[1]].split('>>')
                    second_reaction_reactants = set(second_reaction[0].split('.'))
                    if len(list(first_reaction_products.intersection(second_reaction_reactants))) != 0:
                        if reaction_pair[1] not in dependencies.keys():
                            dependencies[reaction_pair[1]] = set()
                        dependencies[reaction_pair[1]].add(reaction_pair[0])
                compiled_reaction_details[product_key][tree_key_name] = {}
                problem_reaction = 0
                additional_reaction_index = 0
                for reaction_index, reaction in enumerate(smiles_set):
                    try:
                        sub_reaction_context = reaction_contexts[reaction]
                    except KeyError:
                        problem_reaction += 1
                        continue
                    reactants = reaction.split('>>')[0].split('.')
                    products = reaction.split('>>')[1].split('.')
                    if 'reagent_prep_reactions' in incoming_parameters.keys():
                        for reactant in reactants:
                            for special_reactant in incoming_parameters['reagent_prep_reactions'].keys():
                                if Chem.MolFromSmiles(reactant) == Chem.MolFromSmiles(special_reactant) or reactant == special_reactant:
                            #if reactant in incoming_parameters['reagent_prep_reactions'].keys():
                                    for reaction_key in incoming_parameters['reagent_prep_reactions'][special_reactant]['reactions'].keys():
                                        reaction_details = incoming_parameters['reagent_prep_reactions'][special_reactant]['reactions'][reaction_key]
                                        context_details = incoming_parameters['extra_template_conditions'][reaction_details['reaction_contexts']]
                                        subkey_name = 'reaction_%s' % str(additional_reaction_index + reaction_index + 1)
                                        reaction_dependencies = list(compiled_reaction_details[product_key][tree_key_name].keys())
                                        compiled_reaction_details[product_key][tree_key_name][subkey_name] = {'reaction': reaction_details['reaction_smiles'], 
                                                                                                              'reactants': reaction_details['reaction_smiles'].split('>>')[0].split('.'), 
                                                                                                              'products': reaction_details['reaction_smiles'].split('>>')[1].split('.'), 
                                                                                                              'contexts': context_details['context'], 
                                                                                                              'templates': [context_details['template']['reaction_smarts']], 
                                                                                                              'dependencies': list(compiled_reaction_details[product_key][tree_key_name].keys())}
                                        additional_reaction_index += 1

                    subkey_name = 'reaction_%s' % str(additional_reaction_index + reaction_index + 1)
                    reaction_dependencies = []

                    for previous_reaction_index in range(1, additional_reaction_index + reaction_index + 1):
                        previous_key = 'reaction_%s' % str(previous_reaction_index)
                        if subkey_name == previous_key:
                            break
                        previous_reaction_products = compiled_reaction_details[product_key][tree_key_name][previous_key]['products']
                        if len(list(set(previous_reaction_products).intersection(set(reactants)))) != 0:
                            reaction_dependencies.append(previous_key)

                    compiled_reaction_details[product_key][tree_key_name][subkey_name] = {'reaction': reaction, 'reactants': reactants,
                                        'products': products, 'contexts': {}, 'templates': template_set[reaction_index],
                                        'dependencies': reaction_dependencies}

                    for template in template_set[reaction_index]:
                        if template in template_context_corrections.keys():
                            context_correction = template_context_corrections[template]
                            break
                    else:
                        context_correction = {}
                        
                    current_reaction = compiled_reaction_details[product_key][tree_key_name][subkey_name]
                    if context_correction:
                        context_key = 'context_1'
                        current_reaction['contexts'][context_key] = context_correction
                    else:
                        for context_index, context in enumerate(sub_reaction_context['output']):
                            context_key = 'context_%s' % str(context_index + 1)
                            chemical_fields = ['reagent', 'solvent', 'catalyst']
                            for chemical_field in chemical_fields:
                                if context[chemical_field] == '':
                                    continue
                                if type(context[chemical_field]) == list:
                                    field_split = context[chemical_field]
                                else:
                                    field_split = context[chemical_field].split('.')
                                if context_key not in current_reaction['contexts'].keys():
                                    current_reaction['contexts'][context_key] = {'temperature': context['temperature']}
                                current_reaction['contexts'][context_key][chemical_field] = field_split
    return compiled_reaction_details
























