# -*- coding: utf-8 -*-
"""
Created on 05/18/22


@author: Brent Koscher
"""

import glob
import sys
import re
import json
import itertools
import yaml
import os
import time
import traceback

from pprint import pprint
from copy import deepcopy
from rdkit import Chem
import numpy as np

from functions.supporting_functions import print_progress_bar
from functions.supporting_functions import directory_check
from molecular_generation.rationale_model.generate_mols import run_rationale_model
from molecular_generation.site_mapping import atom_map_scaffold
from reaction_planning.api_client import APIClient
from functions.sascore.sascorer_class import SAScorer


def generate_molecules(job_information, job_name):
    # First we need to check to see if we need to generate molecules
    molecular_generation_details = job_information['molecular_generation']
    if 'generate_molecules' not in molecular_generation_details.keys():
        return ['Error', 'Molecular generation details not present: %s' % job_name]
    if not molecular_generation_details['generate_molecules']:
        return ['Success', '']

    # Next we can figure out if we have the keys we need
    if 'model_type' not in molecular_generation_details.keys():
        return ['Error', 'No molecular generation type specified: %s' % job_name]
    generation_method = molecular_generation_details['model_type']

    # Now we can make a molecular generation folder
    molecular_generation_folder = r'.\output_jobs\%s\molecular_generation' % job_name
    return_statement, print_statement = directory_check(molecular_generation_folder)
    if return_statement != 'Success':
        return ['Error', 'Unknown error occurred when making folders']

    # Now select the correct generation method to use
    if generation_method == 'rationale':
        return_statement, return_object = graph_completion_model(molecular_generation_folder,
                                                                 molecular_generation_details,
                                                                 job_name)
    elif generation_method == 'iterative_rationale':
        return_statement, return_object = iterative_graph_completion(molecular_generation_folder,
                                                                     molecular_generation_details,
                                                                     job_name)
    else:
        return ['Error', 'Molecular generation of type %s not implemented' % generation_method]

    if return_statement != 'Success':
        return ['Error', return_object]
    else:
        return ['Success', return_object]


def graph_completion_model(mol_gen_folder, mol_gen_details, job_name):
    # Start by preparing the molecular generation status file
    generation_status_filepath = r'%s\%s_generation_status.yaml' % (mol_gen_folder, job_name)
    if os.path.exists(generation_status_filepath):
        with open(generation_status_filepath, 'r') as yaml_file:
            mol_gen_status = yaml.load(yaml_file, Loader=yaml.Loader)
    else:
        mol_gen_status = {'scaffolds_to_run': [], 'scaffolds_completed': []}

    # Need to check for un-run and new scaffolds
    incoming_scaffolds = set(mol_gen_details['scaffolds'])
    mol_gen_status['scaffolds_to_run'] = incoming_scaffolds.difference(set(mol_gen_status['scaffolds_completed']))
    with open(generation_status_filepath, 'w') as yaml_file:
        yaml.dump(mol_gen_status, yaml_file)

    # Are there more scaffolds to run?
    if len(mol_gen_status['scaffolds_to_run']) == 0:
        return['Success', '']
    else:
        print('\t\tNeed to generate molecules for %s scaffolds' % len(mol_gen_status['scaffolds_to_run']))

    # Now we can prepare the input file for the graph completion model
    mol_gen_input = r'%s\%s_mol_gen_input.txt' % (mol_gen_folder, job_name)
    with open(mol_gen_input, 'w') as outfile:
        for scaffold in mol_gen_status['scaffolds_to_run']:
            outfile.write('%s %s\n' % (scaffold, scaffold))

    # Track down an output filepath that is available
    file_index = 1
    while True:
        mol_gen_output = r'%s\%s_mol_gen_output_%s.txt' % (mol_gen_folder, job_name, file_index)
        if os.path.exists(mol_gen_output):
            file_index += 1
        else:
            break

    # Now we can run the graph completion model
    return_statement, return_object = run_rationale_model(os.path.abspath(mol_gen_input),
                                                          os.path.abspath(mol_gen_output),
                                                          mol_gen_details['num_decodes'])
    if return_statement != 'Success':
        return ['Error', return_object]
    # Update the molecular generation status file
    with open(generation_status_filepath, 'r') as yaml_file:
        mol_gen_status = yaml.load(yaml_file, Loader=yaml.Loader)
    mol_gen_status['scaffolds_completed'] = list(set(mol_gen_status['scaffolds_completed']).union(incoming_scaffolds))
    with open(generation_status_filepath, 'w') as yaml_file:
        yaml.dump(mol_gen_status, yaml_file)

    # Now we can gather the generated molecules
    # The expected input from the files is based on the graph completion model output
    # Expected format: SCAFFOLD, GENERATED_SMILES/n
    output_files = glob.glob(r'%s\%s_mol_gen_output_*.txt' % (mol_gen_folder, job_name))
    smiles_set = set()
    for output_file in output_files:
        with open(output_file, 'r') as infile:
            for line in infile:
                cleaned_line = line.strip().split(',')
                if len(cleaned_line) != 2 or any([item in cleaned_line[0] for item in ['None', '.']]):
                    continue
                try:
                    rdkit_mol = Chem.MolFromSmiles(cleaned_line[1].strip())
                    rdkit_smiles = Chem.MolToSmiles(rdkit_mol)
                except:
                    continue
                smiles_set.add(rdkit_smiles)
    print('\t\tGathered %s generated SMILES for "%s"' % (len(smiles_set), job_name))
    gathered_molecules_path = r'%s\%s_gen_smiles.txt' % (mol_gen_folder, job_name)
    with open(gathered_molecules_path, 'w') as outfile:
        for item in list(smiles_set):
            outfile.write('%s\n' % item)

    return ['Success', gathered_molecules_path]


def iterative_graph_completion(generation_folder, generation_details, job_name, scaffold_filepath=None):
    # Start by preparing the molecular generation status file
    generation_status_filepath = r'%s\%s_generation_status.yaml' % (generation_folder, job_name)
    if os.path.exists(generation_status_filepath):
        with open(generation_status_filepath, 'r') as yaml_file:
            mol_gen_status = yaml.load(yaml_file, Loader=yaml.Loader)
    else:
        mol_gen_status = {'scaffolds_to_run': [], 'scaffolds_completed': []}

    # Need to check for un-run and new scaffolds
    incoming_scaffolds = set(generation_details['scaffolds'])
    mol_gen_status['scaffolds_to_run'] = incoming_scaffolds.difference(set(mol_gen_status['scaffolds_completed']))
    with open(generation_status_filepath, 'w') as yaml_file:
        yaml.dump(mol_gen_status, yaml_file)

    # Are there more scaffolds to run?
    if len(mol_gen_status['scaffolds_to_run']) == 0:
        return ['Success', '']
    else:
        print('\t\tNeed to generate molecules for %s scaffolds' % len(mol_gen_status['scaffolds_to_run']))

    number_of_iterations = generation_details['generation_iterations']
    for iteration_number in range(0, number_of_iterations):
        print('\tWorking on iteration %s' % iteration_number)
        output_directory = r'%s\%s_%s' % (generation_folder, job_name, iteration_number)

        return_statement, print_statement = directory_check(output_directory)
        if return_statement != 'Success':
            return ['Error', 'Unknown error occurred when making folders']

        if iteration_number == 0:
            current_scaffolds = {'scaffolds': mol_gen_status['scaffolds_to_run']}
        else:
            with open(scaffold_filepath, 'r') as input_file:
                current_scaffolds = {'scaffolds': []}
                for line in input_file:
                    cleaned_line = line.strip().split(',')
                    if len(cleaned_line) != 2 or any([item in cleaned_line[0] for item in ['None', '.']]):
                        continue
                    try:
                        rdkit_mol = Chem.MolFromSmiles(cleaned_line[1].strip())
                        rdkit_smiles = Chem.MolToSmiles(rdkit_mol)
                    except:
                        continue
                    current_scaffolds['scaffolds'].append(rdkit_smiles)

        return_statement, mapped_scaffolds_to_use = atom_map_scaffold(current_scaffolds['scaffolds'],
                                                                      generation_details)
        if return_statement != 'Success':
            return ['Error', 'Atom mapping failed to map scaffolds: %s' % job_name]

        file_index = 1
        while True:
            scaffold_input_filepath = r'%s\scaffolds_%s_input.txt' % (output_directory, file_index)
            if not os.path.exists(scaffold_input_filepath):
                break
            file_index += 1
            if file_index > 100:
                return ['Error', 'File name issue in %s' % output_directory]

        with open(scaffold_input_filepath, 'w') as outfile:
            for mapped_scaffold in mapped_scaffolds_to_use:
                outfile.write('%s %s\n' % (mapped_scaffold, mapped_scaffold))
        if 'filtering_method' in generation_details.keys():
            scaffold_filepath = scaffold_input_filepath.replace('_input.txt', '_intermediate.txt')
        else:
            scaffold_filepath = scaffold_input_filepath.replace('_input.txt', '_output.txt')

        return_statement, return_object = run_rationale_model(os.path.abspath(scaffold_input_filepath),
                                                              os.path.abspath(scaffold_filepath),
                                                              generation_details['num_decodes'])
        if return_statement != 'Success':
            return ['Error', return_object]

        # Now we can determine what type of filtering that we want to do
        if 'filtering_method' in generation_details.keys():
            filtered_smiles = []
            with open(scaffold_filepath, 'r') as input_file:
                generated_smiles = []
                for line in input_file:
                    cleaned_line = line.strip().split(',')
                    if len(cleaned_line) != 2 or any([item in cleaned_line[1] for item in ['None', '.']]):
                        continue
                    try:
                        rdkit_mol = Chem.MolFromSmiles(cleaned_line[1].strip())
                        rdkit_smiles = Chem.MolToSmiles(rdkit_mol)
                    except:
                        continue
                    generated_smiles.append(rdkit_smiles)
            if generation_details['filtering_method'] == 'askcos':
                askcos_output_directory = r'%s\askcos_data' % output_directory
                return_statement, filtered_smiles = get_reaction_trees(generated_smiles, askcos_output_directory)
                if return_statement != 'Success':
                    return ['Error', filtered_smiles]
            elif generation_details['filtering_method'] == 'sascore':
                if 'filtering_percentile' not in generation_details.keys():
                    filter_percentile = 0.5
                else:
                    filter_percentile = generation_details['filtering_percentile']
                return_statement, filtered_smiles = get_sa_score(generated_smiles, percentile_cutoff=filter_percentile)
                if return_statement != 'Success':
                    return ['Error', filtered_smiles]
            else:
                return ['Error', 'Iterative filtering method %s not defined' % generation_details['filtering_method']]

            scaffold_filepath = scaffold_input_filepath.replace('_input.txt', '_output.txt')
            with open(scaffold_filepath, 'w') as output_file:
                for filtered_item in filtered_smiles:
                    output_file.write('%s, %s\n' % ('_', filtered_item))

    compiled_generated_molecules = set()
    if number_of_iterations > 0:
        generation_files = glob.glob(r'%s\%s_*\*_output.txt' % (generation_folder, job_name))
        print(generation_files)
        for generation_file in generation_files:
            print('\t%s' % generation_file)
            with open(generation_file, 'r') as incoming_file:
                for line in incoming_file:
                    cleaned_line = line.strip().split(',')
                    if len(cleaned_line) != 2 or any([item in cleaned_line[0] for item in ['None', '.']]):
                        continue
                    try:
                        rdkit_mol = Chem.MolFromSmiles(cleaned_line[1].strip())
                        rdkit_smiles = Chem.MolToSmiles(rdkit_mol)
                    except:
                        continue
                    compiled_generated_molecules.add(rdkit_smiles)

    gathered_molecules_path = r'%s\%s_gen_smiles.txt' % (generation_folder, job_name)
    with open(gathered_molecules_path, 'w') as outfile:
        for item in list(compiled_generated_molecules):
            outfile.write('%s\n' % item)
    # Update the molecular generation status file
    with open(generation_status_filepath, 'r') as yaml_file:
        mol_gen_status = yaml.load(yaml_file, Loader=yaml.Loader)
    mol_gen_status['scaffolds_completed'] = list(set(mol_gen_status['scaffolds_completed']).union(incoming_scaffolds))
    with open(generation_status_filepath, 'w') as yaml_file:
        yaml.dump(mol_gen_status, yaml_file)

    return ['Success', '']


def get_reaction_trees(incoming_smiles, askcos_outpath):
    askcos_url = r'ASKCOS url and port #'
    askcos_client = APIClient(askcos_url, verify=False)
    smiles_to_request = incoming_smiles
    size_of_groups = 4
    request_groups = [smiles_to_request[i * size_of_groups:(i + 1) * size_of_groups] for
                      i in range((len(smiles_to_request) + size_of_groups - 1) // size_of_groups)]
    print('\tNeed to look at %s target SMILES in %s request groups' % (len(smiles_to_request), len(request_groups)))
    with open(r'.\molecular_generation\iterative_askcos_params.yaml', 'r') as yaml_file:
        askcos_details = yaml.load(yaml_file, Loader=yaml.Loader)
    max_askcos_tries = 3
    if not os.path.exists(askcos_outpath):
        os.makedirs(askcos_outpath)
    for product_index, request_group in enumerate(request_groups):
        print('\tWorking on reaction request group %s of %s' % (product_index + 1, len(request_groups)))
        current_ids = []
        for smiles in request_group:
            askcos_params = deepcopy(askcos_details['askcos_params']['tree_builder'])
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
                print('\t\t', task_result['state'], current_ids[current_index][0])
            except:
                continue
            try:
                if task_result['state'] == 'SUCCESS':
                    returned_trees = []
                    if 'output' not in task_result.keys():
                        print('\tASKCOS failed to return something...')
                        continue
                    if len(task_result['output']) != 0:
                        print('\t\t%s returned %s trees' % (current_ids[current_index][0],
                                                            len(task_result['output'])))
                    save_dataset[current_ids[current_index][0]] = task_result
            except Exception as e:
                print(e)
                continue

        file_index = 0
        while os.path.exists(r'%s\trees_%s.json' % (askcos_outpath, file_index)):
            file_index += 1
        file_path = r'%s\trees_%s.json' % (askcos_outpath, file_index)
        with open(file_path, 'w') as outfile:
            json.dump(save_dataset, outfile)

    files = glob.glob(r'%s\*.json' % askcos_outpath)
    smiles_with_paths = set()
    for file in files:
        with open(file, 'r') as jsonfile:
            incoming_dataset = json.load(jsonfile)
        for molecule_key in incoming_dataset.keys():
            if len(incoming_dataset[molecule_key]['output']) == 0:
                continue
            else:
                smiles_with_paths.add(molecule_key)
    return ['Success', smiles_with_paths]


def get_sa_score(incoming_smiles, percentile_cutoff=0.5):
    sa_scorer = SAScorer()

    mean = 2.23044
    sigma = 0.6526

    all_sa_scores = []
    for product_index, product_smiles in enumerate(incoming_smiles):
        try:
            rdkit_mol = Chem.MolFromSmiles(product_smiles)
            rdkit_smiles = Chem.MolToSmiles(rdkit_mol)
        except:
            continue
        sa_score = sa_scorer.calculateScore(Chem.MolFromSmiles(product_smiles))
        if sa_score < mean:
            sa_score_value = 1
        else:
            sa_score_value = np.exp(-(sa_score - mean) ** 2 / (2 * sigma))
        all_sa_scores.append([product_smiles, sa_score_value])

    all_sa_scores.sort(key=lambda x: x[1], reverse=True)
    last_index = int(np.ceil(len(all_sa_scores) * percentile_cutoff))
    remaining_products = all_sa_scores[0:last_index]
    pprint(remaining_products)
    return_products = set()
    for remaining_product in remaining_products:
        return_products.add(remaining_product[0])

    return ['Success', return_products]

