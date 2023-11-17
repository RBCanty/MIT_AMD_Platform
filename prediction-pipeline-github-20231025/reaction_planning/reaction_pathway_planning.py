# -*- coding: utf-8 -*-
"""
Created on 11/28/2022

@author: Brent Koscher
"""

import os
import yaml
import traceback
import json

from functions.supporting_functions import directory_check
from reaction_planning.askcos_functions import get_reaction_contexts, get_reaction_trees, compile_reaction_details


def plan_reaction_pathways(job_information, job_name, smiles_filepath):
    # Which type of reaction planning to use?
    reaction_planning_details = job_information['reaction_planning']
    if 'planning_type' not in reaction_planning_details.keys():
        return ['Error', 'Reaction planning type was not specified', '']
    planning_type = reaction_planning_details['planning_type']

    if planning_type == 'askcos':
        reaction_planning_folder = r'.\output_jobs\%s\reaction_planning' % job_name
        return_statement, print_statement = directory_check(reaction_planning_folder)
        if return_statement != 'Success':
            return ['Error', print_statement, '']
        # Create/open a status file to keep track of ASKCOS planning
        askcos_tracker_filepath = r'%s\%s_status.yaml' % (reaction_planning_folder, job_name)
        if os.path.exists(askcos_tracker_filepath):
            with open(askcos_tracker_filepath, 'r') as yaml_file:
                askcos_status = yaml.load(yaml_file, Loader=yaml.Loader)
        else:
            askcos_status = {
                'molecules_to_run': [],
                'completed_molecules': [],
                'no_trees': [],
                'failed_molecules': [],
                'unique_reactions': [],
                'completed_contexts': []}

        # With the ASKCOS status tracker present we can find API calls we should make
        with open(smiles_filepath, 'r') as yaml_file:
            smiles_details = yaml.load(yaml_file, Loader=yaml.Loader)
        molecules_to_run = set(askcos_status['molecules_to_run']).union(smiles_details['smiles_set'])
        molecules_to_run = molecules_to_run.difference(set(askcos_status['completed_molecules']))
        askcos_status['molecules_to_run'] = list(molecules_to_run)
        with open(askcos_tracker_filepath, 'w') as yaml_file:
            yaml.dump(askcos_status, yaml_file)

        # Now we can move onto making the ASKCOS API calls
        askcos_data_filepath = r'%s\askcos_data' % reaction_planning_folder
        return_statement, print_statement = directory_check(askcos_data_filepath)
        if return_statement != 'Success':
            return ['Error', print_statement, '']
        max_askcos_iterations = 3
        new_askcos_information = 0
        for iteration in range(0, max_askcos_iterations):
            print('\t\tWorking on ASKCOS iteration %s of %s' % (iteration+1, max_askcos_iterations))
            if iteration == max_askcos_iterations-1:
                askcos_cleanup = True
            else:
                askcos_cleanup = False
            return_statement, number_of_products, askcos_status = get_reaction_trees(askcos_status,
                                                                                     job_information,
                                                                                     reaction_planning_folder,
                                                                                     job_name,
                                                                                     askcos_cleanup)
            if return_statement != 'Success':
                continue
            else:
                new_askcos_information += number_of_products
                print('\t\tNeeded to get trees for %s products' % number_of_products)
            if number_of_products == 0:
                print('\t\tNo additional molecules to plan')
                break
        # Rarely there are instances that leave out some reaction contexts
        return_statement, number_of_contexts, askcos_status = get_reaction_contexts(askcos_status,
                                                                                    job_information,
                                                                                    reaction_planning_folder,
                                                                                    job_name)
        if return_statement == 'Success':
            print('\t\tNeeded to finish contexts for %s reactions' % number_of_contexts)
            new_askcos_information += number_of_contexts
        else:
            return ['Error', 'Something failed with ASKCOS contexts', '']

        # With the reaction plans from ASKCOS we can compile/collate the details
        compiled_output_path = r'%s\%s_compiled_askcos.json' % (reaction_planning_folder, job_name)
        if os.path.exists(compiled_output_path) and new_askcos_information == 0:
            print('\t\tCompiled ASKCOS data already exists and no new information present')
        else:
            reaction_trees = compile_reaction_details(askcos_data_filepath, job_information)
            with open(compiled_output_path, 'w') as json_file:
                json.dump(reaction_trees, json_file)

    else:
        return ['Error', 'Reaction planning of type %s is not implemented' % planning_type, '']

    return ['Success', '\t\tCompleted reaction planning using %s' % planning_type, new_askcos_information]

