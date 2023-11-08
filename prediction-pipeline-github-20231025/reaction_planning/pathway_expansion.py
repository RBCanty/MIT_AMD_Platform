
import json
import yaml
import os

from functions.supporting_functions import directory_check
from functions.supporting_functions import print_progress_bar
from chemical_library.chemical_inventory_functions import load_reagents_datasets
from reaction_planning.template_expansion.askcos_template_expansion import askcos_template_expansion


def reaction_pathway_expansion(job_information, job_name):
    new_reaction_plans = False
    expansion_details = job_information['pathway_expansion']
    if 'pathway_expansion' not in expansion_details.keys() or not expansion_details['pathway_expansion']:
        return ['Success', 'Reaction pathway expansion not selected', new_reaction_plans]
    expansion_method = expansion_details['expansion_method']
    if expansion_method == 'templates':
        reaction_planning_folder = r'.\output_jobs\%s\reaction_planning' % job_name
        expansion_folder = r'%s\pathway_expansion' % reaction_planning_folder
        return_statement, return_object = directory_check(expansion_folder)
        if return_statement != 'Success':
            return ['Error', return_object, new_reaction_plans]
        return_statement, return_object, new_reaction_plans = askcos_template_expansion(job_information, job_name,
                                                                                        reaction_planning_folder)
        if return_statement != 'Success':
            return ['Error', return_object, new_reaction_plans]
    else:
        return ['Error', 'Pathway expansion method %s not specified' % expansion_method, new_reaction_plans]
    return ['Success', '', new_reaction_plans]


