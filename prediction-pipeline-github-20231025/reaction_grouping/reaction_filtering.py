
import json
import yaml
import os

from copy import deepcopy

from chemical_library.chemical_inventory_functions import load_reagents_datasets
from functions.supporting_functions import print_progress_bar


def reaction_filtering(job_information, job_name, reactions_filepath, new_reaction_plans):
    grouping_details = job_information['grouping_information']
    if 'pre_filtering' not in grouping_details.keys() or not grouping_details['pre_filtering']:
        return ['Success', 'No pre-filtering needed', reactions_filepath]
    filtered_filepath = r'.\output_jobs\%s\reaction_planning\%s_filtered_reactions.json' % (job_name, job_name)
    if os.path.exists(filtered_filepath) and not new_reaction_plans:
        return ['Success', 'Pre-filtering finished', filtered_filepath]

    # Load the reagents that we have access to
    reagent_details = job_information['reagent_information']
    if reagent_details['reagent_origin'] == 'restricted':
        price_limit = reagent_details['price_cutoff']
        return_statement, reagent_library, chemistry_inventory = load_reagents_datasets(reagent_details,
                                                                                        price_limit=price_limit)
        if return_statement != 'Success':
            return ['Error', reagent_library, '']
        available_chemicals = set(list(reagent_library.keys()))
    elif reagent_details['reagent_origin'] == 'unrestricted':
        reagent_library = {}
        available_chemicals = set()
    else:
        return ['Error', 'Unrecognized reagent source: %s' % reagent_details['reagent_origin'], '']

    # Load the reactions that we have available to consider
    allowed_temperatures = grouping_details['allowed_temperatures']
    with open(reactions_filepath, 'r') as json_file:
        reaction_trees = json.load(json_file)
    filtered_reaction_trees = {}
    for product_index, product_key in enumerate(reaction_trees.keys()):
        print_progress_bar(product_index + 1, len(list(reaction_trees.keys())), prefix='\t\t\tWorking: ', length=50)
        for tree_key in reaction_trees[product_key].keys():
            tree_details = reaction_trees[product_key][tree_key]
            cleaned_tree_details = {}

            # Start by checking the reactants if reagents are restricted
            if reagent_details['reagent_origin'] == 'restricted':
                tree_products = set()
                tree_reactants = set()
                for reaction_key in tree_details.keys():
                    tree_products.update(tree_details[reaction_key]['products'])
                    tree_reactants.update(tree_details[reaction_key]['reactants'])
                tree_reactants = tree_reactants.difference(tree_products)
                if not tree_reactants.issubset(available_chemicals):
                    continue

            # If the reagents are available then we can check the reaction conditions
            # Right now the check is for the reaction temperatures but other checks can be added
            for reaction_key in tree_details.keys():
                allowed_contexts = {}
                reaction_details = tree_details[reaction_key]
                for context_key in tree_details[reaction_key]['contexts'].keys():
                    reaction_context = tree_details[reaction_key]['contexts'][context_key]
                    if not (min(allowed_temperatures) <= reaction_context['temperature'] <= max(allowed_temperatures)):
                        continue
                    else:
                        if reaction_key not in cleaned_tree_details.keys():
                            cleaned_tree_details[reaction_key] = deepcopy(reaction_details)
                            cleaned_tree_details[reaction_key]['contexts'] = {}
                        cleaned_tree_details[reaction_key]['contexts'][context_key] = reaction_context
                if reaction_key not in cleaned_tree_details:
                    break

            # Now if any of the reactions were not valid we will ignore the whole tree
            if len(list(tree_details.keys())) != len(list(cleaned_tree_details.keys())):
                continue
            if product_key not in filtered_reaction_trees:
                filtered_reaction_trees[product_key] = {}
            filtered_reaction_trees[product_key][tree_key] = cleaned_tree_details

    with open(filtered_filepath, 'w') as json_file:
        json.dump(filtered_reaction_trees, json_file)
    return ['Success', 'Filtered reaction trees', filtered_filepath]
