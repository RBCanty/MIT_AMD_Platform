
import yaml
import json
import sys
import itertools

import numpy as np
from pprint import pprint
from copy import deepcopy
from scipy.stats import lognorm
from rdkit import Chem

from chemical_library.chemical_inventory_functions import load_reagents_datasets
from functions.supporting_functions import print_progress_bar
from functions.sascore.sascorer_class import SAScorer


def group_products(job_information, job_name, chemprop_data_filepath, candidate_filepath):
    with open(candidate_filepath, 'r') as json_file:
        candidate_reactions = json.load(json_file)
    with open(chemprop_data_filepath, 'r') as yaml_file:
        candidate_predictions = yaml.load(yaml_file, Loader=yaml.Loader)

    temperature_sequence = None
    if job_information['grouping_information']['method'] == 'normal':
        statement, grouping_bins, molecule_values = evaluate_property_prediction_value(job_information,
                                                                                       candidate_predictions)
        if statement != 'Success':
            return ['Error', grouping_bins, '', '']

        # Load the reagents that we have access to
        reagent_details = job_information['reagent_information']
        price_limit = reagent_details['price_cutoff']
        statement, reagent_library, chemistry_inventory = load_reagents_datasets(reagent_details,
                                                                                 price_limit=price_limit)
        if statement != 'Success':
            return ['Error', reagent_library, '', '']

        # Pre-filter the reaction trees based on reagents available
        statement, candidate_reactions = filter_reaction_trees(job_information,
                                                               candidate_reactions,
                                                               reagent_library,
                                                               molecule_values=molecule_values)
        if statement != 'Success':
            return ['Error', candidate_reactions, '', '']

        statement, temperature_sequence, candidate_reactions = temperature_sequence_finding(job_information,
                                                                                            candidate_reactions,
                                                                                            molecule_values=molecule_values,
                                                                                            sequence_type='prefilter')
        if statement != 'Success':
            return ['Error', temperature_sequence, '', '']



        statement, reagents_needed, candidate_reactions = product_grouping_iterator(job_information,
                                                                                    grouping_bins,
                                                                                    candidate_reactions,
                                                                                    reagent_library,
                                                                                    temperature_sequence,
                                                                                    molecule_values=molecule_values,
                                                                                    candidate_predictions=candidate_predictions)
        if statement != 'Success':
            return ['Error', reagents_needed, '', '']

        statement, temperature_sequence, candidate_reactions = temperature_sequence_finding(job_information,
                                                                                            candidate_reactions,
                                                                                            molecule_values=molecule_values)
        if statement != 'Success':
            return ['Error', temperature_sequence, '', '']

    else:
        return ['Error', 'Grouping method was not specified in the input document', '', '']

    return ['Success', temperature_sequence, candidate_reactions, molecule_values]


def filter_reaction_trees(job_information, candidate_reactions, reagent_library, molecule_values=None):
    if molecule_values is None:
        return ['Success', candidate_reactions]

    missing_chemicals = {}
    filtered_candidate_reactions = {}
    available_chemicals = set(list(reagent_library.keys()))
    for product_key in candidate_reactions:
        for reaction_tree_key in candidate_reactions[product_key]:
            reaction_tree = candidate_reactions[product_key][reaction_tree_key]
            tree_reagents = set()
            tree_products = set()
            for reaction_key in reaction_tree.keys():
                reaction_smiles = reaction_tree[reaction_key]['reaction']
                tree_products.update(reaction_smiles.split('>>')[1].split('.'))
                tree_reagents.update(reaction_smiles.split('>>')[0].split('.'))
            tree_chemicals = tree_reagents.difference(tree_products)
            if len(tree_chemicals) == 0:
                continue
            if tree_chemicals.issubset(available_chemicals):
                if product_key not in filtered_candidate_reactions.keys():
                    filtered_candidate_reactions[product_key] = {}
                filtered_candidate_reactions[product_key][reaction_tree_key] = reaction_tree
            else:
                for tree_chemical in tree_chemicals:
                    if tree_chemical not in available_chemicals:
                        if tree_chemical not in missing_chemicals:
                            missing_chemicals[tree_chemical] = 0
                        missing_chemicals[tree_chemical] += 1
    print('-----------------------------')
    print('Missing chemicals')
    pprint(missing_chemicals)
    print('-----------------------------')

    grouping_information = job_information['grouping_information']
    predicted_values = []
    for product_key in molecule_values.keys():
        if product_key not in filtered_candidate_reactions.keys():
            continue
        predicted_values.append(molecule_values[product_key]['prediction_value'])
    if 'uncertainty' in grouping_information.keys():
        if grouping_information['uncertainty']['method'] == 'exploration':
            threshold = np.percentile(predicted_values, 50)
        elif grouping_information['uncertainty']['method'] == 'exploitation':
            threshold = np.percentile(predicted_values, 75)
        else:
            threshold = np.percentile(predicted_values, 25)
    else:
        threshold = np.percentile(predicted_values, 25)

    product_keys = list(filtered_candidate_reactions.keys())
    for product_key in product_keys:
        if product_key not in molecule_values:
            filtered_candidate_reactions.pop(product_key, None)
        elif molecule_values[product_key]['prediction_value'] < threshold:
            filtered_candidate_reactions.pop(product_key, None)

    return ['Success', filtered_candidate_reactions]


def product_grouping_iterator(job_information, grouping_bins, candidate_reactions, reagent_library, temperature_sequence,
                              molecule_values=None, candidate_predictions=None):
    if molecule_values is None:
        molecule_values = {}
    if candidate_predictions is None:
        candidate_predictions = {}
    grouping_information = job_information['grouping_information']
    temperature_spacing = grouping_information['temperature_spacing']
    sa_scorer = SAScorer()

    # First we need to figure out which reactants are used in the reaction trees
    chemicals_needed = {}
    for grouping_bin in grouping_bins.keys():
        chemicals_needed[grouping_bin] = {}
        for product_key in grouping_bins[grouping_bin]['molecules']:
            # TODO Check for a previous product
            if product_key not in candidate_reactions:
                continue
            product_trees = candidate_reactions[product_key]
            for reaction_tree_key in product_trees.keys():
                reaction_tree = product_trees[reaction_tree_key]
                tree_products = set()
                tree_reagents = set()
                for reaction_key in reaction_tree.keys():
                    reaction_smiles = reaction_tree[reaction_key]['reaction']
                    tree_products.update(reaction_smiles.split('>>')[1].split('.'))
                    tree_reagents.update(reaction_smiles.split('>>')[0].split('.'))
                tree_chemicals = set(tree_reagents).difference(tree_products)
                for tree_chemical in tree_chemicals:
                    if tree_chemical not in chemicals_needed[grouping_bin].keys():
                        chemicals_needed[grouping_bin][tree_chemical] = set()
                    chemicals_needed[grouping_bin][tree_chemical].add(product_key)

    chemicals_added = {'chemicals': set(), 'budget': job_information['reagent_information']['reagents_budget'],
                       'money_spent': 0}
    bins_to_workon = list(grouping_bins.keys())
    for iteration_number in range(grouping_information['max_grouping_iterations']):
        print('Working on iteration %s' % (iteration_number+1))
        for grouping_bin in grouping_bins:
            if grouping_bin not in bins_to_workon:
                continue
            print('\tLooking at bin %s' % grouping_bin)

            available_chemicals = deepcopy(chemicals_added['chemicals'])
            if grouping_bin not in chemicals_needed.keys():
                print('\tThere are no chemicals needed in the bin %s to group...' % grouping_bin)
                if grouping_bin in bins_to_workon:
                    bins_to_workon.remove(grouping_bin)
                continue
            chemicals_to_check = list(set(list(chemicals_needed[grouping_bin].keys())).difference(available_chemicals))
            if len(chemicals_to_check) == 0:
                print('\tNo chemicals left to check...')
                if grouping_bin in bins_to_workon:
                    bins_to_workon.remove(grouping_bin)
                continue
            number_of_chemicals_to_check = len(chemicals_to_check)
            chemical_values = []
            product_values_per_chemical = {}
            for chemical_index, chemical_to_check in enumerate(chemicals_to_check):
                print_progress_bar(chemical_index + 1, number_of_chemicals_to_check, prefix='\tWorking: ', length=50)
                if chemical_to_check not in reagent_library.keys():
                    continue
                vendor_information = reagent_library[chemical_to_check]
                iteration_chemical_inventory = deepcopy(available_chemicals)
                iteration_chemical_inventory.add(chemical_to_check)
                tree_progress_value = []
                new_products = []
                value_added = []
                for potential_target in list(chemicals_needed[grouping_bin][chemical_to_check]):
                    if molecule_values != {} and potential_target not in molecule_values.keys():
                        continue
                    target_reaction_trees = candidate_reactions[potential_target]
                    reagent_value = fraction_tree_progress(target_reaction_trees,
                                                           iteration_chemical_inventory,
                                                           chemical_to_check)
                    tree_progress_value.append(reagent_value)
                    if reagent_value == 1:
                        new_products.append(potential_target)
                    extra_value = []
                    prediction_value = molecule_values[potential_target]['prediction_value']
                    uncertainty_value = molecule_values[potential_target]['uncertainty_value']
                    if 'sa_score' in grouping_information.keys() and grouping_information['sa_score'][0]:
                        sa_score = sa_scorer.calculateScore(Chem.MolFromSmiles(potential_target))
                        extra_value.append(1 / sa_score)
                    value_added.append(np.prod(extra_value) * prediction_value * uncertainty_value) # * reagent_value)
                if len(value_added) == 0:
                    value_added.append((np.sum(tree_progress_value) / len(tree_progress_value) ** 2) * 1E-100)
                reagent_price = vendor_information['prices'][0][0]
                value_added.sort(reverse=True)
                if grouping_information['uncertainty']['method'] == 'exploitation':
                    chemical_value = value_added[0] * (1 / reagent_price)
                else:
                    chemical_value = np.sum(value_added[:10]) * (1 / reagent_price)
                chemical_values.append([chemical_to_check, float(chemical_value), reagent_price])
            chemical_values.sort(key=lambda x: x[1], reverse=True)
            pprint(chemical_values)
            if len(chemical_values) == 0:
                continue
            for chemical_value in chemical_values:
                if chemicals_added['money_spent'] + chemical_value[2] < chemicals_added['budget']:
                    chemicals_added['chemicals'].add(chemical_value[0])
                    chemicals_added['money_spent'] += chemical_value[2]
                    print('\t\tAdded %s to the reagent set ($%s of $%s)' % (chemical_value,
                                                                            chemicals_added['money_spent'],
                                                                            chemicals_added['budget']))
                    break
            else:
                print('\t\tThere was no chemical suitable that fits into the budget specified')
                continue

            products_made = set()
            for product_key in candidate_reactions.keys():
                if product_key not in molecule_values:
                    continue
                for reaction_tree_key in candidate_reactions[product_key].keys():
                    reaction_tree = candidate_reactions[product_key][reaction_tree_key]
                    tree_reagents = set()
                    tree_products = set()
                    for reaction_key in reaction_tree.keys():
                        reaction_smiles = reaction_tree[reaction_key]['reaction']
                        tree_products.update(reaction_smiles.split('>>')[1].split('.'))
                        tree_reagents.update(reaction_smiles.split('>>')[0].split('.'))
                    tree_chemicals = tree_reagents.difference(tree_products)
                    if len(tree_chemicals) == 0:
                        continue
                    if tree_chemicals.issubset(chemicals_added['chemicals']):
                        products_made.add(product_key)
                        break
            print('\t\tAble to access %s molecules' % len(products_made))
            if len(products_made) >= grouping_information['number_of_targets_per_group']:
                print('\t\tReached the terminal number of molecules to make')
                bins_to_workon.remove(grouping_bin)
        if len(bins_to_workon) == 0:
            break

    updated_reagent_set = deepcopy(chemicals_added['chemicals'])
    filtered_candidate_reactions = {}
    for product_key in candidate_reactions.keys():
        if product_key not in molecule_values:
            continue
        for reaction_tree_key in candidate_reactions[product_key].keys():
            reaction_tree = candidate_reactions[product_key][reaction_tree_key]
            tree_reagents = set()
            tree_products = set()
            for reaction_key in reaction_tree.keys():
                reaction_smiles = reaction_tree[reaction_key]['reaction']
                tree_products.update(reaction_smiles.split('>>')[1].split('.'))
                tree_reagents.update(reaction_smiles.split('>>')[0].split('.'))
            tree_chemicals = tree_reagents.difference(tree_products)
            if len(tree_chemicals) == 0:
                continue
            if tree_chemicals.issubset(updated_reagent_set):
                filtered_candidate_reactions[product_key] = {}
                filtered_candidate_reactions[product_key][reaction_tree_key] = reaction_tree
                print(product_key, molecule_values[product_key])
                pprint(candidate_predictions[product_key])
                break

    products_to_remove = set()
    for product_key in candidate_reactions.keys():
        reaction_tree_to_remove = set()
        for reaction_tree_key in candidate_reactions[product_key].keys():
            reaction_tree = candidate_reactions[product_key][reaction_tree_key]
            number_of_reactions = len(list(reaction_tree.keys()))
            tree_products = set()
            for reaction_key in reaction_tree.keys():
                reaction_smiles = reaction_tree[reaction_key]['reaction']
                tree_products.update(reaction_smiles.split('>>')[1].split('.'))
            for reaction_index, reaction_key in enumerate(reaction_tree.keys()):
                reaction_smiles = reaction_tree[reaction_key]['reaction']
                reaction_chemical_set = set(reaction_smiles.split('>>')[0].split('.')).difference(tree_products)
                if not reaction_chemical_set.issubset(updated_reagent_set):
                    reaction_tree_to_remove.add(reaction_tree_key)
                    break
                contexts_to_remove = []
                for context in reaction_tree[reaction_key]['contexts'].keys():
                    context_details = reaction_tree[reaction_key]['contexts'][context]
                    temperature_key = temperature_spacing * round(context_details['temperature'] / temperature_spacing)
                    if temperature_key not in temperature_sequence:
                        contexts_to_remove.append(context)
                for context_to_remove in contexts_to_remove:
                    reaction_tree[reaction_key]['contexts'].pop(context_to_remove, None)
                if len(reaction_tree[reaction_key]['contexts']) == 0:
                    reaction_tree_to_remove.add(reaction_tree_key)
                    break
        for tree_to_remove in list(reaction_tree_to_remove):
            candidate_reactions[product_key].pop(tree_to_remove, None)
            if len(list(candidate_reactions[product_key].keys())) == 0:
                products_to_remove.add(product_key)
    for product_to_remove in products_to_remove:
        candidate_reactions.pop(product_to_remove, None)

    with open(r'.\testing_candidate_reactions.json', 'w') as json_file:
        json.dump(candidate_reactions, json_file)
    with open(r'.\testing_filtered_candidate_reactions.json', 'w') as json_file:
        json.dump(filtered_candidate_reactions, json_file)
    with open(r'.\testing_molecule_values.json', 'w') as json_file:
        json.dump(molecule_values, json_file)

    return ['Success', updated_reagent_set, filtered_candidate_reactions]


def fraction_tree_progress(target_dict, available_chemicals, new_reagent_smiles):
    fractional_value_list = []
    for reaction_tree in target_dict.keys():
        tree_products = set()
        tree_reagents = set()
        for reaction_key in target_dict[reaction_tree].keys():
            reaction_smiles = target_dict[reaction_tree][reaction_key]['reaction']
            tree_products.update(reaction_smiles.split('>>')[1].split('.'))
            tree_reagents.update(list(set(reaction_smiles.split('>>')[0].split('.')).difference(tree_products)))
        if new_reagent_smiles not in tree_reagents:
            continue
        reactions_completed = set()
        for reaction_key in target_dict[reaction_tree].keys():
            reaction_smiles = target_dict[reaction_tree][reaction_key]['reaction']
            reaction_chemical_set = set(reaction_smiles.split('>>')[0].split('.')).difference(tree_products)
            if reaction_chemical_set.issubset(available_chemicals):
                reactions_completed.add(reaction_key)
        if len(list(tree_reagents)) == 0:
            continue
        fractional_value_list.append((len(reactions_completed) / len(list(target_dict[reaction_tree].keys()))))
    if len(fractional_value_list) == 0:
        fractional_value = 0.00001
    else:
        fractional_value = max(fractional_value_list)
    return fractional_value


def evaluate_property_prediction_value(job_information, candidate_predictions):
    grouping_information = job_information['grouping_information']
    model_information = job_information['chemprop_details']

    if 'scaffolds' not in grouping_information.keys():
        scaffolds = ['*']
    else:
        scaffolds = grouping_information['scaffolds']

    # To start we will gather the predicted values into a dictionary
    constraining_models = model_information['constraining_models']
    model_keys = list(constraining_models.keys())
    grouping_values = []
    for model_key in model_keys:
        centers = constraining_models[model_key]['center']
        ranges = constraining_models[model_key]['range']
        group_types = constraining_models[model_key]['type']
        if type(centers) != list:
            grouping_values.append([[centers, ranges, group_types]])
        else:
            if type(group_types) != list:
                group_types = [group_types]
            grouping_values.extend([zip(centers, ranges, group_types)])
    grouping_key_values = itertools.product(*grouping_values)

    grouping_bins = {}
    for scaffold in scaffolds:
        for grouping_key_value in grouping_key_values:
            values = [str(item[0]) for item in grouping_key_value]
            grouping_bin_key = '_'.join([scaffold] + values)
            grouping_bins[grouping_bin_key] = {'models': {}, 'scaffold': scaffold, 'molecules': {}}
            if scaffold != '*':
                grouping_bins[grouping_bin_key]['scaffold_mol'] = Chem.MolFromSmarts(scaffold)
            grouping_bin_models = grouping_bins[grouping_bin_key]['models']
            for model_index, model in enumerate(model_keys):
                grouping_bin_models[model] = {'grouping_detail': {'center': grouping_key_value[model_index][0],
                                                                  'width': grouping_key_value[model_index][1],
                                                                  'type': grouping_key_value[model_index][2]}}

    # We need to figure out which of the molecules fit into which bins
    for candidate_key in candidate_predictions.keys():
        model_predictions = candidate_predictions[candidate_key]
        candidate_mol = Chem.MolFromSmiles(candidate_key)
        for grouping_bin in grouping_bins.keys():
            bin_information = grouping_bins[grouping_bin]
            if bin_information['scaffold'] != '*' and candidate_mol.HasSubstructMatch(bin_information['scaffold_mol']):
                # Check if the molecule fits into the grouping bin
                for model in bin_information['models'].keys():
                    if model not in model_predictions.keys():
                        continue
                    model_bin_values = bin_information['models'][model]['grouping_detail']
                    lower_bound = model_bin_values['center'] - model_bin_values['width']
                    upper_bound = model_bin_values['center'] + model_bin_values['width']
                    if not lower_bound < model_predictions[model]['predicted_value'] < upper_bound:
                        break
                else:
                    bin_information['molecules'][candidate_key] = {}

    # Now evaluate the value of molecules in each of the grouping bins
    for grouping_bin in grouping_bins:
        current_bin = grouping_bins[grouping_bin]
        bin_models = current_bin['models']
        for bin_model in bin_models:
            predicted_values = []
            predicted_uncertainty = []
            for candidate_key in current_bin['molecules']:
                model_predictions = candidate_predictions[candidate_key]
                if bin_model not in model_predictions:
                    continue
                predicted_values.append(model_predictions[bin_model]['predicted_value'])
                predicted_uncertainty.append(model_predictions[bin_model]['prediction_uncertainty'])
            bin_models[bin_model]['stats'] = {'values': {'mean': np.mean(predicted_values),
                                                         'stdev': np.std(predicted_values),
                                                         'max_value': np.max(predicted_values),
                                                         'min_value': np.min(predicted_values)},
                                              'uncertainty': {'mean': np.mean(predicted_uncertainty),
                                                              'stdev': np.std(predicted_uncertainty),
                                                              'max_value': np.max(predicted_uncertainty),
                                                              'min_value': np.min(predicted_uncertainty)}}
            test_values = np.linspace(start=bin_models[bin_model]['stats']['uncertainty']['min_value'],
                                      stop=bin_models[bin_model]['stats']['uncertainty']['max_value'],
                                      num=1000)
            ens_values = lognorm.pdf(x=test_values, scale=bin_models[bin_model]['stats']['uncertainty']['mean'],
                                     s=bin_models[bin_model]['stats']['uncertainty']['stdev'])
            bin_models[bin_model]['stats']['uncertainty']['normalization'] = max(ens_values)
    molecule_values = {}
    for grouping_bin in grouping_bins:
        current_bin = grouping_bins[grouping_bin]
        bin_models = current_bin['models']
        for candidate_index, candidate_key in enumerate(current_bin['molecules']):
            property_values = []
            uncertainty_values = []
            model_predictions = candidate_predictions[candidate_key]
            if any(model not in model_predictions.keys() for model in constraining_models):
                continue
            for model_key in bin_models.keys():
                grouping_details = bin_models[model_key]['grouping_detail']
                property_value_stats = bin_models[model_key]['stats']['values']
                model_prediction = model_predictions[model_key]
                if grouping_details['type'] == 'min':
                    property_value = np.exp(-5 * (model_prediction['predicted_value'] - property_value_stats['min_value']) /
                                            (property_value_stats['max_value'] - property_value_stats['min_value']))
                elif grouping_details['type'] == 'max':
                    property_value = np.exp(-5 * (property_value_stats['max_value'] - model_prediction['predicted_value']) /
                                            (property_value_stats['max_value'] - property_value_stats['min_value']))
                elif grouping_details['type'] == 'target':
                    target_value = grouping_details['center']
                    target_width = grouping_details['width']
                    property_value = np.exp(-2 * ((2 * np.sqrt(2 * np.log(2)) * (model_prediction['predicted_value'] - target_value)) / target_width)^2)
                else:
                    return ['Error', 'Value of type %s is not implemented' % grouping_details['type']]
                if property_value > 1:
                    property_value = 1
                elif property_value < 0:
                    property_value = 0
                property_values.append(property_value)

                uncertainty_value_stats = bin_models[model_key]['stats']['uncertainty']
                if 'uncertainty' in grouping_information.keys() and grouping_information['uncertainty']['consider']:
                    if grouping_information['uncertainty']['method'] == 'exploration':
                        lognorm_value = lognorm.pdf(x=model_prediction['prediction_uncertainty'],
                                                    scale=uncertainty_value_stats['mean'],
                                                    s=2*uncertainty_value_stats['stdev'])
                        uncertainty_value = lognorm_value / uncertainty_value_stats['normalization']
                    elif grouping_information['uncertainty']['method'] == 'exploitation':
                        uncertainty_value = np.exp(-5 * (model_prediction['prediction_uncertainty'] - uncertainty_value_stats['min_value'])/
                                                   (uncertainty_value_stats['mean']+2*uncertainty_value_stats['stdev'] - uncertainty_value_stats['min_value']))
                    else:
                        return ['Error', 'Grouping method of type %s not implemented' % grouping_information['method']]
                    uncertainty_values.append(uncertainty_value)
                else:
                    uncertainty_values.append(1)
            if 'uncertainty' in grouping_information.keys() and grouping_information['uncertainty']['method'] == 'exploitation':
                scaling = 1
            else:
                scaling = 1
            current_bin['molecules'][candidate_key].update({'prediction_value': float(np.product(property_values))**scaling,
                                                            'uncertainty_value': float(np.product(uncertainty_values))})
            if candidate_key not in molecule_values.keys():
                molecule_values[candidate_key] = {'prediction_value': float(np.product(property_values))**scaling,
                                                  'uncertainty_value': float(np.product(uncertainty_values))}
            elif float(np.product(property_values))**scaling > molecule_values[candidate_key]['prediction_value']:
                molecule_values[candidate_key]['prediction_value'] = float(np.product(property_values)**scaling)

    return ['Success', grouping_bins, molecule_values]


def temperature_sequence_finding(job_information, candidate_reactions, molecule_values=None, sequence_type='standard'):
    if molecule_values is None:
        molecule_values = {}
    grouping_information = job_information['grouping_information']
    if sequence_type == 'prefilter':
        number_of_plates = grouping_information['number_of_plates'] + 2
    else:
        number_of_plates = grouping_information['number_of_plates']
    temperature_spacing = grouping_information['temperature_spacing']
    number_of_iterations = grouping_information['max_grouping_iterations']
    reaction_bin_values = list(range(*grouping_information['allowed_temperatures'], temperature_spacing))
    print('\tWorking to narrow down the reaction set to groups for %s wellplates' % number_of_plates)

    # Iteratively test temperature sequences
    # Want to select temperature sequences with the most value
    #   Higher ranked trees and contexts are better
    # Assign value based on molecule values from evaluation of predictions

    if molecule_values:
        number_of_products = len(list(molecule_values.keys()))
    else:
        number_of_products = len(list(candidate_reactions.keys()))
    temperature_sequence = []
    temperature_bins_values = list(range(*grouping_information['allowed_temperatures'], temperature_spacing))
    best_completed_reactions = set()
    best_completed_contexts = set()
    best_completed_products = set()
    best_products_value = 0
    all_products_sequence = 0
    for plate_sequence in range(0, number_of_plates):
        products_captured = []
        for plate_iteration in range(0, number_of_iterations):
            completed_reactions = set()
            completed_reaction_contexts = set()
            completed_products_set = set()
            completed_products = []
            test_temperature_sequence = []
            for test_plate in range(0, number_of_plates):
                temperature_bins = {}
                for temperature_value in temperature_bins_values:
                    temperature_bins[temperature_value] = []
                for target_product in candidate_reactions.keys():
                    if target_product in completed_products_set:
                        continue
                    for reaction_tree in candidate_reactions[target_product].keys():
                        reaction_keys = list(candidate_reactions[target_product][reaction_tree].keys())
                        for r_index, reaction in enumerate(reaction_keys):
                            reaction_key = reaction
                            if '.'.join([target_product, reaction_tree, reaction_key]) in completed_reactions:
                                continue
                            reaction_details = candidate_reactions[target_product][reaction_tree][reaction_key]
                            for context_index, reaction_context in enumerate(reaction_details['contexts'].keys()):
                                context_temperature = reaction_details['contexts'][reaction_context]['temperature']
                                nearest_temperature_bin = temperature_spacing * round(context_temperature / temperature_spacing)
                                try:
                                    if molecule_values and target_product in molecule_values.keys():
                                        product_value = molecule_values[target_product]['prediction_value'] * \
                                                        molecule_values[target_product]['uncertainty_value']
                                    else:
                                        product_value = 1
                                    temperature_bins[nearest_temperature_bin].append([target_product, reaction_tree,
                                                                                      reaction_key, reaction_context,
                                                                                      (r_index + 1) / len(reaction_keys)])
                                except KeyError:
                                    continue
                            break
                temp_bin_ranking = []
                for temp_bin in temperature_bins.keys():
                    if len(temperature_bins[temp_bin]) == 0:
                        continue
                    completion_value = sum([item[4] for item in temperature_bins[temp_bin]])
                    temp_bin_ranking.append([temp_bin, completion_value])
                temp_bin_ranking = sorted(temp_bin_ranking, key=lambda x: x[1], reverse=True)
                if len(temp_bin_ranking) == 0:
                    continue
                if test_plate < plate_sequence != 0:
                    max_temperature_bin = temperature_sequence[test_plate][0]
                else:
                    if plate_iteration >= len(temp_bin_ranking):
                        break
                    max_temperature_bin = temp_bin_ranking[plate_iteration][0]
                test_temperature_sequence.append(max_temperature_bin)
                products_value = 0
                for item_index, item in enumerate(temperature_bins[max_temperature_bin]):
                    completed_reaction_contexts.add('.'.join(item[0:4]))
                    if abs(item[4] - 1) < 0.00005:
                        completed_products_set.add(item[0])
                    if '.'.join(item[0:3]) not in completed_reactions:
                        completed_reactions.add('.'.join(item[0:3]))
                    if molecule_values and item[0] in molecule_values.keys():
                        product_value = molecule_values[item[0]]['prediction_value'] * molecule_values[item[0]]['uncertainty_value']
                        products_value += product_value**2
                    else:
                        products_value += 1
                if molecule_values:
                    if products_value > best_products_value:
                        best_completed_reactions = deepcopy(completed_reactions)
                        best_completed_contexts = deepcopy(completed_reaction_contexts)
                        best_completed_products = deepcopy(completed_products_set)
                else:
                    if len(completed_products_set) > len(best_completed_products):
                        best_completed_reactions = deepcopy(completed_reactions)
                        best_completed_contexts = deepcopy(completed_reaction_contexts)
                        best_completed_products = deepcopy(completed_products_set)
                completed_products.append(len(completed_products_set))
                if len(completed_products_set) == len(list(candidate_reactions.keys())):
                    # all_products_sequence = 1
                    continue
            else:
                if len(test_temperature_sequence) < plate_sequence + 1:
                    continue
                if len(completed_products) == 0:
                    print(test_temperature_sequence)
                    products_captured.append([test_temperature_sequence[plate_sequence], 0])
                else:
                    print(test_temperature_sequence)
                    print(completed_products)
                    print(plate_sequence)
                    products_captured.append([test_temperature_sequence[plate_sequence], max(completed_products)])
            if all_products_sequence == 1:
                continue
        if all_products_sequence == 1:
            print('Temperature sequence captures all of the available products')
            temperature_sequence = [[item, ''] for item in test_temperature_sequence]
            break
        sorted_products_captured = sorted(products_captured, key=lambda x: x[1], reverse=True)
        if len(sorted_products_captured) == 0:
            break
        temperature_sequence.append(sorted_products_captured[0])
        print('\t\tWorking on temperature sequence: %s' % temperature_sequence)

    final_temperature_sequence = [item[0] for item in temperature_sequence]

    products_to_remove = set()
    for product_key in candidate_reactions.keys():
        reaction_tree_to_remove = set()
        for reaction_tree_key in candidate_reactions[product_key].keys():
            reaction_tree = candidate_reactions[product_key][reaction_tree_key]
            number_of_reactions = len(list(reaction_tree.keys()))
            for reaction_index, reaction_key in enumerate(reaction_tree.keys()):
                contexts_to_remove = []
                for context in reaction_tree[reaction_key]['contexts'].keys():
                    context_details = reaction_tree[reaction_key]['contexts'][context]
                    temperature_key = temperature_spacing * round(context_details['temperature'] / temperature_spacing)
                    if temperature_key not in final_temperature_sequence:
                        contexts_to_remove.append(context)
                for context_to_remove in contexts_to_remove:
                    reaction_tree[reaction_key]['contexts'].pop(context_to_remove, None)
                if len(list(reaction_tree[reaction_key]['contexts'].keys())) == 0:
                    reaction_tree_to_remove.add(reaction_tree_key)
        for tree_to_remove in list(reaction_tree_to_remove):
            candidate_reactions[product_key].pop(tree_to_remove, None)
            if len(list(candidate_reactions[product_key].keys())) == 0:
                products_to_remove.add(product_key)
    for product_to_remove in products_to_remove:
        candidate_reactions.pop(product_to_remove, None)

    return ['Success', final_temperature_sequence, candidate_reactions]



