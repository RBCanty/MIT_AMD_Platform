
import json
import os
import yaml
import glob
import sys

from copy import deepcopy
from pprint import pprint
from rdkit import Chem
from rdkit.Chem import AllChem

from chemical_library.chemical_inventory_functions import load_reagents_datasets
from functions.supporting_functions import print_progress_bar


def askcos_template_expansion(job_information, job_name, reaction_planning_folder):
    # Start by preparing an expansion status file
    new_reaction_plans = False
    expansion_directory = r'%s\pathway_expansion' % reaction_planning_folder
    expansion_status_filepath = r'%s\%s_expansion_status.yaml' % (expansion_directory, job_name)
    if os.path.exists(expansion_status_filepath):
        print('\t\tOpening the expansion status file')
        with open(expansion_status_filepath, 'r') as yaml_file:
            expansion_status = yaml.load(yaml_file, Loader=yaml.Loader)
    else:
        print('\t\tNeed to make a new expansion status file')
        expansion_status = {'products': [],
                            'completed_products': [],
                            'recompile_expansion_files': True}

    # Need to figure out which reaction pathways can be expanded
    reaction_details_path = r'%s\%s_compiled_askcos.json' % (reaction_planning_folder, job_name)
    with open(reaction_details_path, 'r') as json_file:
        reaction_details = json.load(json_file)
    if not reaction_details:
        return ['Error', 'No reaction details present for askcos template expansion', new_reaction_plans]
    products = list(reaction_details.keys())
    expansion_status['products'] = list(set(expansion_status['products'] + products))
    with open(expansion_status_filepath, 'w') as yaml_file:
        yaml.dump(expansion_status, yaml_file)

    # A small exit to avoid loading anything if no products to expand
    if set(expansion_status['products']) == set(expansion_status['completed_products']):
        if expansion_status['recompile_expansion_files']:
            compiled_output_path = r'%s\%s_compiled_products.json' % (reaction_planning_folder, job_name)
            return_statement, return_object = gather_expanded_products(expansion_directory)
            if return_statement != 'Success':
                return ['Error', return_object]
            with open(compiled_output_path, 'w') as json_file:
                json.dump(return_object, json_file)
            expansion_status['recompile_expansion_files'] = False
            with open(expansion_status_filepath, 'w') as yaml_file:
                yaml.dump(expansion_status, yaml_file)
            new_reaction_plans = True
        return ['Success', 'No additional products need template expansion', new_reaction_plans]

    # With the bookkeeping finished we can load reagents and expand reaction trees
    reagent_details = job_information['reagent_information']
    expansion_price_limit = job_information['pathway_expansion']['expansion_price_cutoff']
    return_statement, reagent_library, chemistry_inventory = load_reagents_datasets(reagent_details,
                                                                                    price_limit=expansion_price_limit)
    if return_statement != 'Success':
        return ['Error', reagent_library, new_reaction_plans]

    # Load the ASKCOS template set
    templates_filepath = r'.\reaction_planning\template_expansion\askcos_template_set.json'
    with open(templates_filepath, 'r') as json_file:
        askcos_template_set = json.load(json_file)

    # We will start expanding the reaction_trees
    products_to_expand = set(expansion_status['products']).difference(set(expansion_status['completed_products']))
    print('\tBeginning template expansion for %s products' % len(products_to_expand))
    save_dataset = {}
    current_products_to_save = []
    for product_index, product_smiles in enumerate(products_to_expand):
        current_products_to_save.append(product_smiles)
        # There is a leftover problem from molecular generation and salts
        if '.' in product_smiles:
            continue
        if product_smiles not in reaction_details.keys():
            continue
        print('\t\tWorking on (%s/%s): %s' % (product_index + 1, len(products_to_expand), product_smiles))
        product_details = reaction_details[product_smiles]
        for tree_index, reaction_tree_key in enumerate(reaction_details[product_smiles].keys()):
            print('\t\t\tWorking on (%s/%s)' % (tree_index + 1, len(reaction_details[product_smiles].keys())))
            reaction_tree = product_details[reaction_tree_key]
            return_statement, expanded_product = expand_reaction_tree(reaction_tree, reagent_library,
                                                                      askcos_template_set, job_information)
            if return_statement != 'Success':
                return ['Error', expanded_product, new_reaction_plans]

            for new_product_smiles in expanded_product.keys():
                if new_product_smiles not in save_dataset.keys():
                    save_dataset[new_product_smiles] = {}
                for reaction_tree in expanded_product[new_product_smiles].keys():
                    tree_index = 1
                    while True:
                        tree_key = '%s_tree_%s' % (new_product_smiles, tree_index)
                        if tree_key not in save_dataset[new_product_smiles].keys():
                            break
                        tree_index += 1
                        if tree_index > 50:
                            break
                    save_dataset[new_product_smiles][tree_key] = expanded_product[new_product_smiles][reaction_tree]

        # Now we can think about saving the expanded trees
        if (product_index + 1) % 5 == 0:
            file_index = 0
            while os.path.exists(r'%s\%s_expanded_products_%s.json' % (expansion_directory, job_name, file_index)):
                file_index += 1
            file_outpath = r'%s\%s_expanded_products_%s.json' % (expansion_directory, job_name, file_index)
            with open(file_outpath, 'w') as json_file:
                json.dump(save_dataset, json_file)
            with open(expansion_status_filepath, 'r') as yaml_file:
                expansion_status = yaml.load(yaml_file, Loader=yaml.Loader)
            expansion_status['completed_products'].extend(current_products_to_save)
            expansion_status['recompile_expansion_files'] = True
            with open(expansion_status_filepath, 'w') as yaml_file:
                yaml.dump(expansion_status, yaml_file)
            save_dataset = {}
            current_products_to_save = []

    # Upon leaving the expansion loop, save any leftover trees
    if len(save_dataset.keys()) != 0:
        file_index = 0
        while os.path.exists(r'%s\%s_expanded_products_%s.json' % (expansion_directory, job_name, file_index)):
            file_index += 1
        file_outpath = r'%s\%s_expanded_products_%s.json' % (expansion_directory, job_name, file_index)
        with open(file_outpath, 'w') as json_file:
            json.dump(save_dataset, json_file)
        with open(expansion_status_filepath, 'r') as yaml_file:
            expansion_status = yaml.load(yaml_file, Loader=yaml.Loader)
        expansion_status['completed_products'].extend(current_products_to_save)
        expansion_status['recompile_expansion_files'] = True
        with open(expansion_status_filepath, 'w') as yaml_file:
            yaml.dump(expansion_status, yaml_file)

    if expansion_status['recompile_expansion_files']:
        compiled_output_path = r'%s\%s_compiled_products.json' % (reaction_planning_folder, job_name)
        return_statement, return_object = gather_expanded_products(expansion_directory)
        if return_statement != 'Success':
            return ['Error', return_object, new_reaction_plans]
        with open(compiled_output_path, 'w') as json_file:
            json.dump(return_object, json_file)
        expansion_status['recompile_expansion_files'] = False
        new_reaction_plans = True
        with open(expansion_status_filepath, 'w') as yaml_file:
            yaml.dump(expansion_status, yaml_file)

    return ['Success', 'Template expanded for %s original products' % len(products_to_expand), new_reaction_plans]


def expand_reaction_tree(reaction_tree, reagent_library, askcos_template_set, job_information):
    # Start by finding the reaction sequence that we can expand due to dependencies
    expansion_details = job_information['pathway_expansion']
    reaction_sequence = []
    for iteration in range(0, 5):
        for reaction_key in reaction_tree.keys():
            if reaction_key in reaction_sequence:
                continue
            dependencies = reaction_tree[reaction_key]['dependencies']
            if set(dependencies).issubset(set(reaction_sequence)):
                reaction_sequence.append(reaction_key)
        if len(reaction_sequence) == len(list(reaction_tree.keys())):
            break
    else:
        return ['Error', 'Unable to determine a reaction sequence from dependencies']

    # Start working through the existing reactions
    original_reaction_tree = deepcopy(reaction_tree)
    for reaction_index, reaction_key in enumerate(reaction_sequence):
        # Find the predicted templates
        predicted_templates = reaction_tree[reaction_key]['templates']
        template_splits = []
        for predicted_template in predicted_templates:
            if '>>' in predicted_template:
                template_split = predicted_template.split('>>')
                template_splits.append([item.split('.') for item in template_split])
            else:
                template_split = askcos_template_set[predicted_template]['reaction_smarts'].split('>>')
                template_splits.append([item.split('.') for item in template_split])

        reagent_map = {}
        for reactant in reaction_tree[reaction_key]['reactants']:
            reactant_mol = Chem.MolFromSmiles(reactant)
            relevant_templates = []
            for template_split in template_splits:
                for template_reactant in template_split[1]:
                    template_reactant_pattern = Chem.MolFromSmarts(template_reactant)
                    if reactant_mol.HasSubstructMatch(template_reactant_pattern):
                        relevant_templates.append(template_reactant)
            reagent_map[reactant] = {'original_smiles': [reactant, reaction_tree[reaction_key]['reaction']],
                                     'derivative_smiles': [[reactant, reaction_tree[reaction_key]['reaction']]],
                                     'previous_product': False, 'templates': relevant_templates}

        # We need to find products from previous reactions that match a reagent template
        dependent_reactions = reaction_tree[reaction_key]['dependencies']
        for dependent_reaction in dependent_reactions:
            dependent_reaction_details = reaction_tree[dependent_reaction]
            if 'other_products' in dependent_reaction_details.keys():
                other_products = dependent_reaction_details['other_products']
            else:
                other_products = []
            for previous_product_reaction in other_products:
                previous_product = previous_product_reaction.split('>>')[1]
                previous_product_mol = Chem.MolFromSmiles(previous_product)
                for reagent in reagent_map:
                    template_mols = [Chem.MolFromSmarts(template) for template in reagent_map[reagent]['templates']]
                    if all([previous_product_mol.HasSubstructMatch(template) for template in template_mols]):
                        reagent_map[reagent]['derivative_smiles'].append([previous_product, previous_product_reaction])
                        reagent_map[reagent]['previous_product'] = True

        # Then we can fill in reagents into the other reagents
        for reagent in reagent_map:
            if reagent_map[reagent]['previous_product']:
                continue
            original_reagent_mol = Chem.MolFromSmiles(reagent_map[reagent]['original_smiles'][0])
            template_mols = [Chem.MolFromSmarts(template) for template in reagent_map[reagent]['templates']]
            for reagent_key in reagent_library.keys():
                chemical_mol = reagent_library[reagent_key]['mol']
                # Check to make sure that the new chemical has the original molecule as a substructure
                if not chemical_mol.HasSubstructMatch(original_reagent_mol):
                    continue
                if chemical_mol.GetNumAtoms() > original_reagent_mol.GetNumAtoms() + 6:
                    continue
                if all(chemical_mol.HasSubstructMatch(template) for template in template_mols):
                    reagent_map[reagent]['derivative_smiles'].append([reagent_key, ''])
                if len(reagent_map[reagent]['derivative_smiles']) > expansion_details['new_reagents_to_consider']:
                    break

        # Get combinations of reagents to use
        reagent_keys = list(reagent_map.keys())
        if len(reagent_keys) == 1:
            combinations = [[x] for x in reagent_map[reagent_keys[0]]['derivative_smiles']]
        elif len(reagent_keys) == 2:
            combinations = [[x, y] for x in reagent_map[reagent_keys[0]]['derivative_smiles']
                            for y in reagent_map[reagent_keys[1]]['derivative_smiles']]
        elif len(reagent_keys) == 3:
            combinations = [[x, y, z] for x in reagent_map[reagent_keys[0]]['derivative_smiles']
                            for y in reagent_map[reagent_keys[1]]['derivative_smiles']
                            for z in reagent_map[reagent_keys[2]]['derivative_smiles']]
        else:
            print('Broken.... %s' % reaction_tree[reaction_key]['reaction'])
            return ['', '']
        predicted_template = template_splits[0]
        rxn = AllChem.ReactionFromSmarts('.'.join(predicted_template[1]) + '>>' + '.'.join(predicted_template[0]))
        unique_reactions = set()
        for combination_index, combination in enumerate(combinations):
            if len(combination) != len(predicted_template[1]):
                continue
            print_progress_bar(combination_index + 1, len(combinations), prefix='\t\t\t\tWorking: ', length=50)
            unique_results = set()
            reagent_combination = []
            if combination_index > 1000:
                continue
            for item in combination:
                reagent_combination.append(item[0])
            results = rxn.RunReactants([Chem.MolFromSmiles(reagent) for reagent in reagent_combination])
            for item in results:
                unique_results.add(Chem.MolToSmiles(item[0]))
            if len(reagent_combination) != 1:
                results = rxn.RunReactants([Chem.MolFromSmiles(reagent) for reagent in reagent_combination[::-1]])
                for item in results:
                    unique_results.add(Chem.MolToSmiles(item[0]))
            for unique_result in list(unique_results):
                if 'other_products' not in reaction_tree[reaction_key].keys():
                    reaction_tree[reaction_key]['other_products'] = []
                new_reaction_smiles = '.'.join(reagent_combination) + '>>' + unique_result
                if new_reaction_smiles not in unique_reactions:
                    reaction_tree[reaction_key]['other_products'].append(new_reaction_smiles)
                    unique_reactions.add(new_reaction_smiles)

    # With the expanded reaction trees we can now rebuild the reaction trees
    returned_reaction_trees = {}
    for reaction_key in reaction_tree.keys():
        if 'other_products' not in reaction_tree[reaction_key].keys():
            reaction_tree[reaction_key]['other_products'] = [reaction_tree[reaction_key]['reaction']]
    for other_product in reaction_tree[reaction_sequence[-1]]['other_products']:
        product_smiles = other_product.split('>>')[1]
        tree_chemicals = set(other_product.split('>>')[0].split('.'))
        tree_chemicals.add(product_smiles)
        rebuilt_reaction_tree = [[reaction_sequence[-1], other_product]]
        for previous_reaction_key in reaction_sequence[::-1][1::]:
            for other_reaction in reaction_tree[previous_reaction_key]['other_products']:
                overlap = tree_chemicals.intersection(other_reaction.split('>>')[1].split('.'))
                if len(overlap) > 0:
                    tree_chemicals.update(other_reaction.split('>>')[0].split('.'))
                    rebuilt_reaction_tree.append([previous_reaction_key, other_reaction])
        if len(rebuilt_reaction_tree) != len(reaction_sequence):
            continue
        if product_smiles not in returned_reaction_trees.keys():
            returned_reaction_trees[product_smiles] = {}
        reaction_tree_key = '%s_tree_%s' % (product_smiles, len(returned_reaction_trees[product_smiles].keys()) + 1)
        if reaction_tree_key not in returned_reaction_trees[product_smiles].keys():
            returned_reaction_trees[product_smiles][reaction_tree_key] = {} # original_reaction_tree
        for rebuilt_reaction in rebuilt_reaction_tree[::-1]:
            returned_reaction_trees[product_smiles][reaction_tree_key][rebuilt_reaction[0]] = {}
            current_reaction = returned_reaction_trees[product_smiles][reaction_tree_key][rebuilt_reaction[0]]
            current_reaction['reaction'] = deepcopy(rebuilt_reaction[1])
            current_reaction['products'] = rebuilt_reaction[1].split('>>')[1].split('.')
            current_reaction['reactants'] = rebuilt_reaction[1].split('>>')[0].split('.')
            for key in original_reaction_tree[rebuilt_reaction[0]]:
                if key in current_reaction:
                    continue
                else:
                    current_reaction[key] = deepcopy(original_reaction_tree[rebuilt_reaction[0]][key])

        #if reaction_tree_key not in returned_reaction_trees[product_smiles].keys():
        #    returned_reaction_trees[product_smiles][reaction_tree_key] = original_reaction_tree
        #    for rebuilt_reaction in rebuilt_reaction_tree:
        #        print(rebuilt_reaction)
        #        returned_reaction_trees[product_smiles][reaction_tree_key][rebuilt_reaction[0]]['reaction'] = deepcopy(rebuilt_reaction[1])
        #        current_reaction = returned_reaction_trees[product_smiles][reaction_tree_key][rebuilt_reaction[0]]
        #        current_reaction['products'] = rebuilt_reaction[1].split('>>')[1].split('.')
        #        current_reaction['reactants'] = rebuilt_reaction[1].split('>>')[0].split('.')
        #        current_reaction['reaction'] = rebuilt_reaction[1]
        #        # print(current_reaction)
        #        
        #        #current_reaction.update({'products': rebuilt_reaction[1].split('>>')[1].split('.'),
        #        #                         'reactants': rebuilt_reaction[1].split('>>')[0].split('.'),
        #        #                         'reaction': rebuilt_reaction[1]})
        #        #returned_reaction_trees[product_smiles][reaction_tree_key][rebuilt_reaction[0]] = deepcopy(current_reaction)
    print('\t\t\t\tExpansion yielded %s unique products' % len(returned_reaction_trees.keys()))
    return ['Success', returned_reaction_trees]


def gather_expanded_products(expanded_products_output_path):
    expanded_product_dataset = {}
    expanded_files = glob.glob(r'%s\*_expanded_products_*.json' % expanded_products_output_path)
    if len(expanded_files) == 0:
        return ['Error', 'No expansion files found...']
    for file_index, expanded_file in enumerate(expanded_files):
        print_progress_bar(file_index + 1, len(expanded_files), prefix='\t\tWorking: ', length=50)
        with open(expanded_file, 'r') as jsonfile:
            incoming_trees = json.load(jsonfile)
        for product_key in incoming_trees.keys():
            if product_key not in expanded_product_dataset.keys():
                expanded_product_dataset[product_key] = {}
            for tree_key in incoming_trees[product_key].keys():
                tree_index = 1
                while True:
                    tree_name = product_key + '_Tree_%s' % tree_index
                    if tree_name not in expanded_product_dataset[product_key].keys():
                        break
                    tree_index += 1
                    if tree_index > 100:
                        break
                if tree_index > 100:
                    continue
                expanded_product_dataset[product_key][tree_name] = incoming_trees[product_key][tree_key]
    return ['Success', expanded_product_dataset]
