#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 14 10:35:14 2022

function to correct context

@author: mattmcdonald
"""

import json
from reaction_planning.api_client import APIClient
import copy
import os
from pprint import pprint
import yaml
import traceback
import sys
import time

run_from_spyder = True

if not run_from_spyder:
    from rdkit.Chem import AllChem as ac

# Put your file routes here
rxn_rules_file = r'.\reaction_planning\context_correction\rules_file.json'
substitutions_file = r'.\reaction_planning\context_correction\ASKCOS_substitution_list.json'
name_subs_file = r'.\reaction_planning\context_correction\ASKCOS_name_only_subs.json'
#disallowed_file = r'.json'
#replacement_chems_file = r'.json'

with open(rxn_rules_file, 'r') as json_in:
    rxn_rules = json.load(json_in)    
with open(substitutions_file, 'r') as json_in:
    subs_dict = json.load(json_in)        
with open(name_subs_file, 'r') as json_in:
    name_subs_dict = json.load(json_in)
#with open(disallowed_file, 'r') as json_in:
#    banned_reagents = json.load(json_in)
#with open(replacement_chems_file, 'r') as json_in:
#    chem_subs = json.load(json_in)

rxn_rules_status = {}

# rxn_rules = see: https://docs.google.com/document/d/1VQfULLNfWxZGKKstPMIFiIMF8ZwqK5stl8O7E8R5XDg/edit
# rxn_class = see above
# banned_reagents = dictionary of banned reagents and substitutes
# chem_subs = dictionary of preferred chemicals to be substituted for uncommon/expensive/etc. chemicals
# subs_dict = dictionary of smiles ASKCOS parses incorrectly and their correct smiles for substitution
# name_subs_dict = dictionary of ASKCOS name_only chemicals and their correct smiles


def contexts_handler(rxn_smiles: str, ASKCOS_condition_recommendation_return: dict):
    """
    given a reaction smiles and ASKCOS context recommendations return re-ranked, corrected contexts

    Parameters
    ----------
    contexts : list
        raw return from ASKCOS context recommender.
    rxn_smiles : str
        reaction smiles that generated context.

    Returns
    -------
    an executable context dictionary

    """
    contexts = ASKCOS_condition_recommendation_return['output']
    cleaned_contexts = preclean_contexts(contexts)
    
    rxn_class = predict_rxn_class(rxn_smiles)
    # If there are no rules enumerated for that class then return the re-ranked contexts
    if rxn_class not in rxn_rules.keys():
        print("Unable to classify reaction. Contexts re-ranked by forward predictor only (no corrections applied)")
        ASKCOS_condition_recommendation_return['output'] = forward_prediction_rerank(cleaned_contexts, rxn_smiles)
        return ASKCOS_condition_recommendation_return
    
    print(f"Working on a {rxn_class} reaction: \n\t{rxn_smiles}")
    valid_contexts = []
    rxn_rules_status[rxn_class] = []
    for context in cleaned_contexts:
        # call to clean a raw ASKCOS context based on logic in precleaning.py, also adds in default concs
        check, augmented_context = check_context(context, rxn_class, rxn_smiles)
        if check:
            valid_contexts.append(augmented_context)
    
    if not valid_contexts:
        valid_contexts = cleaned_contexts
        print(f"no valid contexts found matching rules for a {rxn_class}. Rerank based on forward prediction only")
        #valid_contexts.append(rxn_rules[rxn_class]['default_context'])
                                              
    ASKCOS_condition_recommendation_return['output'] = forward_prediction_rerank(valid_contexts, rxn_smiles)
    return ASKCOS_condition_recommendation_return
    
def predict_rxn_class(rxn_smiles: str):
    """
    predict which class the reaction belongs to with ASKCOS reaction-classification
    Assess the probability to determine if its an accurate prediction

    Parameters
    ----------
    rxn_smiles : str
        reaction smiles

    Returns
    -------
    None.

    """
    hostname = r'https://askcos.mit.edu:7000/'
    client = APIClient(hostname, verify=False)
    
    reactants = rxn_smiles.split('>>')[0]
    products = rxn_smiles.split('>>')[1]
    
    params = {'reactants': reactants,
              'products': products}
    result = client.post('reaction-classification', data = params)
    try:
        rxn_class_result = result['task_id']
    except KeyError:
        print(f'error getting task_id on {rxn_smiles}')
        # this seems to be caused by reagents that can't be parsed by rdkit
    
    try:
        task_result = client.get_result(rxn_class_result, timeout=600, interval=5)
    except json.JSONDecodeError:
        print('\tJSON DECODE ERROR!!! Problem with it...')
    
    try:
        rxn_class = task_result['output']['result'][0]['reaction_name']
        probability = task_result['output']['result'][0]['prediction_certainty']
    except KeyError:
        print(f'failed to classify reaction: {rxn_smiles}')
    except TypeError:
        probability = 0
    cutoff = 0.5
    if probability > cutoff:
        return rxn_class
    else:
        return 'No confident reaction class predicted'
    
def forward_prediction_rerank(contexts: list, rxn_smiles: str):
    """
    rank a list of reaction contexts based on probability of producing the desired product 
    estimated by the forward reaction predictor in askcos

    Parameters
    ----------
    contexts : list
        DESCRIPTION.
    rxn_smiles : str
        DESCRIPTION.

    Returns
    -------
    a list of reaction contexts in rank order

    """
    product_smiles = rxn_smiles.split('>>')[1]
    reactant_smiles = rxn_smiles.split('>>')[0]
    reactions_to_request = []
    for i, context in enumerate(contexts):
        if context['reagent'] or context['catalyst']:
            reagents = ""
            for reagent in context['reagent']:                
                reagents = reagents + reagent + "."
            for catalyst in context['catalyst']:
                reagents = reagents + catalyst + "."
            reagents = reagents[:-1]
        else:
            reagents = ""
        rxn_dict = {'reaction': str(i) + "_" + rxn_smiles,
                    'reactants': reactant_smiles,
                    'products': product_smiles,
                    'reagents': reagents,
                    'solvent': context['solvent'],
                    'temperature': context['temperature'],
                    'score': context['score']}
        reactions_to_request.append(rxn_dict)
    
    hostname = r'https://askcos.mit.edu:7000/'
    client = APIClient(hostname, verify=False)
    size_of_groups = 10
    request_groups = [reactions_to_request[i*size_of_groups:(i+1)*size_of_groups] for
                      i in range((len(reactions_to_request)+size_of_groups-1)//size_of_groups)]
    
    forwards = {}
    for idx, group in enumerate(request_groups):
        current_ids = []
        for rxn in group:
            if rxn['reaction'] not in forwards.keys():
                params = {'reactants':rxn['reactants'],
                          'reagents':rxn['reagents'],
                          'solvent':rxn['solvent'],
                          'num_results': 10}
                result = client.post('forward', data = params)
                try:
                    current_ids.append([rxn['reaction'], result['task_id']])
                except KeyError:
                    print(f'\terror getting task_id on {rxn_smiles}:')
                    print(f'\t\tLikely unidentifiable chemical, check for _XXX: {params["reagents"]}')
                    # this seems to be caused by reagents that can't be parsed by rdkit
        try:
            task_results = [client.get_result(task_id[1], timeout=600, interval=1) for task_id in current_ids]
        except json.JSONDecodeError:
            print('\tJSON decode error')
            continue

        for current_index, task_result in enumerate(task_results):
            forwards[current_ids[current_index][0]] = task_result
            if task_result is None:
                continue
            products = task_result['output']
            for product in products:
                if not isinstance(product, dict):
                    continue
                if product['smiles'] == product_smiles:
                    contexts[int(current_ids[current_index][0].split('_')[0])]['probability'] = product['prob']
                    break
                contexts[int(current_ids[current_index][0].split('_')[0])]['probability'] = 0.0
            #print('\tCompleted %s, returned %s predictions' % (current_ids[current_index][0], len(task_result['output'])))   
    try:
        return_contexts = sorted(contexts, key=lambda i: i['probability'], reverse=True)  
    except:
        #pprint(contexts)
        return_contexts = contexts
        
    return return_contexts

def check_context(context: dict, rxn_class: str, rxn_smiles: str):
    """
    check_context checks that a given reaction context matches the reaction
    template used to generate the retrosynthetic reaction.
    How this is executed is up to you...
 
    Parameters
    ----------
    context : dict
        a dictionary with fields:
            reagents : list (a list of smiles strings)
            catalysts : list (a list of smiles strings)
            solvents : list (a list of smiles strings)
            temperature : int (the recommended temperature)
    rxn_class : str
        a string corresponding to an ASKCOS reaction classification.
    rxn_smiles : str
        a string representing the reaction smiles
 
    Returns
    -------
    complete : bool
        True if it is a valid context
    augmented_context : dict
        a dictionary with all the info for an augmented context, see:
        https://docs.google.com/document/d/1VQfULLNfWxZGKKstPMIFiIMF8ZwqK5stl8O7E8R5XDg/edit
 
    """
    try:
        rxn_sets = rxn_rules[rxn_class]['reaction_sets']
    except KeyError:
        rxn_sets = [rxn_rules[rxn_class]]
        #print("warning: no sets defined, assuming only one valid reaction set")
    complete = False
        
    for rxn_set in rxn_sets:
        
        reagents_complete = True
        missing_reagents = []
        for req_reagent_set in rxn_set['reagents'].values():
            if req_reagent_set:
                in_context_r = bool(set(req_reagent_set.split(',')) & set(context['reagent']))
                in_context_c = bool(set(req_reagent_set.split(',')) & set(context['catalyst']))
                in_context = in_context_r or in_context_c
                #^^this should be changed when we have cats_complete working correctly
                reagents_complete = reagents_complete and in_context
                if not in_context:
                    missing_reagents.append(req_reagent_set.split(',')[0])
        
        solvs_complete = bool(set(rxn_set['solvent'].split(',')) & set(context['solvent']))
        if not solvs_complete:
            missing_solv = rxn_set['solvent'].split(',')[0]
            
        # crude implementation of catalysts as being separate from reagents, check that the element is there
        if 'catalyst' not in rxn_set.keys():
            cats_complete = True
        else:
            try:
                cat_metal = rxn_set['catalyst'].split(',')[0]
                cats_complete = (cat_metal in ''.join(context['catalyst'])) or (cat_metal in ''.join(context['reagent']))
            except KeyError:
                #print(f"unable to check catalyst: {context['catalyst']}")
                pass

        if not cats_complete:
            missing_cat = rxn_set['catalyst'].split(',')[1]
        
        temp_complete = True
        if (context['temperature'] >= int(rxn_set['T'].split(',')[0])) and \
        (context['temperature'] <= int(rxn_set['T'].split(',')[1])):
            temp_complete = True
        else:
            temp_complete = False
            
        rxn_rules_status[rxn_class].append([reagents_complete, cats_complete, 
                                            solvs_complete, temp_complete]) 
            
        augmented_reagents = context['reagent']
        augmented_cats = context['catalyst']
        augmented_solvs = context['solvent']
        augmented_temp = context['temperature']
        
        # the case where only one of the categories is incomplete... fix it
        sum_of_conditions = reagents_complete + cats_complete + solvs_complete + temp_complete
        if sum_of_conditions == 3:
            if not reagents_complete:
                if len(missing_reagents) > 1:
                    # print(missing_reagents)
                    continue
                augmented_reagents.append(missing_reagents[0])
            if not solvs_complete:
                augmented_solvs.append(missing_solv)
            if not cats_complete:
                augmented_cats.append(missing_cat)
            if not temp_complete:
                augmented_temp = int(int(rxn_set['T'].split(',')[0]) + int(rxn_set['T'].split(',')[1])/2)
            complete = True
            break
        elif sum_of_conditions == 4:
            complete = True
            break
    
    # some code to check that the reactants match the expected substructure for the reaction class
    # not yet implemented in the reaction rules dictionary...
    augmented_context = {}
    """
    if complete:
        reactants = rxn_smiles.split('>>')[0].split(".")
        augmented_reactants = []
        for reactant in reactants:
            if not run_from_spyder:
                mol = ac.MolFromSmiles(reactants)
                for req_reactant in rxn_set['reactants'].split(','):
                    p = ac.MolFromSmarts(req_reactant)
                    if mol.HasSubstructMatch(p):
                        augmented_reactants.append(reactant)
                        break
            else:
                augmented_reactants.append(reactant)
    """
    if complete:        
        augmented_context = {'reaction_smiles': rxn_smiles,
                             'reagent': augmented_reagents,
                             'catalyst': augmented_cats,
                             'solvent': augmented_solvs,
                             'temperature': augmented_temp,
                             'ordered': 'ordered' in rxn_set['notes'],
                             'air_free': rxn_set['air free'],
                             'score': context['score']}
    
    return complete, augmented_context

def preclean_contexts(contexts_raw: list):
    """
    preclean_contexts takes contexts as received from ASKCOS and returns a cleaned up dictionary

    Parameters
    ----------
    context_files : list
        DESCRIPTION.

    Returns
    -------
    contexts_cleaned : TYPE
        DESCRIPTION.

    """
    error_parsing = []
    # need to add section for banned chemicals and allowed substitutions
    # put reagents, catalysts, and solvents into a list
    name_only = ['catalyst_name_only', 'reagent_name_only', 'solvent_name_only']
    smiles_only = ['catalyst', 'reagent', 'solvent']
    
    for context in contexts_raw:
        i = 0
        for category in name_only:
            if context[category]:
                try:
                    if context[smiles_only[i]]:
                        context[smiles_only[i]] = context[smiles_only[i]] + '.' + name_subs_dict[context[category]]
                    else:
                        context[smiles_only[i]] = name_subs_dict[context[category]]
                except (ValueError, KeyError) as e:
                    with open(r'.\context_problems.csv', 'a') as csvfile:
                        csvfile.write(str(traceback.format_exc()) + '\n')
                    pass
            i += 1
    
    cleaned_contexts = []
    for context in contexts_raw:
        try:
            cleaned_context = {'temperature': round(context['temperature']),
                               'score': round(context['score'],6)}
        except KeyError:
            cleaned_context = {'temperature': round(context['temperature']),
                               'score': 'not calculated'}
        for category in smiles_only:
            chems = context[category]
            # substitute out bad SMILES
            for sub in subs_dict.keys():
                try:
                    start_sub = chems.index(sub)
                    chems = chems[0:start_sub] + subs_dict[sub] + chems[start_sub+len(sub):len(chems)]
                except ValueError:
                    pass
            # group salts and add chemicals as individual items in list    
            index = 0
            charge = 0
            cation = False
            anion = False
            in_salt = False
            chem_start = 0
            reagents = []
            for char in chems:
                if cation:
                    try:
                        charge += int(char)
                    except ValueError:
                        charge +=1
                    cation = False
                if anion:
                    try:
                        charge -= int(char)
                    except ValueError:
                        charge -=1
                    anion = False
                    
                if char == "+":
                    cation = True
                    in_salt = True
                if char == "-":
                    anion = True
                    in_salt = True
                
                if in_salt and (charge == 0):
                    if char == ".":
                        reagents.append(chems[chem_start:index])
                        chem_start = index+1
                        in_salt = False
                
                if (char == ".") and not in_salt:
                    reagents.append(chems[chem_start:index])
                    chem_start = index+1
                    
                index +=1
            if in_salt and (charge == 0):
                reagents.append(chems[chem_start:index])
            elif charge == 0:
                reagents.append(chems[chem_start:index])
            else:
                if chems not in error_parsing:
                    print(f"error parsing: {chems}")
                    error_parsing.append(chems)
            try:
                reagents.remove('')
            except:
                pass
            # rdkit canonicalization
            if not run_from_spyder:
                rdkit_reagents = []
                error_rdkit = []
                for reagent_chem in reagents:
                    try:
                        mol = ac.MolFromSmiles(reagent_chem)
                        rdkit_reagents.append(ac.MolToSmiles(mol))
                    except:
                        error_rdkit.append(reagent_chem)
            
            cleaned_context[category] = reagents
        cleaned_contexts.append(cleaned_context)
        
    return cleaned_contexts

if __name__ == "__main__":
    test_rxn = 'C=Cc1ccccc1.Nc1nc2c(cc(I)c3[nH]c4ncccc4c(=O)c32)s1>>Nc1nc2c(cc(C=Cc3ccccc3)c3[nH]c4ncccc4c(=O)c32)s1'
    test_rxn = 'C=[N+]=[N-].O=c1c2c(F)ccnc2[nH]c2ccc3sc(-c4ccc(O)cc4)nc3c12>>COc1ccc(-c2nc3c(ccc4[nH]c5nccc(F)c5c(=O)c43)s2)cc1'
    test_rxn = 'C1CCNCC1.Fc1ccc(Br)cc1>>Fc1ccc(N2CCCCC2)cc1'
    rxn_dict = {'reactants': test_rxn.split('>>')[0],
                'products': test_rxn.split('>>')[1],
                'num_results': 10}

    hostname = r'https://askcos.mit.edu:7000/'
    client = APIClient(hostname, verify=False)
    current_ids = []
    result = client.post('context', data = rxn_dict)
    try:
        current_ids.append([test_rxn, result['task_id']])
    except KeyError:
        print(f'\terror getting task_id on {test_rxn}')
        # this seems to be caused by reagents that can't be parsed by rdkit
    try:
        task_results = [client.get_result(tid[1], timeout=600, interval=5) for tid in current_ids]
    except json.JSONDecodeError:
        print('\tJSON DECODE ERROR!!! Problem with it...')
    original = copy.deepcopy(task_results[0])
    reranked = contexts_handler(test_rxn, task_results[0])


