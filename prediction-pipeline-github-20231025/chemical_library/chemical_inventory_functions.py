# -*- coding: utf-8 -*-
"""
Created on 01/21/22

Condensing of the chemical inventory related functions into this file
 -- load_reagents_datasets
 
@author: brent
"""

from rdkit import Chem
import traceback
import sys
import yaml


def load_reagents_datasets(reagent_details, price_limit=100):
    # Establish the price bounds that will be considered
    print('\t\tWorking on building the reagent library')
    if 'price_cutoff' not in reagent_details.keys():
        reagent_details['price_cutoff'] = price_limit
    price_cutoff = min([price_limit, reagent_details['price_cutoff']])

    # Now work on building the reagents_dataset
    reagents_dataset = {}
    if reagent_details['reagent_origin'] == 'restricted':
        chemical_sources = reagent_details['reagent_sources']
        if 'ambeed' in chemical_sources:
            print('\t\tAdding the Ambeed Chemical Catalog to the reagents working set')
            with open(r'.\chemical_library\reagent_sources\ambeed_chemical_list.txt', 'r') as infile:
                for line in infile:
                    line_split = line.strip('\n').split(',')
                    if line_split[0] == '' or '@' in line_split[0]:
                        continue
                    if float(line_split[1].split(';')[0]) > price_cutoff:
                        continue
                    else:
                        try:
                            rdkit_mol = Chem.MolFromSmiles(line_split[0])
                            rdkit_smiles = Chem.MolToSmiles(rdkit_mol)
                        except:
                            continue
                        if line_split[0] not in reagents_dataset.keys():
                            reagents_dataset[line_split[0]] = {'prices': [], 'mol': rdkit_mol}
                        for reagent_price in line_split[1::]:
                            reagent_price_split = reagent_price.split(';')
                            reagents_dataset[line_split[0]]['prices'].append([float(reagent_price_split[0]),
                                                                              *reagent_price_split[1::]])
                        reagents_dataset[line_split[0]]['prices'].sort(key=lambda x: x[0], reverse=False)
        if 'in_lab' in chemical_sources:
            print('\t\tAdding all of the in-lab chemicals to the reagents working set')
            with open(r'.\chemical_library\reagent_sources\inlab_chemical_list.txt', 'r') as infile:
                for line in infile:
                    line_split = line.strip('\n')
                    if line_split == '':
                        continue
                    if '@' in line_split:
                        continue
                    try:
                        rdkit_mol = Chem.MolFromSmiles(line_split)
                        rdkit_smiles = Chem.MolToSmiles(rdkit_mol)
                    except:
                        continue
                    if line_split not in reagents_dataset.keys():
                        reagents_dataset[line_split] = {'prices': [], 'mol': rdkit_mol}
                    reagents_dataset[line_split]['prices'].append([1, '', 'InLab'])
                    reagents_dataset[line_split]['prices'].sort(key=lambda x: x[0], reverse=False)
        if 'chemistry' in chemical_sources:
            print('\t\tAdding all of the chemistry chemicals to the reagents working set')
            with open(r'.\chemical_library\reagent_sources\chemistry_chemicals_cleaned.yaml', 'r') as yamlfile:
                chemistry_inventory = yaml.load(yamlfile, Loader=yaml.Loader)
        else:
            chemistry_inventory = {}
    else:
        sys.exit('Reagent origin of %s is not implemented' % reagent_sources['reagent_origin'])
    banned_chemicals = load_banned_chemicals()
    for chemical in banned_chemicals:
        if chemical in reagents_dataset.keys():
            del reagents_dataset[chemical]
    return ['Success', reagents_dataset, chemistry_inventory]
            

def load_inlab_chemicals():
    reagents_dataset = set()
    with open(r'.\chemical_library\reagent_sources\inlab_chemical_list_general.txt', 'r') as infile:
        for line in infile:
            line_split = line.strip('\n')
            if line_split == '':
                continue
            try:
                rdkit_mol = Chem.MolFromSmiles(line_split)
                rdkit_smiles = Chem.MolToSmiles(rdkit_mol)
            except:
                continue
            reagents_dataset.add(line_split)
    with open(r'.\chemical_library\reagent_sources\banned_chemicals.txt', 'r') as infile:
        for line in infile:
            line_split = line.strip('\n')
            if line_split == '':
                continue
            if line_split in reagents_dataset:
                reagents_dataset.remove(line_split)
    return reagents_dataset


def load_banned_chemicals():
    banned_chemicals = []
    with open(r'.\chemical_library\reagent_sources\banned_chemicals.txt', 'r') as infile:
        for line in infile:
            banned_chemicals.append(line.strip())
    return banned_chemicals








