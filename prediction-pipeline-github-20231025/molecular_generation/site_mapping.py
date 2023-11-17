
import random

import numpy as np

from copy import deepcopy
from rdkit import Chem
from functions.supporting_functions import print_progress_bar


def normalize(item):
    if isinstance(item, list):
        gross = sum(item)
        return [i / gross for i in item]
    elif isinstance(item, dict):
        gross = sum(item.values())
        return {k: v / gross for k, v in item.items()}
    else:
        raise ValueError(f"Normalize expects a dict {{k: number}} or list [number] not {type(item)}")


def atom_map_scaffold(scaffolds_to_map, molecular_generation_details, derivative_sites=3, number_of_derivatives=20):
    if 'number_of_mapped_derivatives' in molecular_generation_details.keys():
        number_of_derivatives = molecular_generation_details['number_of_mapped_derivatives']
    if 'derivative_sites' in molecular_generation_details.keys():
        derivative_sites = molecular_generation_details['derivative_sites']
    
    special_sites = {'halogen': {'pattern': Chem.MolFromSmarts('[#9,#17,#35,#53]'),
                                 'value': 20,
                                 'substitution': '[#9,#17,#35,#53:2]'},
                     'carboxylic_acid': {'pattern': Chem.MolFromSmarts('[#8H1]-[#6]=O'),
                                         'value': 1,
                                         'substitution': '[#8H1]-[#6]=O'},
                     'carboxylic_acid_2': {'pattern': Chem.MolFromSmarts('[#8H1]-[#6]=O'),
                                           'value': 5,
                                           'substitution': '[#8H1]'},
                     'aldehyde': {'pattern': Chem.MolFromSmarts('[#1]-[#6]=O'),
                                  'value': 1,
                                  'substitution': '[#1]-[#6]=O'},
                     'alcohol': {'pattern': Chem.MolFromSmarts('[#6]-[#8H1]'),
                                 'value': 1,
                                 'substitution': '[#8H1]'},
                     'primary-amine': {'pattern': Chem.MolFromSmarts('[#6]-[#7H2]'),
                                       'value': 1,
                                       'substitution': '[#7H2]'},
                     'diazonium': {'pattern': Chem.MolFromSmarts('[$([#6]=[N+]=[N-]),$([#6-]-[N+]#[N])]'),
                                   'value': 50,
                                   'substitution': '[$([#6]=[N+]=[N-]),$([#6-]-[N+]#[N])]'},
                     'perfluoroalkylsulfonates': {
                         'pattern': Chem.MolFromSmarts('[#6]-[#8]S(=O)(=O)C([F,Cl,Br,I])([F,Cl,Br,I])[F,Cl,Br,I]'),
                         'value': 50,
                         'substitution': '[#8]S(=O)(=O)C([F,Cl,Br,I])([F,Cl,Br,I])[F,Cl,Br,I]'},
                     'nitro_1': {'pattern': Chem.MolFromSmarts('[#6][NX3+](=O)[O-]'),
                                 'value': 1,
                                 'substitution': '[NX3+](=O)[O-]'},
                     'nitro_2': {'pattern': Chem.MolFromSmarts('[#6][NX3](=O)=O'),
                                 'value': 1,
                                 'substitution': '[NX3](=O)=O'},
                     'boronic_acid': {'pattern': Chem.MolFromSmarts('[#6]-[#5;X3](-[#8H1])-[#8H1]'),
                                      'value': 1,
                                      'substitution': '[#5;X3](-[#8H1])-[#8H1]'},
                     }

    mapped_scaffolds_to_use = set()
    for scaffold_index, scaffold in enumerate(scaffolds_to_map):
        print_progress_bar(scaffold_index + 1, len(scaffolds_to_map), prefix='\t\tMapping: ', length=50)
        try:
            scaffold_mol = Chem.MolFromSmiles(scaffold)
            scaffold_smiles = Chem.MolToSmiles(scaffold_mol)
        except:
            continue

        # First we can look at the A-H sites
        mapping_sites = {}
        for atom in scaffold_mol.GetAtoms():
            if atom.GetImplicitValence() != 0:
                if atom.GetIdx() not in mapping_sites.keys():
                    mapping_sites[atom.GetIdx()] = {'value': [], 'type': []}
                value = 0
                if atom.GetAtomicNum() == 6:
                    if atom.GetIsAromatic():
                        value = 0.5
                    elif atom.GetNumImplicitHs() == 3:
                        value = 0.05
                    elif atom.GetNumImplicitHs() == 2:
                        value = 0.1
                    else:
                        value = 0.01
                elif atom.GetAtomicNum() == 7:
                    value = 1
                elif atom.GetAtomicNum() == 8:
                    value = 1
                elif atom.GetAtomicNum() == 16:
                    value = 1
                if value != 0:
                    mapping_sites[atom.GetIdx()]['value'].append(value)
                    mapping_sites[atom.GetIdx()]['type'].append('normal')

        # Then at a set of special sites with useful potential leaving groups
        for special_site in special_sites:
            current_site = special_sites[special_site]
            pattern_matches = scaffold_mol.GetSubstructMatches(current_site['pattern'])
            sub_matches = scaffold_mol.GetSubstructMatches(Chem.MolFromSmarts(current_site['substitution']))
            if len(sub_matches) == 0:
                continue
            for sub_match in sub_matches:
                for pattern_match in pattern_matches:
                    if set(sub_match).issubset(pattern_match):
                        break
                else:
                    continue
                other_atoms = set()
                current_atoms = list(sub_match)
                for current_atom in current_atoms:
                    atom = scaffold_mol.GetAtomWithIdx(current_atom)
                    connected_bonds = atom.GetBonds()
                    for connected_bond in connected_bonds:
                        bond_atoms = [connected_bond.GetBeginAtomIdx(),
                                      connected_bond.GetEndAtomIdx()]
                        other_atoms.update(bond_atoms)
                connecting_atoms = other_atoms.difference(set(current_atoms))
                if len(connecting_atoms) != 1:
                    continue
                other_atom = list(connecting_atoms)[0]
                if other_atom not in mapping_sites.keys():
                    mapping_sites[other_atom] = {'value': [], 'type': []}
                mapping_sites[other_atom]['value'].append(current_site['value'])
                mapping_sites[other_atom]['type'].append(special_site)

        potential_sites = list(mapping_sites.keys())
        weights = [max(mapping_sites[site]['value']) for site in mapping_sites.keys()]
        norm_weights = normalize(weights)
        mapped_scaffolds = set()
        for i in range(0, number_of_derivatives):
            num_sites = random.choice(range(1, derivative_sites+1))
            sites_to_map = np.random.choice(potential_sites, p=norm_weights, size=num_sites, replace=False)
            sites_list = list(map(int, sites_to_map))

            current_mol = deepcopy(scaffold_mol)
            for site in sites_list:
                current_mol.GetAtomWithIdx(site).SetAtomMapNum(1)
            for site in sites_list:
                sub_site = random.choices(mapping_sites[site]['type'], weights=mapping_sites[site]['value'])[0]
                if sub_site != 'normal':
                    substitution_type = special_sites[sub_site]
                    substitution = substitution_type['substitution']
                    mod_mol = Chem.DeleteSubstructs(current_mol,
                                                    Chem.MolFromSmarts(substitution.split('>>')[0]))
                    current_mol = deepcopy(mod_mol)
            mapped_scaffolds_to_use.add(Chem.MolToSmiles(current_mol))

    return ['Success', list(mapped_scaffolds_to_use)]
