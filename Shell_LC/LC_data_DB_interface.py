#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 11:32:05 2022

database_interface.py but only meant to be accessed by LC
  The interface is slightly different and we never got around to merging them

@author: mattmcdonald
"""

import datetime
import json
import re
from mcn_status import MCN_CFG, S_DATABASES
import requests

DATABASE_URL, DATABASE_PORT = MCN_CFG[S_DATABASES]['Data']


def change_network_settings(url, port):
    global DATABASE_URL, DATABASE_PORT
    DATABASE_URL = str(url)
    DATABASE_PORT = str(int(port))


def data_DB_request(request: dict):
    """ Convenience function for making database requests

    """
    incoming_request = {'auth': {'port': DATABASE_PORT, 'user': 'xxxxxxxx', 'password': 'xxxxxxxx'},
                        'request': request}

    try:
        resp = requests.post(DATABASE_URL, json=incoming_request, headers={'content-type': 'application/authjson'})
        return json.loads(resp.content.decode('utf-8'))
    except:
        print("data DB error, either bad connection or bad request")
        return "Error"


def query_reaction(reaction_smiles: str):
    """
    query_reaction returns the reaction dictionary from the database

    Parameters
    ----------
    reaction_smiles : str
        the reaction smiles

    Returns
    -------
    reaction_dict.

    """
    request = {'request_type': 'query_data',
               'search_type': 'reaction',
               'search_term': ['reaction_smiles', reaction_smiles]}
    db_return = data_DB_request(request)
    # TODO: check that the success return is capitalized
    if db_return[0] == 'Success':
        return db_return[1]
    else:
        print("Reaction query did not return success:")
        print(db_return[0])
        return db_return[1]


def add_reaction(contents: dict, reaction_info: dict, operations: dict):
    """
    Adds a reaction to the data database

    Parameters
    ----------
    contents : dict
        standard contents dictionary
    reaction_info : dict
        reaction information as returned by lc_data_analysis
    operations : dict
        operations from the queue document used to perform the reaction

    Returns
    -------
    Result of database request.... but maybe should turn this into a simple boolean

    """
    reaction_smiles = contents['reaction_smiles']        
    reactant_smiles = reaction_smiles.split('>>')[0].split('.')
    product_smiles = contents['target_product'][0][0]
    reagent_smiles = contents['predicted_reagents']
    if not reagent_smiles:
        reagent_smiles = []
    catalyst_smiles = contents['predicted_catalyst']
    if not catalyst_smiles:
        catalyst_smiles = []

    if 'templates' in contents.keys():
        templates = contents['templates']
    else:
        templates = []
    reactants = []
    reagents = []
    catalysts = []
    for reagent in contents['reagents']:
        if reagent[2] in reactant_smiles:
            reactants.append(reagent)
        elif reagent[2] in reagent_smiles:
            reagents.append(reagent)
        elif reagent[2] in catalyst_smiles:
            catalysts.append(reagent)
            
    method = 'LC PDA estimate... in development'
    product_type = 'final' if contents['final_product'][1] == 'yes' else 'intermediate'
    product_info = [[product_smiles, reaction_info['yield'], method, product_type]]
    
    temperature = 25.0
    duration = ''
    location = ''
    first_heat_shake = True
    for i in range(len(operations)):
        op = operations[str(i+1)]
        if op['operation'] == 'start_stop_heater_shaker' and first_heat_shake:
            match = re.search('reaction_plate_*.*', op['container'])
            if match:
                if op['details']['power'] == 'on':
                    temperature = op['details']['temperature']
                    try:
                        start_time = datetime.datetime.strptime(op['start_time'], '%m/%d/%Y %H:%M:%S')
                    except:
                        start_time = datetime.datetime.now()
                    location = 'AMD_heater_shaker'
                elif op['details']['power'] == 'off':
                    try:
                        end_time = datetime.datetime.strptime(op['start_time'], '%m/%d/%Y %H:%M:%S')
                    except:
                        end_time = datetime.datetime.now()
                    first_heat_shake = False
        if op['operation'] == 'Th.run':
            temperature = op['details']['heating_profile']['heat_temperature']
            duration = op['details']['heating_profile']['heat_hold_time']
            duration = datetime.timedelta(seconds=duration)
            location = 'AMD_thermal_reactor'
            break
    if not duration:
        try:
            duration = end_time - start_time
        except:
            print('There appear to be no reaction steps in this queue')
    duration = str(duration).split('.')[0]
    
    # currently use metadata to give filepaths to HPLCMS data (MOCCA reports)
    metadata = {'PDA_Data': reaction_info['PDA_Data'],
                'MS_Data': reaction_info['MS_Data'],
                'templates': reaction_info['templates']}
    
    db_formatted_details = {'reaction_reactants': reactants,
                            'reaction_reagents': reagents,
                            'reaction_catalysts': catalysts,
                            'reaction_products': product_info,
                            'reaction_solvents': contents['solvents'],
                            'reaction_temperature': temperature,
                            'reaction_duration': duration,
                            'reaction_template': templates,
                            'reaction_location': ['measured', location],
                            'reaction_date': datetime.datetime.now().strftime('%m/%d/%Y'),
                            'reaction_metadata': metadata
                            }
    
    my_request = {'request_type': 'add_reaction',
                  'reaction_type': 'synthesis',
                  'reaction_smiles': reaction_smiles,
                  'reaction_details': db_formatted_details}
    db_return = data_DB_request(my_request)
    if db_return[0] == 'Success':
        return True
    else:
        print("issue adding reaction to data DB:")
        print(db_return[1])
        return False


def add_extinction_coef_prediction(smiles: str, value: dict, campaign_name: str):
    """
    add an extinction coefficient prediction to the data DB

    Parameters
    ----------
    smiles : str
        SMILES string for predicted molecule.
    value : dict
        A dictionary with keys = solvents (must include acetonitrile) and values = [log_E, log_E unc].
    campaign_name : str
        the name of the campaign that generated the data

    Returns
    -------
    The data DB return for adding data, not sure what that looks like yet but probably a tuple

    """
    if 'CC#N' not in value.keys():
        print("Acetonitrile (CC#N) not in log_e values!")
        return False
    my_request = {'request_type': 'add_data',
                  'molecule': [smiles, 'smiles'],
                  'add_data': True,
                  'data_to_add': {'data_name': smiles+' ext_coef_prediction',
                                  'type_tag': 'calc_loge_chemprop',
                                  'data_origin': ['calculated', 'chemprop'],
                                  'origin_date': datetime.datetime.now().strftime('%m/%d/%Y'),
                                  'loge_value': value,
                                  'metadata': {}}}
    request_return = data_DB_request(my_request)
    if request_return[0] == 'Success':
        return True
    else:
        return request_return[1]


def add_HPLC_logD_measurement(smiles: str, value: list, campaign_name: str, metadata=None):
    """
    add a measured logP value to the data DB

    Parameters
    ----------
    smiles : str
        SMILES string for predicted molecule.
    value : list
        [logp, uncertainty]
    campaign_name : str
        the name of the campaign that generated the data

    Returns
    -------
    The data DB return for adding data, not sure what that looks like yet but probably a tuple.
    """
    if not metadata:
        metadata = {}
    print(smiles)
    my_request = {'request_type': 'add_data',
                  'molecule': [smiles, 'smiles'],
                  'add_data': True,
                  'data_to_add': {'data_name': '%s_HPLC_logP_measurement' % smiles,
                                  'campaign_name': campaign_name,
                                  'type_tag': 'exp_logp_hplc',
                                  'data_origin': ['measured', 'amd_platform_HPLC'],
                                  'origin_date': datetime.datetime.now().strftime('%m/%d/%Y'),
                                  'logp_value': [[7.0, value[0], value[1]]],
                                  'solvent': ['O','CC#N'],
                                  'metadata': metadata}}
    request_return = data_DB_request(my_request)
    if request_return[0] == 'Success':
        return True
    else:
        return request_return[1]


def test_connection():
    test_request = {'request_type': 'add_reaction_bad',
                    'reaction_type': 'synthesis',
                    'reaction_smiles': 'garbage',
                    'reaction_details': 'garbage'}
    print(data_DB_request(test_request))


if __name__ == '__main__':
    test_connection()
