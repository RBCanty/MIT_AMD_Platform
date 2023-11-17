# -*- coding: utf-8 -*-
"""
Extra functions that are useful for the AMD Database
Created on Tue Dec 17 15:43:11 2019

@author: Brent
"""
import datetime
import math
import string
import collections
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from schema_container import get_schema
from jsonschema import validate
import jsonschema
import copy
import re


def without_keys(d,keys):
    return {k: v for k, v in d.items() if k not in keys}


def update(dct, merge_dct):
    for k, v in merge_dct.items():
        if (k in dct and isinstance(dct[k], dict)
                and isinstance(merge_dct[k], collections.Mapping)):
            update(dct[k], merge_dct[k])
        else:
            dct[k] = merge_dct[k]
    return dct        


def rdkit_id_sanitize(cmpd_id, id_type):
    # Here take the information from the incoming request and give it to RDKit
    # to sanitize the incoming string
    if str(id_type) == 'smiles':
        rdk_mol = Chem.MolFromSmiles(cmpd_id)
    elif str(id_type) == 'inchi':
        rdk_mol = Chem.MolFromInchi(cmpd_id)     
    else:
        return str('Molecule id_type (' + str(id_type) + ') for incoming request not recognized')
    if rdk_mol == None:
        return str('Incoming Molecule (' + str(cmpd_id) + ') is not usable in RDKit as a ' + str(id_type) + ' entry')
    create_date = datetime.date.today()
    kept_data = {}
    kept_data['compound_inchi_key'] = Chem.MolToInchiKey(rdk_mol)
    kept_data['date_created'] = create_date.strftime('%m/%d/%Y')
    kept_data['last_modified'] = create_date.strftime('%m/%d/%Y')
    kept_data['basic_properties'] = {}
    kept_data['basic_properties']['compound_smiles'] = Chem.MolToSmiles(rdk_mol)
    kept_data['basic_properties']['compound_inchi'] = Chem.MolToInchi(rdk_mol)
    kept_data['basic_properties']['compound_inchi_fixedh'] = Chem.MolToInchi(rdk_mol, options='/fixedH')
    kept_data['basic_properties']['compound_mformula'] = rdMolDescriptors.CalcMolFormula(rdk_mol)
    kept_data['basic_properties']['compound_mweight'] = round(Descriptors.MolWt(rdk_mol), 2)
    return kept_data


def spectral_data_entry(data_dict,exist_dict):
    # Right now the allowed rules switches between calculated and experimental datasets
    # Then check to make sure that one of type tags is present in the incoming request
    spectral_rules = ['exp_abs_spec' in data_dict['type_tag'],
                      'exp_pl_spec' in data_dict['type_tag'],
                      'exp_abs_stick' in data_dict['type_tag'],
                      'exp_pl_stick' in data_dict['type_tag'],
                      'exp_ir_spec' in data_dict['type_tag'],
                      'exp_ir_stick' in data_dict['type_tag']]
    calculated_spectral_rules = ['calc_abs_spec' in data_dict['type_tag'],
                        'calc_abs_stick' in data_dict['type_tag'],
                        'calc_pl_spec' in data_dict['type_tag'],
                        'calc_pl_stick' in data_dict['type_tag'],
                        'calc_ir_spec' in data_dict['type_tag'],
                        'calc_ir_stick' in data_dict['type_tag']]
    if any(spectral_rules):
        container = ['experimental_data','spectral_data']
    elif any(calculated_spectral_rules):
        container = ['calculated_data','spectral_data']
    else:
        return str('Specified type tag ' + data_dict['type_tag'] + ' is not recognized'),''
    current_time = datetime.datetime.now()
    data_dict['date_added'] = current_time.strftime('%m/%d/%Y')
    # Check to make sure that the incoming data request contains the required properties
    try:
        schema_check = validate(instance = data_dict, schema = get_schema('spectral_data'))
    except jsonschema.exceptions.ValidationError as e:
        return e, 'Missing required keys to add spectral data to the database'
    # Check to make sure that the incoming ydata is given as a list of lists
    # This can be changed to be more flexible but for consistent formatting this
    # statement block was left behind
    try:
        if type(data_dict['ydata'][0]) is list:
            data_name = re.sub(r'\D', "", "%s" % max(data_dict['ydata'][0])) if data_dict['ydata'][0] else ""
            data_id_string = data_dict['type_tag'] + '_' + \
            str(data_dict['origin_date']).translate(str.maketrans('','',string.punctuation)) + \
            '_' + data_name + exist_dict['basic_properties']['compound_smiles']
        else:
            return 'Please give ydata as a list of lists (even for one entry) for consistent formatting',''
    except KeyError:
        return 'Missing required keys to add spectral data to the database',''
    
    # Now we take the existing entry and determine if the dataset is already present
    # or not by checking the systematically generated data_id_strings
    try:
        spectra_number = len(exist_dict[container[0]][container[1]].items())            
    except KeyError:
        spectra_number = 0
    if spectra_number > 0:
        for element in exist_dict[container[0]][container[1]].keys():
            if data_id_string == element:
                return str('Entry ' + data_id_string + ' already exists'),'duplicate'
    else:
        exist_dict[container[0]] = {}
        exist_dict[container[0]][container[1]] = {}
    
    # Now we can prepare the data storage database entry which is the full entry
    # and prepare a dictionary without the xdata, ydata, and metadata to append
    # to the compound entry listing
    exist_dict[container[0]][container[1]][data_id_string] = without_keys(data_dict, ['xdata','ydata', 'metadata'])
    create_date = datetime.datetime.now()
    exist_dict['last_modified'] = create_date.strftime('%m/%d/%Y')
    storage_data = {}
    storage_data['data_id_string'] = data_id_string
    storage_data['data_entry'] = data_dict
    if "campaign_name" not in data_dict.keys():
        campaign_name = ''
    else:
        campaign_name = data_dict['campaign_name']
    storage_data['campaign_name'] = campaign_name
    storage_data['compound_inchi_key'] = exist_dict['compound_inchi_key']
    storage_data['data_entry']['compound_inchi'] = exist_dict['basic_properties']['compound_inchi']
    storage_data['data_entry']['compound_smiles'] = exist_dict['basic_properties']['compound_smiles']
    return (exist_dict, storage_data)


def kinetic_database_entry(data_dict, exist_dict):
    spectral_rules = ['exp_kinetic_abs' in data_dict['type_tag'],
                      'exp_kinetic_pl' in data_dict['type_tag']]
    calc_rules = ['calc_kinetic_abs' in data_dict['type_tag'],
                  'calc_kinetic_pl' in data_dict['type_tag']]
    if any(spectral_rules):
        container = ['experimental_data','spectral_data']
    elif any(calc_rules):
        container = ['calculated_data','spectral_data']
    else:
        return str('Specified type tag ' + data_dict['type_tag'] + ' is not recognized'),''
    current_time = datetime.datetime.now()
    data_dict['date_added'] = current_time.strftime('%m/%d/%Y')
    # Check to make sure that the incoming data request contains the required properties
    try:
        schema_check = validate(instance = data_dict, schema = get_schema('kinetic_data'))
    except jsonschema.exceptions.ValidationError:
        return 'Schema Failure Missing required keys to add kinetic data to the database', ''
    # Check to make sure that the incoming kinetic data is not redundant
    try:
        kinetic_rate = re.sub(r'\D', "", "%s" % data_dict['kinetic_rate'])
        data_id_string = data_dict['type_tag'] + '_' + \
            str(data_dict['origin_date']).translate(str.maketrans('','',string.punctuation)) \
            + '_' + kinetic_rate
    except KeyError as e:
        return 'Missing required keys to add kinetic data to the database', ''
    
    try:
        spectra_number = len(exist_dict[container[0]][container[1]].items())
    except KeyError:
        spectra_number = 0
    if spectra_number > 0:
        for element in exist_dict[container[0]][container[1]].keys():
            if data_id_string == element:
                return str('Entry ' + data_id_string + ' already exists'),'duplicate'
    else:
        exist_dict[container[0]] = {}
        exist_dict[container[0]][container[1]] = {}
    
    # Now we can prepare the data storage database entry which is the full entry
    # and prepare a dictionary without the xdata, ydata, and metadata to append
    # to the compound entry listing
    exist_dict[container[0]][container[1]][data_id_string] = data_dict
    create_date = datetime.datetime.now()
    exist_dict['last_modified'] = create_date.strftime('%m/%d/%Y')
    storage_data = {}
    storage_data['data_id_string'] = data_id_string
    storage_data['data_entry'] = data_dict
    if "campaign_name" not in data_dict.keys():
        campaign_name = ''
    else:
        campaign_name = data_dict['campaign_name']
    storage_data['campaign_name'] = campaign_name
    storage_data['compound_inchi_key'] = exist_dict['compound_inchi_key']
    storage_data['data_entry']['compound_inchi'] = exist_dict['basic_properties']['compound_inchi']
    storage_data['data_entry']['compound_inchi_fixedh'] = exist_dict['basic_properties']['compound_inchi_fixedh']
    storage_data['data_entry']['compound_smiles'] = exist_dict['basic_properties']['compound_smiles']
    return (exist_dict, storage_data)


def logp_database_entry(data_dict,exist_dict):
    exp_rules = ['exp_logp_abs' in data_dict['type_tag'],
                 'exp_logp_hplc' in data_dict['type_tag'],
                 'exp_logp_lit' in data_dict['type_tag']]
    calc_rules = ['calc_logp' in data_dict['type_tag']]
    if any(exp_rules):
        container = ['experimental_data','spectral_data']
    elif any(calc_rules):
        container = ['calculated_data','spectral_data']
    else:
        return str('Specified type tag ' + data_dict['type_tag'] + ' is not recognized'),''
    current_time = datetime.datetime.now()
    data_dict['date_added'] = current_time.strftime('%m/%d/%Y')
    # Check to make sure that the incoming data request contains the required properties
    try:
        schema_check = validate(instance = data_dict, schema = get_schema('logp_data'))
    except jsonschema.exceptions.ValidationError:
        return 'Schema Failure Missing required keys to add logP data to the database', ''
    # Check to make sure that the incoming logP value is a list and a redundant
    # check to make sure that there are no dictionary key errors
    try:
        if type(data_dict['logp_value']) is list:
            try:
                logp_list = copy.deepcopy(data_dict['logp_value'])
                logp_list.sort(key = lambda x: x[1], reverse = True)
                max_logp = logp_list[0][1]
                data_id_string = data_dict['type_tag'] + '_' + \
                str(data_dict['origin_date']).translate(str.maketrans('','',string.punctuation)) \
                    + '_' + '{:.12f}'.format(float(max_logp)).replace(".","").replace("-","")[:8]
            except TypeError:
                return 'List Failure Please give logP as a list of lists (even for one entry) for consistent formatting', ''
        else:
            return 'Please give logP as a list of lists (even for one entry) for consistent formatting', ''
    except KeyError as e:
        return 'Missing required keys to add logP data to the database', ''
    
    try:
        spectra_number = len(exist_dict[container[0]][container[1]].items())            
    except KeyError:
        spectra_number = 0
    if spectra_number > 0:
        for element in exist_dict[container[0]][container[1]].keys():
            if data_id_string == element:
                return str('Entry ' + data_id_string + ' already exists'),'duplicate'
    else:
        exist_dict[container[0]] = {}
        exist_dict[container[0]][container[1]] = {}
        
    # Now we can prepare the data storage database entry which is the full entry
    # and prepare a dictionary without the xdata, ydata, and metadata to append
    # to the compound entry listing
    exist_dict[container[0]][container[1]][data_id_string] = data_dict
    create_date = datetime.datetime.now()
    exist_dict['last_modified'] = create_date.strftime('%m/%d/%Y')
    storage_data = {}
    storage_data['data_id_string'] = data_id_string
    storage_data['data_entry'] = data_dict
    if "campaign_name" not in data_dict.keys():
        campaign_name = ''
    else:
        campaign_name = data_dict['campaign_name']
    storage_data['campaign_name'] = campaign_name
    storage_data['compound_inchi_key'] = exist_dict['compound_inchi_key']
    storage_data['data_entry']['compound_inchi'] = exist_dict['basic_properties']['compound_inchi']
    storage_data['data_entry']['compound_inchi_fixedh'] = exist_dict['basic_properties']['compound_inchi_fixedh']
    storage_data['data_entry']['compound_smiles'] = exist_dict['basic_properties']['compound_smiles']
    return (exist_dict, storage_data)


def loge_database_entry(data_dict,exist_dict):
    exp_rules = ['exp_loge_mass' in data_dict['type_tag']]
    calc_rules = ['calc_loge_chemprop' in data_dict['type_tag'],
                  'calc_loge_dft' in data_dict['type_tag']]
    if any(exp_rules):
        container = ['experimental_data','spectral_data']
    elif any(calc_rules):
        container = ['calculated_data','spectral_data']
    else:
        return str('Specified type tag ' + data_dict['type_tag'] + ' is not recognized'),''
    current_time = datetime.datetime.now()
    data_dict['date_added'] = current_time.strftime('%m/%d/%Y')
    # Check to make sure that the incoming data request contains the required properties
    try:
        schema_check = validate(instance = data_dict, schema = get_schema('loge_data'))
    except jsonschema.exceptions.ValidationError:
        return 'Schema Failure Missing required keys to add logP data to the database', ''
    # Check to make sure that the incoming logP value is a list and a redundant
    # check to make sure that there are no dictionary key errors
    try:
        acn_value = data_dict['loge_value']['CC#N'][0]
        data_id_string = data_dict['type_tag']+ '_' + \
                str(data_dict['origin_date']).translate(str.maketrans('','',string.punctuation)) \
                    + '_' + '{:.12f}'.format(float(acn_value)).replace(".","").replace("-","")[:8]
    except KeyError as e:
        return 'Missing required keys to add logE data to the database', ''
    
    try:
        spectra_number = len(exist_dict[container[0]][container[1]].items())            
    except KeyError:
        spectra_number = 0
    if spectra_number > 0:
        for element in exist_dict[container[0]][container[1]].keys():
            if data_id_string == element:
                return str('Entry ' + data_id_string + ' already exists'),'duplicate'
    else:
        exist_dict[container[0]] = {}
        exist_dict[container[0]][container[1]] = {}
        
    # Now we can prepare the data storage database entry which is the full entry
    # and prepare a dictionary without the xdata, ydata, and metadata to append
    # to the compound entry listing
    exist_dict[container[0]][container[1]][data_id_string] = data_dict
    create_date = datetime.datetime.now()
    exist_dict['last_modified'] = create_date.strftime('%m/%d/%Y')
    storage_data = {}
    storage_data['data_id_string'] = data_id_string
    storage_data['data_entry'] = data_dict
    if "campaign_name" not in data_dict.keys():
        campaign_name = ''
    else:
        campaign_name = data_dict['campaign_name']
    storage_data['campaign_name'] = campaign_name
    storage_data['compound_inchi_key'] = exist_dict['compound_inchi_key']
    storage_data['data_entry']['compound_inchi'] = exist_dict['basic_properties']['compound_inchi']
    storage_data['data_entry']['compound_inchi_fixedh'] = exist_dict['basic_properties']['compound_inchi_fixedh']
    storage_data['data_entry']['compound_smiles'] = exist_dict['basic_properties']['compound_smiles']
    return (exist_dict, storage_data)
