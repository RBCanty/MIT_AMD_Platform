# -*- coding: utf-8 -*-
"""
Some code built to add data to the AMD Database (Compounds and Associated Data)
Created on Tue Dec 17 15:16:00 2019
Version 1
@author: Brent
"""
import pymongo
import sys
from jsonschema import validate
from check_pubchem import get_pubchem_data
from data_entry import spectral_data_entry
from data_entry import kinetic_database_entry
from data_entry import update
from data_entry import rdkit_id_sanitize
from data_entry import logp_database_entry
from data_entry import loge_database_entry
from data_entry import without_keys
from schema_container import get_schema
from pprint import pprint
import jsonschema

# -------------------------------------------------------------------------
# Below is the add_data subroutines for the AMD Database
def add_data_mongo(incoming_request, mongo_database):
    code_returns = {}
    # Brings in the mongo collection names for the AMD Database and initializes
    # them for use in the code
    myclient = mongo_database
    amddatabase = myclient['amd_database']
    cmpddatabase = amddatabase['compound_database']
    datadatabase = amddatabase['data_database']
    # -------------------------------------------------------------------------
    code_returns['print_statements'] = []
    
    # Separates the incoming request inschemato each of the components for the code
    # to use, returns an error if required keys are not present in the request
    # For add_data we need the molecule and if you want to add data to the entry
    try:
        incoming_molecule = incoming_request['molecule']
        if 'add_data' in incoming_request.keys():
            add_data_to_database = incoming_request['add_data']
            if incoming_request['add_data'] == True:
                try:
                    new_data_entry = incoming_request['data_to_add']
                except KeyError:
                    code_returns['print_statements'].append(['No data was specified in add_data request but add_data = True'])
                    return ['Error', code_returns]
            else:
                code_returns['print_statements'].append(['No data to add to the molecule request'])
        else:
            add_data_to_database = False
    except KeyError:
        code_returns['print_statements'].append(['Incoming request not understood, appropriate keys not supplied'])
        return ['Error', code_returns]
    if 'pubchem_check' in incoming_request.keys():
        pubchem_check = incoming_request['pubchem_check']
    else:
        pubchem_check = False

    
    # First sanitize entry with rdkit to ensure uniformity and then check to see
    # if that compound is already in the database, add basic data if it is not
    incoming_molecule_data = rdkit_id_sanitize(incoming_molecule[0],incoming_molecule[1])
    if type(incoming_molecule_data) is not dict:
        code_returns['print_statements'].append(['Something went wrong with your entry and rdkit: ' \
                    + str(incoming_molecule[0]) + ', ' + str(incoming_molecule[1])])
        code_returns['print_statements'].append([incoming_molecule_data])
        return ['Error', code_returns]
    else:
        # Search through the existing compound database to see if there is an
        # entry for the incoming request already by searching for an InChI Key
        file_count =  myclient.amd_database.compound_database.count_documents({"compound_inchi_key": \
                                                                               incoming_molecule_data['compound_inchi_key']})
        if file_count == 0:
            code_returns['print_statements'].append('Inserting a new compound stub for: ' \
                        + incoming_molecule_data['compound_inchi_key'])
            database_entry = cmpddatabase.insert_one(incoming_molecule_data)
        elif file_count == 1:
            code_returns['print_statements'].append('A database document exists for: ' + \
                        incoming_molecule[0] + ' -> ' + incoming_molecule_data['compound_inchi_key'])
        else:
            code_returns['print_statements'].append('Not sure what to do about: ' + \
                        incoming_molecule[0] + ' -> ' + incoming_molecule_data['compound_inchi_key'] + \
                        ' (Entries = ' + str(file_count) + ')')
            return ['Error', code_returns]
    
    
    # If you want to pull some experimental properties about the molecule from PubChem
    # include a statement 'pubchem_check' = True in the http request
    if pubchem_check == True:
        code_returns['print_statements'].append('Grabbed data from pubchem for: ' + incoming_molecule[0] + \
                    ' -> ' + incoming_molecule_data['compound_inchi_key'])
        existing_dict = myclient.amd_database.compound_database.find_one({"compound_inchi_key": \
                                                                          incoming_molecule_data['compound_inchi_key']})
        pubchem_dict = get_pubchem_data(incoming_molecule_data['compound_inchi_key'])
        update(existing_dict, pubchem_dict)
        # This schema check was put here just in case but not really necessary
        # it is most useful for troubleshooting the code really to make sure that
        # things remain in the predefined schema when messing with PubChem Entries
        initial_schema = get_schema('initial_entry')
        try:
            schema_check = validate(instance = existing_dict, schema = initial_schema)
            code_returns['print_statements'].append('Database entry for ' + incoming_molecule_data['compound_inchi_key'] \
                        + ' -> ' + incoming_molecule[0] + ' complies with schema')
        except jsonschema.exceptions.ValidationError:
            code_returns['print_statements'].append('Database entry for ' + incoming_molecule[0] + ' does not comply with schema')
            return ['Error', code_returns]
        # Basically as long as it is successful it moves on which is almost always
        # since the code and database build consistent entries
        cmpddatabase.update_one({"compound_inchi_key": incoming_molecule_data['compound_inchi_key']}\
                                 ,{"$set": existing_dict})
        database_entry = myclient.amd_database.compound_database.find_one({"compound_inchi_key": \
                                                                           incoming_molecule_data['compound_inchi_key']})    
    else:
        database_entry = myclient.amd_database.compound_database.find_one({"compound_inchi_key": \
                                                                           incoming_molecule_data['compound_inchi_key']})
        initial_schema = get_schema('initial_entry')
        try:
            schema_check = validate(instance = database_entry, schema = initial_schema)
            code_returns['print_statements'].append('Database entry for ' + incoming_molecule_data['compound_inchi_key'] \
                        + ' -> ' + incoming_molecule[0] + ' complies with schema')
        except jsonschema.exceptions.ValidationError:
            code_returns['print_statements'].append('Database entry for ' + incoming_molecule[0] + ' does not comply with schema')
            return ['Error', code_returns]

    code_returns['initial_database_entry'] = without_keys(database_entry,['_id'])
    code_returns['initial_database_entry']['_id'] = str(database_entry['_id'])
    
    #pprint(database_entry)
    # -------------------------------------------------------------------------
    # With the basic compound entry completed, it is time to add specific data
    # to the compound entry which is selected based on pre-defined tags
    if add_data_to_database == True:
        # Here is the type_tag selector to go through spectral addition or logP
        # addition because each type of entry has their own schema
        # TODO: Include entry rules/schema for other types of data
        spec_rules = ['_abs_s' in new_data_entry['type_tag'],
                      '_pl_' in new_data_entry['type_tag'],
                      '_ir_' in new_data_entry['type_tag']]
        logp_rules = ['_logp_' in new_data_entry['type_tag']]
        extinct_rules = ['_loge_' in new_data_entry['type_tag']]
        kinetic_rules = ['_kinetic_' in new_data_entry['type_tag']]
        
        if any(spec_rules):
            new_data_struct, data_storage_struct = spectral_data_entry(new_data_entry, database_entry)
        elif any(logp_rules):
            new_data_struct, data_storage_struct = logp_database_entry(new_data_entry, database_entry)
        elif any(extinct_rules):
            new_data_struct, data_storage_struct = loge_database_entry(new_data_entry, database_entry)
        elif any(kinetic_rules):
            new_data_struct, data_storage_struct = kinetic_database_entry(new_data_entry, database_entry)
        else:
            code_returns['print_statements'].append('Specified data type tag ' +  new_data_entry['type_tag'] + ' not recognized')
            return ['Error', code_returns]
        try:
            if type(new_data_struct) is not dict:
                code_returns['print_statements'].append(str(new_data_struct) + ' for ' + new_data_entry['data_name'] + \
                            ' in ' + incoming_molecule[0] + ' -> ' + incoming_molecule_data['compound_inchi_key'])
                if isinstance(data_storage_struct, str) and data_storage_struct == 'duplicate':
                    return ['Success', code_returns]
                return ['Error', code_returns]
            else:
                code_returns['print_statements'].append('Adding data entry for ' + new_data_entry['data_name'] +\
                            ' in ' + incoming_molecule_data['compound_inchi_key'] + ' -> ' + incoming_molecule[0])
                cmpddatabase.find_one_and_update({"compound_inchi_key": incoming_molecule_data['compound_inchi_key']},\
                                                  {"$set": new_data_struct})
                file_counter = myclient.amd_database.data_database.count_documents({"data_id_string": data_storage_struct['data_id_string']})
                if file_counter == 0:
                    datadatabase.insert_one(data_storage_struct)
                elif file_counter == 1:
                    datadatabase.update_one({"data_id_string": data_storage_struct['data_id_string']},\
                                             {"$set": data_storage_struct})
                else:
                    code_returns['print_statements'].append('Not sure what to do about: ' + data_storage_struct['data_id_string'] +\
                                ' -> ' + incoming_molecule[0] + ' (Entries = ' + str(file_count) + ')')
                    return ['Error', code_returns]
                database_entry = myclient.amd_database.data_database.find_one({"data_id_string": data_storage_struct['data_id_string']})
                code_returns['data_entry'] = without_keys(database_entry,['_id'])
                code_returns['data_entry']['_id'] = str(database_entry['_id'])
        except NameError:
            code_returns['print_statements'].append('Data type not understood: ' + new_data_entry['type_tag'] +\
                        ' -> ' + incoming_molecule[0])
            return ['Error', code_returns]
    else:
        code_returns['print_statements'].append('No data to add to the database entry: ' \
                    + incoming_molecule_data['compound_inchi_key'])
        return ['Error', code_returns]
    database_entry = myclient.amd_database.compound_database.find_one({"compound_inchi_key": incoming_molecule_data['compound_inchi_key']})
    code_returns['updated_database_entry'] = without_keys(database_entry,['_id'])
    code_returns['updated_database_entry']['_id'] = str(database_entry['_id'])
    return ['Success', code_returns]    

