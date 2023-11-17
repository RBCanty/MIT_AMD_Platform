# -*- coding: utf-8 -*-
"""
Some code built to query the AMD Database
Created on Thu Dec 19 19:48:10 2019
Version 1
@author: Brent
"""
import pymongo
import sys
import re
from jsonschema import validate
from data_entry import rdkit_id_sanitize
from schema_container import get_schema
from data_entry import without_keys
from pprint import pprint
from bson.objectid import ObjectId

# -------------------------------------------------------------------------
# Below is the query_data subroutines for the AMD Database
def query_data_mongo(incoming_request, mongo_database):
    code_returns = {}
    # Brings in the mongo collection names for the AMD Database and initializes them
    myclient = mongo_database
    amddatabase = myclient['amd_database']
    cmpddatabase = amddatabase['compound_database']
    datadatabase = amddatabase['data_database']
    rxndatabase = amddatabase['reaction_database']
    # -------------------------------------------------------------------------
    code_returns['print_statements'] = []
    
    #schema_check = validate(instance = incoming_request, schema = get_schema('query_data_database'))
    #if schema_check is not None:
    #    code_returns['print_statements'].append(['Incoming request not understood'])
        #return code_returns
    try:
        search_type = incoming_request['search_type']
    except KeyError:
        code_returns['print_statements'].append(['Incoming request not understood please check the request Keys'])
        return ['Error', code_returns]
    if search_type == 'molecule':
        try:
            incoming_molecule = incoming_request['molecule']
            return_data_check = incoming_request['return_data']
        except KeyError:
            code_returns['print_statements'].append(['Incoming request not understood please check the request Keys'])
            return ['Error', code_returns]
        incoming_molecule_data = rdkit_id_sanitize(incoming_molecule[0],incoming_molecule[1])
        if type(incoming_molecule_data) is not dict:
            code_returns['print_statements'].append(['Something went wrong with your entry and rdkit: ' \
                        + str(incoming_molecule[0]) + ', ' + str(incoming_molecule[1])])
            code_returns['print_statements'].append([incoming_molecule_data])
            return ['Error', code_returns]
        else:
            file_count =  myclient.amd_database.compound_database.count_documents({"compound_inchi_key": \
                                                                                   incoming_molecule_data['compound_inchi_key']})
            if file_count == 0:
                code_returns['print_statements'].append('No database entries exist for: ' \
                            + incoming_molecule_data['compound_inchi_key'])
                return ['Error', code_returns]
            elif file_count == 1:
                code_returns['print_statements'].append('A database entry exists for: ' \
                            + incoming_molecule[0] + ' -> ' + incoming_molecule_data['compound_inchi_key'])
            else:
                code_returns['print_statements'].append('Not sure what to do about: ' + incoming_molecule_data['compound_inchi_key']\
                            + ' -> ' + incoming_molecule[0] + ' (Entries = ' + str(file_count) + ')')
                return ['Error', code_returns]
            database_entry = myclient.amd_database.compound_database.find_one({"compound_inchi_key": \
                                                                               incoming_molecule_data['compound_inchi_key']})
            code_returns['compound_database_entry'] = without_keys(database_entry,['_id'])
            code_returns['compound_database_entry']['_id'] = str(database_entry['_id'])
            if return_data_check == True:
                code_returns['print_statements'].append('Looking for data entries for ' + incoming_molecule[0]\
                            + ' -> ' + incoming_molecule_data['compound_inchi_key'])
                data_entries = myclient.amd_database.data_database.find({"compound_inchi_key": incoming_molecule_data['compound_inchi_key']})
                index = 1
                code_returns['data_database_entry'] = {}
                for element in data_entries:
                    #print(element)
                    code_returns['data_database_entry'][element['data_id_string']] = without_keys(element,['_id'])
                    code_returns['data_database_entry'][element['data_id_string']]['id'] = str(element['_id'])
                    index = index + 1
                if index == 1:
                    code_returns['print_statements'].append('No data entries for ' + incoming_molecule_data['compound_inchi_key'])
                    code_returns['data_database_entry'] = None
                else:
                    code_returns['print_statements'].append('Gathered data entries for ' + incoming_molecule_data['compound_inchi_key']\
                                + ' -> ' + str(index - 1) + ' entries')
            return ['Success', code_returns]
        
        
    elif search_type == 'special_search':
        try:
            if type(incoming_request['search_term']) == list:
                if incoming_request['search_term'][0].lower() == '$exists':
                    search_term = {"$exists" : bool(incoming_request['search_term'][1])}
                elif incoming_request['search_term'][0].lower() == '$regex':
                    search_term = {"$regex" : str(incoming_request['search_term'][1])}
                elif len(incoming_request['search_term']) == 1:
                    search_term = str(incoming_request['search_term'][0])
            elif type(incoming_request['search_term']) == dict:
                search_term = incoming_request['search_term']
            else:
                code_returns['print_statements'].append('Unable to interpret incoming search term')
                return ['Error', code_returns]
            try:
                if incoming_request['collection'] == 'data':
                    cursor_object = myclient.amd_database.data_database.find({incoming_request['document_property'] : \
                                                                              search_term},\
                                                                            {'compound_inchi_key': 1, 'data_id_string': 1,\
                                                                             'data_entry.type_tag': 1, 'data_entry.compound_inchi': 1,\
                                                                             'data_entry.compound_smiles': 1})
                elif incoming_request['collection'] == 'compound':
                    cursor_object = myclient.amd_database.compound_database.find({incoming_request['document_property'] : \
                                                                              search_term},\
                                                                            {'compound_inchi_key': 1, 'data_id_string': 1,\
                                                                             'data_entry.type_tag': 1, 'data_entry.compound_inchi': 1,\
                                                                             'data_entry.compound_smiles': 1}) 
                elif incoming_request['collection'] == 'reaction':
                    cursor_object = myclient.amd_database.reaction_database.find({incoming_request['document_property'] : \
                                                                              search_term},\
                                                                            {'reaction_smiles': 1, 'data_id_string': 1,\
                                                                             'data_entry.type_tag': 1, 'data_entry.compound_inchi': 1,\
                                                                             'data_entry.compound_smiles': 1}) 
                else:
                    code_returns['print_statements'].append('Search term (' + str(incoming_request['collection']) + ') not recognized')
                    return ['Error', code_returns]
            except UnboundLocalError:
                code_returns['print_statements'].append('Unable to interpret incoming search term')
                return ['Error', code_returns]
            raw_search_results = [element for element in cursor_object]
            search_results = []
            for item in raw_search_results:
                temp = without_keys(item,['_id'])
                temp['_id'] = str(item['_id'])
                search_results.append(temp)
            if len(search_results) == 0:
                code_returns['print_statements'].append('No entries consistent with database property search, check search parameters')
                return ['Error', code_returns]
            else:
                code_returns['print_statements'].append('Entries consistent with database special molecule search were found: ' +\
                            str(len(search_results)))
                code_returns['search_results'] = search_results
                return ['Success', code_returns]
        except KeyError:
            code_returns['print_statements'].append('Request search keys were not properly sent, please check formatting')
            return ['Error', code_returns]

    elif search_type == 'object_id':
        search_term = incoming_request['object_id']
        try:
            if incoming_request['collection'] == 'data':
                cursor_object = myclient.amd_database.data_database.find({'_id' : ObjectId(str(search_term))})
            elif incoming_request['collection'] == 'compound':
                cursor_object = myclient.amd_database.compound_database.find({'_id' : ObjectId(str(search_term))})
            elif incoming_request['collection'] == 'reaction':
                cursor_object = myclient.amd_database.reaction_database.find({'_id' : ObjectId(str(search_term))})
            else:
                code_returns['print_statements'].append('Search term (' + str(incoming_request['collection']) + ') not recognized')
                return ['Error', code_returns]
            raw_search_results = [element for element in cursor_object]
            search_results = []
            for item in raw_search_results:
                temp = without_keys(item,['_id'])
                temp['_id'] = str(item['_id'])
                search_results.append(temp)
            if len(search_results) == 0:
                code_returns['print_statements'].append('No entries consistent with object id search, check search parameters')
                return ['Error', code_returns]
            else:
                code_returns['print_statements'].append('Entries consistent with object id search were found: ' +\
                            str(len(search_results)))
                code_returns['search_results'] = search_results
                return ['Success', code_returns]
        except KeyError:
            code_returns['print_statements'].append('Request search keys were not properly sent, please check formatting')
            return ['Error', code_returns]
    elif search_type == 'reaction':
        if 'search_term' not in incoming_request:
            code_returns['print_statements'].append('Reaction search key were not properly sent, please check formatting')
            return ['Error', code_returns]
        if incoming_request['search_term'][0] == 'reaction_smiles':
            split_smiles = incoming_request['search_term'][0].split('>>')
            reagent_list = []
            for i in split_smiles[0].split('.'):
                rdkit_cleaned = rdkit_id_sanitize(i, 'smiles')
                reagent_list.append(rdkit_cleaned['basic_properties']['compound_smiles'])
            product_list = []
            for i in split_smiles[1].split('.'):
                rdkit_cleaned = rdkit_id_sanitize(i, 'smiles')
                product_list.append(rdkit_cleaned['basic_properties']['compound_smiles'])
            cleaned_search_term = '.'.join(reagent_list) + '>>' + '.'.join(product_list)
            reaction_entries = myclient.amd_database.reaction_database.find(\
                                       {incoming_request['search_term'][1] : cleaned_search_term})
        elif incoming_request['search_term'][1] == 'reagents' or incoming_request['search_term'][1] == 'products':
            cleaned_search_term = rdkit_id_sanitize(incoming_request['search_term'][0], 'smiles')       
            reaction_entries = myclient.amd_database.reaction_database.find(\
                                       {incoming_request['search_term'][1] : cleaned_search_term['basic_properties']['compound_smiles']})
        #rxndatabase.delete_many({incoming_request['search_term'][0] : cleaned_search_term['basic_properties']['compound_smiles']})
        index = 1
        code_returns['reaction_database_entry'] = {}
        for element in reaction_entries:
            code_returns['reaction_database_entry'][element['reaction_smiles']] = {}
            #print(element.keys())
            #print(element)
            for item in element.keys():
                if item != '_id':
                    code_returns['reaction_database_entry'][element['reaction_smiles']][item] = element[item]
                elif item == '_id':
                    code_returns['reaction_database_entry'][element['reaction_smiles']]['id'] = str(element['_id'])
            #code_returns['reaction_database_entry'][element['reaction_smiles']] = without_keys(element,['_id'])
            #code_returns['reaction_database_entry'][element['reaction_smiles']]['id'] = str(element['_id'])
            index = index + 1
        if index == 1:
            code_returns['print_statements'].append('No data entries for ' + '='.join(incoming_request['search_term']))
            code_returns['reaction_database_entry'] = None
        else:
            code_returns['print_statements'].append('Gathered reaction entries for ' + '='.join(incoming_request['search_term'])\
                        + ' -> ' + str(index - 1) + ' entries')
        return ['Success', code_returns]
    
    
    elif search_type == 'model_information':
        if 'search_term' not in incoming_request:
            code_returns['print_statements'].append('Model information search not properly sent')
            return ['Error', code_returns]
        
        if incoming_request['search_term'][0] == 'model_type':
            cursor_object = myclient.amd_database.model_information.find({incoming_request['search_term'][0]:
                                                                            incoming_request['search_term'][1]})
            code_returns['model_information_entry'] = {}
            for element in cursor_object:
                code_returns['model_information_entry'][str(element['_id'])] = element
                code_returns['model_information_entry'][str(element['_id'])]['_id'] = str(element['_id'])
            if len(list(code_returns['model_information_entry'].keys())) == 0:
                code_returns['print_statements'].append('No new model information documents for %s' % ' '.join(incoming_request['search_term']))
            else:
                code_returns['print_statements'].append('Gathered %s model information documents for %s' % (
                        len(list(code_returns['model_information_entry'].keys())), incoming_request['search_term']))
            return ['Success', code_returns]
        else:
            code_returns['print_statements'].append('Search term %s not recognized' % incoming_request['search_term'][0])
            return ['Error', code_returns]
        
    
    else:
        code_returns['print_statements'].append('Search type identifier not understood')
        return ['Error', code_returns]

