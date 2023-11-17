# -*- coding: utf-8 -*-
"""
Some code built to add reactions to the AMD Database
Created on Mon Jan 13 12:47:27 2020
Version 0
@author: Brent
"""

import pymongo
from pprint import pprint
from data_entry import rdkit_id_sanitize
from data_entry import update
from schema_container import get_schema
from jsonschema import validate
import datetime
import re
import jsonschema
from copy import deepcopy


def add_reaction_mongo(incoming_request, mongo_database):
    code_returns = {}
    myclient = mongo_database
    amddatabase = myclient['amd_database']
    rxndatabase = amddatabase['reaction_database']
    cmpddatabase = amddatabase['compound_database']
    current_time = datetime.datetime.now()
    code_returns['print_statements'] = []
    
    if incoming_request['reaction_type'] == 'synthesis':
        reaction_smiles = incoming_request['reaction_smiles']
        #rxndatabase.delete_one({"reaction_smiles": reaction_smiles})
        #cmpddatabase.delete_one({"compound_inchi_key": "IICHURGZQPGTRD-UHFFFAOYSA-N"})
        
        # Check to make sure that the reaction SMILES is somewhat valid, this 
        # passes the reagents and products to RDKit to see if they are valid
        # chemical species. There is also a check to see if >> is within the
        # supplied reaction SMILES
        rules = ['>>' in reaction_smiles, type(reaction_smiles) == str]
        if all(rules):
            components = re.split('\.|>{2}', reaction_smiles)
            result = sum([-1 if type(rdkit_id_sanitize(element, 'smiles')) \
                          != dict else 0 for element in components])
            if result != 0:
                code_returns['print_statements'].append('Reaction SMILES contains one or more invalid species: ' + \
                            str(result))
                return ['Error', code_returns]
            else:
                code_returns['print_statements'].append('Reaction SMILES passed basic validity checks')
                reaction_entry = {}
                reaction_entry['date_added'] = current_time.strftime('%m/%d/%Y')
                reaction_entry['last_modified'] = current_time.strftime('%m/%d/%Y')
                split_smiles = reaction_smiles.split('>>')
                reagent_list = []
                for i in split_smiles[0].split('.'):
                    rdkit_cleaned = rdkit_id_sanitize(i, 'smiles')
                    reagent_list.append(rdkit_cleaned['basic_properties']['compound_smiles'])
                product_list = []
                for i in split_smiles[1].split('.'):
                    rdkit_cleaned = rdkit_id_sanitize(i, 'smiles')
                    product_list.append(rdkit_cleaned['basic_properties']['compound_smiles'])
                reaction_entry['reagents'] = reagent_list
                reaction_entry['products'] = product_list
                cleaned_reaction_smiles = '.'.join(reagent_list) + '>>' + '.'.join(product_list)
                if 'template' in incoming_request.keys():
                    reaction_entry['template'] = incoming_request['template']
                reaction_entry['reaction_smiles'] = cleaned_reaction_smiles
                reaction_entry['reaction_type'] = incoming_request['reaction_type']
                #rxndatabase.delete_many({"reaction_smiles": cleaned_reaction_smiles})
        else:
            code_returns['print_statements'].append('Something is wrong with the submitted reaction SMILES: ' + \
                        reaction_smiles)
            return ['Error', code_returns]

        # Since the reaction SMILES passed some basic validity checks we will
        # check to make sure that the incoming request data contains a valid set
        # of object/dictionary keys
        reacton_enter_schema = get_schema('add_reaction_details')
        #print(incoming_request['reaction_details'])
        try:
            validate(instance = incoming_request['reaction_details'], schema = reacton_enter_schema)
            code_returns['print_statements'].append('Reaction details entry for ' +\
                        reaction_smiles + ' complies with schema')
        except jsonschema.exceptions.ValidationError as e:
            code_returns['print_statements'].append('Reaction details entry for ' +\
                        reaction_smiles + ' does not comply with schema')
            code_returns['print_statements'].append(str(e))
            return ['Error', code_returns]
        
        # Next we will look through the AMD Databasae to see if there is a reaction stub
        file_count = myclient.amd_database.reaction_database.count_documents({
                "reaction_smiles": cleaned_reaction_smiles})
        if file_count == 0:
            code_returns['print_statements'].append('Inserting a new reaction stub for: ' \
                        + reaction_smiles)
            reaction_entry = rxndatabase.insert_one(reaction_entry)
        elif file_count == 1:
            code_returns['print_statements'].append('A reaction database document exists for: ' + \
                        reaction_smiles)
        else:
            code_returns['print_statements'].append('Not sure what to do about: ' + \
                        reaction_smiles + ' (Entries = ' + str(file_count) + ')')
            return ['Error', code_returns]
        reaction_entry = myclient.amd_database.reaction_database.find_one({
                "reaction_smiles": cleaned_reaction_smiles})
        # With the basic reaction stub made, it is now time to take the rest of
        # incoming request and make an entry for the reaction table in the document
        nested_values_to_check = [element for element in incoming_request['reaction_details'] \
               if element in ['reaction_reagents','reaction_solvent','reaction_catalyst']]
        for item in nested_values_to_check:
            for i in range(0, len(incoming_request['reaction_details'][item])):
            #for element in incoming_request['reaction_details'][item]:
                element = incoming_request['reaction_details'][item][i]
                rdkit_cleaned = rdkit_id_sanitize(element[2], 'smiles')
                if type(rdkit_cleaned) != dict:
                    code_returns['print_statements'].append('Invalid chemcial species RDKit: ' + element[0])
                    return ['Error', code_returns]
                incoming_request['reaction_details'][item][i][0] = rdkit_cleaned['basic_properties']['compound_smiles']
        if 'reaction_table' not in reaction_entry:
                reaction_entry['reaction_table'] = {}
        unique_name = incoming_request['reaction_details']['reaction_date'].replace('/','') + '_' + \
                        incoming_request['reaction_details']['reaction_location'][0]
        if unique_name not in reaction_entry['reaction_table']:
            reaction_entry['last_modified'] = current_time.strftime('%m/%d/%Y')
            reaction_entry['reaction_table'][unique_name] = incoming_request['reaction_details']
            rxndatabase.find_one_and_update({"reaction_smiles": cleaned_reaction_smiles},\
                                                  {"$set": reaction_entry})
        else:
            incoming_reaction_details = {}
            incoming_reaction_details['reaction_table'] = {}
            incoming_reaction_details['reaction_table'][unique_name] = incoming_request['reaction_details']
            update(reaction_entry, incoming_reaction_details)
            reaction_entry['last_modified'] = current_time.strftime('%m/%d/%Y')
            rxndatabase.find_one_and_update({"reaction_smiles": cleaned_reaction_smiles},\
                                                  {"$set": reaction_entry})
            code_returns['print_statements'].append('The incoming reaction entry details already exists')
        #pprint(reaction_entry)
        
        # Now we need to associate the reaction with the product entries
        reaction_hash = str(reaction_entry['_id'])
        associated_reaction_entry = {}
        associated_reaction_entry['associated_reactions'] = {}
        associated_reaction_entry['associated_reactions'][reaction_hash] = {'id':str(reaction_entry['_id']),
                                                                            'date_added': reaction_entry['date_added'],
                                                                            'reaction_smiles': reaction_entry['reaction_smiles'],
                                                                            'reaction_type': reaction_entry['reaction_type']}
        product_list = incoming_request['reaction_details']['reaction_products']
        for product in product_list:
            element_rdkit_info = rdkit_id_sanitize(product[0], 'smiles')
            associated_reaction_entry_loop = deepcopy(associated_reaction_entry)
            associated_reaction_entry_loop['associated_reactions'][reaction_hash]['product_type'] = product[3]
            if type(element_rdkit_info) is not dict:
                code_returns['print_statements'].append(['Something went wrong with your product entry: ' \
                            + str(product[0]) + ', ' + str('smiles')])
                code_returns['print_statements'].append([product_list])
                return ['Error', code_returns]
            else:
                file_count =  myclient.amd_database.compound_database.count_documents({"compound_inchi_key": \
                                                                               element_rdkit_info['compound_inchi_key']})
                if file_count == 0:
                    code_returns['print_statements'].append('Inserting a new compound stub for: ' \
                                + element_rdkit_info['compound_inchi_key'])
                    database_entry = cmpddatabase.insert_one(element_rdkit_info)
                elif file_count == 1:
                    code_returns['print_statements'].append('A database document exists for: ' + \
                                product[0] + ' -> ' + element_rdkit_info['compound_inchi_key'])
                else:
                    code_returns['print_statements'].append('Not sure what to do about: ' + \
                                product[0] + ' -> ' + element_rdkit_info['compound_inchi_key'] + \
                                ' (Entries = ' + str(file_count) + ')')
                    return ['Error', code_returns]
            database_entry = myclient.amd_database.compound_database.find_one({"compound_inchi_key": \
                                                                           element_rdkit_info['compound_inchi_key']})
            update(database_entry, associated_reaction_entry_loop)
            initial_schema = get_schema('initial_entry')
            #pprint(database_entry)
            try:
                validate(instance = database_entry, schema = initial_schema)
                code_returns['print_statements'].append('Database entry for ' + element_rdkit_info['compound_inchi_key'] \
                            + ' -> ' + product[0] + ' complies with schema')
            except jsonschema.exceptions.ValidationError:
                code_returns['print_statements'].append('Database entry for ' + product[0] + ' does not comply with schema')
                return ['Error', code_returns]
            cmpddatabase.update_one({"compound_inchi_key": element_rdkit_info['compound_inchi_key']}\
                                 ,{"$set": database_entry})
            code_returns['print_statements'].append('Updating entry for ' + element_rdkit_info['compound_inchi_key'] \
                            + ' -> ' + product[0] + ' to include associated reaction')

    elif incoming_request['reaction_type'] == 'degradation':
        code_returns['print_statements'].append('Degradation reactions still need to be implemented')
        return ['Error', code_returns]
    else:
        code_returns['print_statements'].append('Reaction Type not understood')
        ['Error', code_returns]
    #rxndatabase.delete_one({"reaction_smiles": cleaned_reaction_smiles})
    return ['Success', code_returns]


