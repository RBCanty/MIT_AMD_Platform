# -*- coding: utf-8 -*-
"""
Some code built to update database documents in the AMD Database
Created on Thu Dec 19 19:48:10 2019
Version 1
@author: Brent
"""
import pymongo
import sys
import re
from copy import deepcopy
from jsonschema import validate
from data_entry import rdkit_id_sanitize
from schema_container import get_schema
from data_entry import without_keys
from pprint import pprint
from bson.objectid import ObjectId
import traceback
from pprint import pprint

def update_data_mongo(incoming_request, mongo_database):
    # Fields that are needed: collection, document_id, previous_document, new_document
    code_returns = {}
    # Brings in the mongo collection names for the AMD Database and initializes them
    myclient = mongo_database
    amddatabase = myclient['amd_database']
    # -------------------------------------------------------------------------
    code_returns['print_statements'] = []
    
    valid_collections = ['model_information']
    if 'collection' not in incoming_request.keys() or incoming_request['collection'] not in valid_collections:
        code_returns['print_statements'].append(['Incoming request not understood please check the collection key-value'])
        return ['Error', code_returns]
    relevant_collection = amddatabase[incoming_request['collection']]
    
    try:
        db_doc_id = incoming_request['previous_document']['_id']
    except KeyError:
        code_returns['print_statements'].append(['Previous document does not contain a valid _id field'])
        return ['Error', code_returns]
    if incoming_request['previous_document']['_id'] != incoming_request['new_document']['_id']:
        code_returns['print_statements'].append(['Previous document and new document _id field do not match'])
        return ['Error', code_returns]
    filecount = relevant_collection.count_documents({'_id': ObjectId(db_doc_id)})
    if filecount == 0 or filecount > 1:
        code_returns['print_statements'].append('Bad number of instances of "%s" in the database: %s instances' % (db_doc_id, filecount))
        return ['Error', code_returns]
    relevant_collection.update_one({'_id': ObjectId(db_doc_id)}, {"$set": without_keys(incoming_request['new_document'], '_id')})
    code_returns['print_statements'].append('Updated document with id %s' % db_doc_id)
    return ['Success', code_returns]


def add_model_information(incoming_request, mongo_database):
    code_returns = {}
    myclient = mongo_database
    amddatabase = myclient['amd_database']
    relevant_collection = amddatabase['model_information']
    code_returns['print_statements'] = []
    incoming_document = deepcopy(incoming_request['incoming_document'])
    file_count = relevant_collection.count_documents({"filename": incoming_document['filename']})
    if file_count == 0:
        code_returns['print_statements'].append('Inserting a new compound stub for: ' \
                            + incoming_document['filename'])
        relevant_collection.insert_one(incoming_document)
    elif file_count == 1:
        code_returns['print_statements'].append('A database document exists for: ' + incoming_document['filename'])
        return ['Error', code_returns]
    else:
        code_returns['print_statements'].append('Bad number of instances of %s, %s' % (incoming_document['filename'], file_count))
        return ['Error', code_returns]
    return ['Success', code_returns]
    
    
def complete_data_campaign(incoming_request, mongo_database):
    code_returns = {}
    myclient = mongo_database
    amddatabase = myclient['amd_database']
    relevant_collection = amddatabase['model_information']
    code_returns['print_statements'] = []
    incoming_document = deepcopy(incoming_request)
    problems = []
    
    try:
        if 'campaign_name' not in incoming_document.keys():
            return ['Error', 'Campaign name was not provided, unable to update document']
        if 'model_details' not in incoming_document.keys():
            return ['Error', 'Model documents to update were not provided']
        model_docs_to_update = incoming_document['model_details']
        
        for model_doc_to_update in model_docs_to_update.keys():
            file_count = relevant_collection.count_documents({"filename": model_doc_to_update})
            if file_count != 1:
                code_returns['print_statements'].append('Bad number of instances of %s, %s' % (model_doc_to_update, file_count))
                return ["Error", code_returns]
            relevant_document = relevant_collection.find_one({"filename": model_doc_to_update})
            
            model_retraining_dataset = relevant_document['model_retraining_dataset']['dataset']
            if incoming_document['campaign_name'] in model_retraining_dataset.keys():
                campaign_status = model_docs_to_update[model_doc_to_update]['status']
                model_retraining_dataset[incoming_document['campaign_name']]['status'] = campaign_status
                code_returns['print_statements'].append('Marked %s complete for campaign %s' % (model_doc_to_update,
                                                                                             incoming_document['campaign_name']))
            elif relevant_document['model_location'] is None:
                campaign_status = model_docs_to_update[model_doc_to_update]['status']
                model_retraining_dataset[incoming_document['campaign_name']] = {'status': campaign_status}
                code_returns['print_statements'].append('Marked %s complete for campaign %s' % (model_doc_to_update,
                                                                                             incoming_document['campaign_name']))
            else:
                problems.append('Model has already been retrained for "%s", campaign "%s"' (model_doc_to_update,
                                                                                        incoming_document['campaign_name']))
                continue
            num_in_progress = 0
            finished_status = ['completed', 'failed']
            for campaign in model_retraining_dataset.keys():
                if model_retraining_dataset[campaign]['status'] not in finished_status:
                    num_in_progress += 1
            
            if num_in_progress == 0 and 'experiments' in  relevant_document['job_status']:
                relevant_document['job_status'] = 'experiments_completed'
                code_returns['print_statements'].append('Marked %s complete, model retraining ready' % model_doc_to_update)
            relevant_collection.update_one({"filename": model_doc_to_update}, 
                                           {"$set": relevant_document})
            code_returns['print_statements'].append('Updated document with filename %s' % model_doc_to_update)
        
        if problems:
            code_returns['Error Statements'] = problems
            return ['Error', code_returns]
    except:
        code_returns['Error Statements'] = traceback.format_exc()
        pprint(traceback.format_exc())
        return ['Error', code_returns]
    return ['Success', code_returns]
    
    
    
    
    
    
    