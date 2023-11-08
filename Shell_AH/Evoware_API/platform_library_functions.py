# -*- coding: utf-8 -*-
"""
Updated 8/26/2021

@author: Brent
"""

import pymongo, time, json, datetime
import jsonschema
from jsonschema import validate
from pprint import pprint
from copy import deepcopy
from deepdiff import DeepDiff
from bson import ObjectId


def get_schema(type_tag):
    # TODO: Add schema for well-plates and pipette tips
    if type_tag in ['solvents', 'reagents', 'wellplates']:
        schema_template = {
            'type': 'object',
            'properties': {
                'container_name': {'type': 'string'},
                'labware_type': {'type': 'string'},
                'contents': {'type': ['object', 'string']},
                'date_created': {'type': 'string'},
                'date_updated': {'type': 'string'},
                'location': {'type': 'array'}
            },
            'additionalProperties': False
        }
    elif type_tag in ['consumables']:
        schema_template = {
            'type': 'object',
            'properties': {
                'on_platform': {'type': 'object'},
                'off_platform': {'type': 'object'}
            },
            'additional_properties': False
        }
    elif type_tag in ['queue']:
        schema_template = {
            'type': 'object',
            'properties': {
                'queue_name': {'type': 'string'},
                'date_created': {'type': 'string'},
                'containers': {'type': 'object'},
                'contents': {'type': 'object'},
                'operations': {'type': 'object'}
            },
            'additional_properties': False
        }
    else:
        return 'Not a currently implemented schema %s... sorry...' % type_tag
    return schema_template


def without_keys(d, keys):
    return {k: v for k, v in d.items() if k not in keys}


def labware_to_well_names(labware, cl_dict):
    try:
        nrows = cl_dict['labware'][labware]['nrows']
        ncols = cl_dict['labware'][labware]['ncols']
        nwells = nrows * ncols
    except:
        return ['Error', 'Something is not correct in carrier_library for "%s" with nrows and ncols' % labware]
    if nwells == 1:
        return ['XX']
    well_list = []
    row, col = 1, 1
    for i in range(0, nwells):
        well_list.append(chr(ord('A') + row - 1) + str(col))
        if row == nrows:
            row = 1
            col += 1
        else:
            row += 1
    return well_list


# ----------------------------------------------------------------------------
# --------------------- Useful functions for the database --------------------
# ----------------------------------------------------------------------------


def query_collection(collection, mongo_database):
    platform_database = mongo_database['platform']
    valid_collections = ['consumables', 'solvents', 'wellplates', 'reagents', 'queue', 'historical']
    success_statements = []
    if type(collection) != str:
        success_statements.append('The "collection" field requires string entry, not %s' % type(collection))
        return ['Error', success_statements]
    if collection not in valid_collections:
        success_statements.append('Requested collection %s not valid in the database' % collection)
        return ['Error', success_statements]
    relevant_collection = platform_database[collection]
    filecount = relevant_collection.count_documents({})
    if filecount >= 1:
        documents_cursor = relevant_collection.find({})
        success_statements.append('Found "%s" documents in the "%s" collection' % (str(filecount), collection))
    else:
        success_statements.append('No documents found in collection: %s' % collection)
        return ['Error', success_statements]
    documents_to_return = {}
    for collection_document in documents_cursor:
        document_to_return = collection_document
        document_to_return['_id'] = str(document_to_return['_id'])
        documents_to_return[str(document_to_return['_id'])] = document_to_return
    return [documents_to_return, success_statements]


def query_document(collection, search_term, search_field, mongo_database):
    platform_database = mongo_database['platform']
    valid_collections = ['queue', 'consumables', 'solvents', 'wellplates', 'reagents', 'historical']
    success_statements = []
    if search_field == '_id':
        search_term = ObjectId(search_term)
        success_statements.append('Query document understood for "%s" search field' % search_field)
    if type(collection) != str:
        success_statements.append('The "collection" field requires string entry, not %s' % type(collection))
        return ['Error', success_statements]
    if collection not in valid_collections:
        success_statements.append('Requested collection %s not valid in the database' % collection)
        return ['Error', success_statements]
    success_statements.append('The requested collection "%s" is a valid database collection' % collection)
    relevant_collection = platform_database[collection]
    if search_field in ['location']:
        success_statements.append(
            'The targeted search field is valid for multi-document returns, expect a dictionary of dictionaries')
        filecount = relevant_collection.count_documents({search_field: search_term})
        success_statements.append('Found "%s" documents matching search term: %s' % (str(filecount), search_term))
        document_cursor = relevant_collection.find({search_field: search_term})
        documents_to_return = {}
        for collection_document in document_cursor:
            document_to_return = collection_document
            document_to_return['_id'] = str(document_to_return['_id'])
            documents_to_return[str(document_to_return['_id'])] = document_to_return
    elif search_field in ['labware_type']:
        success_statements.append(
            'The targeted search field is valid for multi-document returns, expect a dictionary of dictionaries')
        filecount = relevant_collection.count_documents({search_field: search_term})
        success_statements.append('Found "%s" documents matching search term: %s' % (str(filecount), search_term))
        document_cursor = relevant_collection.find({search_field: search_term})
        documents_to_return = {}
        for collection_document in document_cursor:
            document_to_return = collection_document
            document_to_return['_id'] = str(document_to_return['_id'])
            documents_to_return[str(document_to_return['_id'])] = document_to_return
    elif search_field in ['campaign_name']:
        success_statements.append(
            'The targeted search field is valid for multi-document returns, expect a dictionary of dictionaries')
        filecount = relevant_collection.count_documents({search_field: search_term})
        success_statements.append('Found "%s" documents matching search term: %s' % (str(filecount), search_term))
        document_cursor = relevant_collection.find({search_field: search_term})
        documents_to_return = {}
        for collection_document in document_cursor:
            document_to_return = collection_document
            document_to_return['_id'] = str(document_to_return['_id'])
            documents_to_return[str(document_to_return['_id'])] = document_to_return
    else:
        success_statements.append(
            'The targeted search field is valid for single-document returns, expect a single dictionary')
        filecount = relevant_collection.count_documents({search_field: search_term})
        if filecount == 0 or filecount > 1:
            success_statements.append(
                'Bad number of instances of "%s" in the database: %s instances' % (search_term, filecount))
            return ['Error', success_statements]
        else:
            documents_to_return = relevant_collection.find_one({search_field: search_term})
            documents_to_return['_id'] = str(documents_to_return['_id'])
            success_statements.append('Found "%s" documents matching search term: %s' % (str(filecount), search_term))
    return [documents_to_return, success_statements]


def move_database_document(db_doc_id, current_collection, new_collection, mongo_database):
    platform_database = mongo_database['platform']
    valid_collections = ['consumables', 'solvents', 'wellplates', 'reagents', 'queue', 'historical']
    success_statements = []
    if type(current_collection) != str or type(new_collection) != str:
        success_statements.append('The "collection" fields requires string entry')
        return ['Error', success_statements]
    if current_collection not in valid_collections or new_collection not in valid_collections:
        success_statements.append('Requested collection not valid in the database')
        return ['Error', success_statements]
    success_statements.append(
        'The requested collections "%s->%s" are valid collections' % (current_collection, new_collection))
    if current_collection == new_collection:
        success_statements.append(
            'The move "%s->%s" occurs in the same collection' % (current_collection, new_collection))
        return ['Error', success_statements]
    current_db_collection = platform_database[current_collection]
    filecount = current_db_collection.count_documents({'_id': ObjectId(db_doc_id)})
    if filecount == 0 or filecount > 1:
        success_statements.append('Bad number of instances of "%s" in the %s collection: %s instances' % (
        db_doc_id, current_collection, filecount))
        return ['Error', success_statements]
    success_statements.append(
        'Found %s instance of %s in the %s collection' % (filecount, db_doc_id, current_collection))
    document_to_move = current_db_collection.find_one({'_id': ObjectId(db_doc_id)})
    new_db_collection = platform_database[new_collection]
    filecount = new_db_collection.count_documents({'_id': ObjectId(db_doc_id)})
    if filecount > 0:
        success_statements.append('Bad number of instances of "%s" in the %s collection: %s instances' % (
        db_doc_id, new_collection, filecount))
        return ['Error', success_statements]
    success_statements.append('Found %s instances of %s in the %s collection' % (filecount, db_doc_id, new_collection))
    new_db_collection.insert_one(document_to_move)
    current_db_collection.delete_one({'_id': ObjectId(db_doc_id)})
    new_db_document = new_db_collection.find_one({'_id': ObjectId(db_doc_id)})
    new_db_document['_id'] = str(new_db_document['_id'])
    success_statements.append('Moved document %s from %s to %s' % (db_doc_id, current_collection, new_collection))
    return ['Success', success_statements]


def update_database_document(previous_document, new_document, collection, mongo_database):
    platform_database = mongo_database['platform']
    valid_collections = ['consumables', 'solvents', 'wellplates', 'reagents', 'queue']
    success_statements = []
    if type(collection) != str:
        success_statements.append('The "collection" field requires string entry, not %s' % type(collection))
        return ['Error', success_statements]
    if collection not in valid_collections:
        success_statements.append('Requested collection "%s" not valid in the database' % collection)
        return ['Error', success_statements]
    success_statements.append('The requested collection "%s" is a valid database collection' % collection)
    relevant_collection = platform_database[collection]
    try:
        db_doc_id = previous_document['_id']
    except KeyError:
        success_statements.append('Previous document does not contain a valid _id field')
        return ['Error', success_statements]
    if previous_document['_id'] != new_document['_id']:
        success_statements.append('The two documents do not have consistent id values')
        return ['Error', success_statements]
    document_changes = DeepDiff(new_document, previous_document, ignore_order=True)
    if document_changes == False:
        success_statements.append('No valid changes to make to the document')
        return ['Error', success_statements]
    filecount = relevant_collection.count_documents({'_id': ObjectId(db_doc_id)})
    if filecount == 0 or filecount > 1:
        success_statements.append(
            'Bad number of instances of "%s" in the database: %s instances' % (new_document["_id"], filecount))
        return ['Error', success_statements]
    relevant_collection.update_one({'_id': ObjectId(db_doc_id)}, {"$set": without_keys(new_document, '_id')})
    success_statements.append('Updated document: %s' % db_doc_id)
    return ['Success', success_statements]


def query_platform_location(loc_tag, sub_location, mongo_database):
    platform_database = mongo_database['platform']
    valid_collections = ['consumables', 'solvents', 'wellplates', 'reagents']
    success_statements = []
    gathered_documents = {'partial_matches': {}, 'full_matches': {}}
    for collection in valid_collections:
        current_db_collection = platform_database[collection]
        file_cursor = current_db_collection.find({'location': {"$elemMatch": {"$in": [sub_location.split(',')]}}})
        for file in file_cursor:
            if file['location'] == [loc_tag, sub_location]:
                gathered_documents['full_matches'][str(file['_id'])] = file
                gathered_documents['full_matches'][str(file['_id'])]['_id'] = str(
                    gathered_documents['full_matches'][str(file['_id'])]['_id'])
            else:
                gathered_documents['partial_matches'][str(file['_id'])] = file
                gathered_documents['partial_matches'][str(file['_id'])]['_id'] = str(
                    gathered_documents['partial_matches'][str(file['_id'])]['_id'])
    success_statements.append('Gathered things found in consumables, solvents, wellplates and reagents')
    return [gathered_documents, success_statements]


def add_to_database(incoming_dict, collection, carriers_labware_dict, mongo_database):
    platform_database = mongo_database['platform']
    valid_collections = ['consumables', 'solvents', 'wellplates', 'reagents', 'queue']
    success_statements = []
    if type(collection) != str:
        success_statements.append('The "collection" field requires string entry, not %s' % type(collection))
        return ['Error', success_statements]
    if collection not in valid_collections:
        success_statements.append('Requested collection "%s" not valid in the database' % collection)
        return ['Error', success_statements]
    relevant_collection = platform_database[collection]
    success_statements.append('The requested collection "%s" is a valid database collection' % collection)

    # We check to see if the entry is consistent with expected schema
    relevant_schema = get_schema(collection)
    if type(relevant_schema) != dict:
        success_statements.append('Issue getting the schema for "%s", check entries' % collection)
        return ['Error', success_statements]
    try:
        schema_check = validate(instance=incoming_dict, schema=relevant_schema)
    except jsonschema.exceptions.ValidationError as e:
        return ['Error', e]

    # Now that the code meets some basic type checks from the jsonschema we then check some
    # of the specifics about the labware and carriers and such (very much a WIP at this point)
    if collection in ['solvents', 'wellplates', 'reagents']:
        try:
            labware_type = incoming_dict['labware_type']
            lib_entry = carriers_labware_dict['labware'][labware_type]
            nrows = carriers_labware_dict['labware'][labware_type]['nrows']
            ncols = carriers_labware_dict['labware'][labware_type]['ncols']
            nwells = nrows * ncols
        except:
            success_statements.append(
                'Something is not correct in carrier_library for "%s"' % incoming_dict['labware_type'])
            return ['Error', success_statements]
        labware_return = labware_to_well_names(labware_type, carriers_labware_dict)
        if labware_return[0] == 'Error':
            success_statements.append(labware_return[1])
            return ['Error', success_statements]
        # Now we can check to make sure that the well names are consistent with the expected names
        problem_wells = []
        if incoming_dict['contents'] == 'Empty' and collection == 'wellplates':
            success_statements.append('Preparing to insert a new empty wellplate into the wellplates collection')
            pass
        else:
            for key in incoming_dict['contents'].keys():
                if key not in labware_return:
                    problem_wells.append(key)
            if len(problem_wells) != 0:
                success_statements.append(
                    'Unexpected well names for "%s": %s' % (labware_type, ','.join(problem_wells)))
                return ['Error', success_statements]
        filecount = relevant_collection.count_documents({"container_name": incoming_dict['container_name']})
        if filecount == 0:
            success_statements.append(['Inserting document with title "%s"' % incoming_dict['container_name']])
            relevant_collection.insert_one(incoming_dict)
        else:
            success_statements.append(
                'Document already exists with container name "%s"' % incoming_dict['container_name'])
            return ['Error', success_statements]
        filecount = relevant_collection.count_documents({"container_name": incoming_dict['container_name']})
    elif collection == 'consumables':
        filecount = relevant_collection.count_documents({"container_name": incoming_dict['container_name']})
        if filecount == 0:
            success_statements.append(['Inserting document into the consumables collection'])
            relevant_collection.insert_one(incoming_dict)
        else:
            success_statements.append(
                'Document already exists for the consumable container_name "%s"' % incoming_dict['container_name'])
            return ['Error', success_statements]
    elif collection == 'queue':
        filecount = relevant_collection.count_documents({"queue_name": incoming_dict["queue_name"]})
        if filecount == 0:
            success_statements.append('Inserting document "%s" into the queue collection' % incoming_dict['queue_name'])
            relevant_collection.insert_one(incoming_dict)
        else:
            success_statements.append(
                'Document "%s" already exists in the queue collection' % incoming_dict['queue_name'])
            return ['Error', success_statements]
    else:
        success_statements.append('Collection type "%s" has not been implemented yet' % collection)
        return ['Error', success_statements]
    success_statements.append('Added document to the %s collection' % collection)
    return ['Success', success_statements]


def update_database_field(db_plate_name, update_field, old_value, new_value, mongo_database):
    platform_database = mongo_database['platform']
    valid_collections = ['consumables', 'solvents', 'wellplates', 'reagents']
    success_statements = []
    gathered_documents = []
    for collection in valid_collections:
        current_db_collection = platform_database[collection]
        file_cursor = current_db_collection.find({'container_name': db_plate_name})
        for file in file_cursor:
            gathered_documents.append([collection, file])
    current_db_collection = platform_database['queue']
    file_cursor = current_db_collection.find({'queue_name': db_plate_name})
    for file in file_cursor:
        gathered_documents.append(['queue', file])
    if len(gathered_documents) != 1:
        success_statements.append(
            'Bad number of instances of "%s" in the platform library: %s instances: %s' % (db_plate_name,
                                                                                           len(gathered_documents),
                                                                                           ' > '.join(
                                                                                               [str(file['_id']) for
                                                                                                file in
                                                                                                gathered_documents])))
        return ['Error', success_statements]
    else:
        success_statements.append('Found %s document instance in the database' % len(gathered_documents))
    relevant_document = gathered_documents[0][1]
    document_collection = gathered_documents[0][0]
    nested_key_list = update_field.split('.')
    relevant_field = relevant_document
    try:
        for nested_key in nested_key_list:
            relevant_field = relevant_field[nested_key]
    except Exception as e:
        success_statements.append(str(e))
        success_statements.append('Issue with the nested key lookup')
        return ['Error', success_statements]

    if relevant_field != old_value:
        success_statements.append(
            'Inconsistency between database value (%s) and sent old value (%s)' % (relevant_field, old_value))
        return ['Error', success_statements]

    update_collection = platform_database[document_collection]
    return_value = update_collection.update_one({"_id": relevant_document['_id']}, {"$set": {update_field: new_value}})
    success_statements.append(
        'Updated the old value "%s" in document "%s" to new value "%s"' % (old_value, db_plate_name, new_value))
    return ['Success', success_statements]


# Below is not likely to be used at the moment... it requires index searching...
def query_database(chemical_name, collection, mongo_database):
    platform_database = mongo_database['platform']
    valid_collections = ['consumables', 'solvents', 'wellplates', 'queue']
    success_statements = []
    if type(collection) != str:
        return ['Error', 'The "collection" field requires string entry, not %s' % type(collection)]
    if collection not in valid_collections:
        return ['Error', 'Requested collection not valid in the database']
    success_statements.append('Request collection %s is valid' % collection)
    relevant_collection = platform_database[collection]
    relevant_collection.create_index([('$**', pymongo.TEXT)],
                                     name='search_index', default_language='english')
    filecount = relevant_collection.count_documents({"$text": {"$search": chemical_name}})
    if filecount >= 1:
        documents_cursor = relevant_collection.find({"$text": {"$search": chemical_name}})
        success_statements.append('Found "%s" documents matching search term: %s' % (str(filecount), chemical_name))
    else:
        return ['Error', 'No documents found matching search term: %s' % chemical_name]
    documents_to_return = {}
    for item in documents_cursor:
        documents_to_return[str(item['_id'])] = item
    return [documents_to_return, success_statements]


def update_database(previous_collection, new_collection, collection, mongo_database):
    # Take all of the files from the collection and then update each of the documents in the database (!)
    # This way any change only happens at the local level before being pushed to the database (!)
    platform_database = mongo_database['platform']
    valid_collections = ['consumables', 'solvents', 'wellplates', 'reagents', 'queue']
    success_statements = []
    if type(collection) != str:
        return ['Error', 'The "collection" field requires string entry, not %s' % type(collection)]
    if collection not in valid_collections:
        return ['Error', 'Requested collection not valid in the database']
    relevant_collection = platform_database[collection]
    for new_doc_key in new_collection.keys():
        new_document = new_collection[new_doc_key]
        for old_doc_key in previous_collection.keys():
            if new_doc_key == old_doc_key:
                old_document = previous_collection[old_doc_key]
                break
        else:
            return ['Error', 'Document "%s" was not in the previous_collection' % new_document["_id"]]
        document_changes = DeepDiff(new_document, old_document, ignore_order=True)
        if document_changes == False:
            continue
        filecount = relevant_collection.count_documents({'_id': ObjectId(new_document["_id"])})
        if filecount == 0 or filecount > 1:
            return ['Error',
                    'Bad number of instances of "%s" in the database: %s instances' % (ObjectId(new_document["_id"]), filecount)]
        relevant_collection.update_one({'_id': ObjectId(new_document["_id"])}, {"$set": without_keys(new_document, '_id')})
        success_statements.append('Updated document: %s' % str(ObjectId(new_document["_id"])))
        document_cursor = relevant_collection.find({})
        # for item in document_cursor:
        #    pprint(item)
    return ['', success_statements]
