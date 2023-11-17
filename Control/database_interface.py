""" database_interface.py

Used to communicate with the Database for task definition

@author: Ben C, based on DB interface schemas developed by Brent K
"""
import datetime
import json
from copy import deepcopy
from io import StringIO
from pprint import pformat
from typing import List, Union

import requests
import yaml

import custom_exceptions as cexc
import functionals
import mcn_status as mcs
from custom_classes import Narg, safe_open, pathfinder
from database_constants import *

DATABASE_URL, DATABASE_PORT = mcs.MCN_CFG[mcs.S_DATABASES]['Platform']


def update_database(url, port):
    """ Changes the URL and PORT values for the interface

    :param url: The new url for the database (e.g. 'http://01.23.45.678:123453')
    :param port: The new port for the database (e.g. '678907')
    :return: None
    """
    if (not url) or (not port):
        return
    global DATABASE_URL, DATABASE_PORT
    try:
        DATABASE_URL = str(url)
        DATABASE_PORT = str(int(port))
    except (TypeError, ValueError):
        DATABASE_URL, DATABASE_PORT = mcs.MCN_CFG[mcs.S_DATABASES]['Platform']


def get_available_queue_names():
    """ Returns a list of the names of all queues in the Queues collection

    :return: list(generator(queue_names))
    """

    document, err_details, resp_code = query_collection('MC', COL_QUEUE)

    def _grab(doc):
        for doc_id in doc.keys():
            yield doc[doc_id][Q_NAME]

    try:
        return list(_grab(document))
    except KeyError as ke:
        raise cexc.BadQueueFormatError(*ke.args)


def get_historical_queue_names():
    """ Returns a list of the names of all queues in the Historical collection

    :return: list(generator(queue_names))
    """
    document, _, _ = query_collection('MC', COL_HISTORICAL)

    def _grab(doc):
        for doc_id in doc.keys():
            name = doc[doc_id].get(Q_NAME, None)
            if name:
                yield name

    try:
        return list(_grab(document))
    except KeyError as ke:
        raise cexc.BadQueueFormatError(*ke.args)


def database_request(name, request: dict, get_details=False):
    """ Convenience function for making database requests

    Supports the following 'request_type's:

    query_document
      * collection --- acceptable strings: reagents, wellplates, queue, consumables, solvents, historical
      * search_term --- string input for the search
      * search_field --- requires understanding the structure of the document
    query_collection
      * collection --- acceptable strings: reagents, wellplates, queue, consumables, solvents, historical
    move_database_document
      * doc_id --- document ID hash string from Mongo
      * current_collection --- acceptable strings: reagents, wellplates, queue, consumables, solvents, historical
      * new_collection --- acceptable strings: reagents, wellplates, queue, consumables, solvents, historical
    update_database_document
      * previous_document --- returned Mongo database document
      * new_document --- updated document for the Mongo database
      * collection --- acceptable strings: reagents, wellplates, queue, consumables, solvents, historical
    add_to_database
      * incoming_dict --- this needs to be formatted consistent with other MongoDB documents
      * collection --- acceptable strings: reagents, wellplates, queue, consumables, solvents, historical
    update_database_field
      * document_name --- string that matches container_name or queue_name
      * update_field --- string using dot notation for nested keys
      * old_value --- exact match to the expected database document value
      * new_value --- new value matching the expected field formatting
    query_platform_location
      * location_tag --- the main location (location: [location_tag, sublocation])
      * sublocation --- the location (matching format) relative to main (location: [location_tag, sublocation])

    :param name: Name of system (or subsystem) loggin in
    :param request: A request dictionary
    :param get_details: True-Returns request_return[0], request_return[1], resp.status_code, False-Raise exception
    :return: The decoded response, the error details, and the response status code
    """
    incoming_request = {'auth': mongo_auth(name),
                        'request': request}

    try:
        resp = requests.post(DATABASE_URL, json=incoming_request, headers={'content-type': 'application/authjson'})
    except requests.exceptions.ConnectionError:
        if get_details:
            return "Connection Failed", "", ""
        else:
            raise cexc.DatabaseRequestError(f"Failed to connect to database, check connection")

    request_return = json.loads(resp.content.decode('utf-8'))['request_return']

    resp_code = resp.status_code
    try:
        req_doc, err_det = request_return[0], request_return[1]
    except IndexError:
        try:
            req_doc, err_det = request_return[0], "N.A."
        except IndexError:
            req_doc, err_det = "N.A.", "N.A."

    if (resp_code != 200) or (req_doc == "Error"):
        if get_details:
            return request_return[0], request_return[1], resp.status_code
        else:
            # This bypass should be moved to the query_collection() method as something that is
            #   handled by a try-catch rather that it bleeding into this method.
            # In short, this message is unique to the request being a query_collection on the queue collection
            #   and this says, "instead of raising the exception, just return that the queue collection is empty"
            # Is query_collection() raising an exception when a collection (other than queue) is empty required/expected
            #   by any existing code?
            if 'No documents found in collection: queue' in str(err_det):
                return {}, request_return[1], resp.status_code
            raise cexc.DatabaseRequestError(f"Request {request} by {name} failed:\n"
                                            f'----Response Object: Start----\n'
                                            f'{pformat(req_doc)}\n\n'
                                            f'{pformat(err_det)}\n'
                                            f'----Response Object: End----')

    return request_return[0], request_return[1], resp.status_code


def query_document(name, collection, search_field, search_term, get_details=False):
    """ Requests a document from Database

    :raises DatabaseRequestError: When DB socket fails

    :param name: Name of requester (for log-in)
    :param collection: A DB collection (reagents, wellplates, queue, consumables, solvents, historical)
    :param search_field: Field being searched (use dot to nest, requires understanding the structure of the document)
    :param search_term: string input for the search
    :param get_details:
    :return: The decoded response, the error details, and the response status code
    """
    my_request = {'request_type': 'query_document',
                  'collection': collection,
                  'search_field': search_field,
                  'search_term': search_term}
    return database_request(name, my_request, get_details)


def query_collection(name, collection):
    """ Requests a collection from Database

    :raises DatabaseRequestError: When DB socket fails

    :param name: Name of requester (for log-in)
    :param collection: A DB collection (reagents, wellplates, queue, consumables, solvents, historical)
    :return: The decoded response, the error details, and the response status code
    """
    my_request = {'request_type': 'query_collection',
                  'collection': collection}
    return database_request(name, my_request)

    # document, error_details, request_code = database_request(name, my_request, get_details=True)
    #
    # if (request_code != 200) or (document == "Error"):
    #     if 'No documents found in collection:' in str(error_details):
    #         return dict(), error_details, request_code
    #     raise cexc.DatabaseRequestError(f"Request {my_request} by {name} failed:\n"
    #                                     f'----Response Object: Start----\n'
    #                                     f'{pformat(document)}\n\n'
    #                                     f'{pformat(error_details)}\n'
    #                                     f'----Response Object: End----')
    # return document, error_details, request_code


def move_database_document(name, doc_id, current_location, new_location):
    """ Moves a document within the Database

    :raises DatabaseRequestError: When DB socket fails

    :param name: Name of requester (for log-in)
    :param doc_id: Document ID hash string from Mongo
    :param current_location: A DB collection (reagents, wellplates, queue, consumables, solvents, historical)
    :param new_location: A DB collection (reagents, wellplates, queue, consumables, solvents, historical)
    :return: The decoded response, the error details, and the response status code
    """
    my_request = {'request_type': 'move_database_document',
                  'doc_id': doc_id,
                  'current_collection': current_location,
                  'new_collection': new_location}
    return database_request(name, my_request)


def update_document(name, collection, previous_document, new_document):
    """ Updates a Database document

    :raises DatabaseRequestError: When DB socket fails

    :param name: Name of requester (for log-in)
    :param collection: A DB collection (reagents, wellplates, queue, consumables, solvents, historical)
    :param previous_document: The previous version of the document (used for error checking)
    :param new_document: The edited version of the document being pushed to the DB
    :return: The decoded response, the error details, and the response status code
    """
    my_request = {'request_type': 'update_database_document',
                  'previous_document': previous_document,
                  'new_document': new_document,
                  'collection': collection}
    return database_request(name, my_request)


def add_document(name, collection, incoming_dict):
    """ Adds a document to the Database

    :raises DatabaseRequestError: When DB socket fails

    :param name: Name of requester (for log-in)
    :param collection: A DB collection (reagents, wellplates, queue, consumables, solvents, historical)
    :param incoming_dict: A document matching the format of the specified collection to be added to the DB
    :return: The decoded response, the error details, and the response status code
    """
    my_request = {'request_type': 'add_to_database',
                  'incoming_dict': incoming_dict,
                  'collection': collection}
    return database_request(name, my_request)


def update_field(name, field, document_name, old_value, new_value):
    """ Edit documents by field rather than all at once

    :raises DatabaseRequestError: When DB socket fails

    :param name: Name of requester (for log-in)
    :param field: Field being updated (use dot to nest)
    :param document_name: Name of the document being edited
    :param old_value: The previous value (used for error checking)
    :param new_value: The new value
    :return: The decoded response, the error details, and the response status code
    """
    my_request = {'request_type': 'update_database_field',
                  'document_name': document_name,
                  'update_field': field,
                  'old_value': old_value,
                  'new_value': new_value}
    return database_request(name, my_request)


def query_location(name, location_tag, sublocation=None):
    """ Search by location request on Database

    :raises DatabaseRequestError: When DB socket fails

    :param name: Name of requester (for log-in)
    :param location_tag: Location being searched
    :param sublocation: A sublocation of the specified location (location: [location_tag, sublocation])
    :return: The decoded response, the error details, and the response status code
    """
    sublocation = sublocation if sublocation else location_tag
    my_request = {'request_type': 'query_platform_location',
                  'location_tag': location_tag,
                  'sublocation': sublocation}
    return database_request(name, my_request)


def update_mongo_queue(name, queuename, item_num, new_value):
    """ Changes the status field for a queued collection's operations' status

    :param name: System logging in to make edit
    :param queuename: Name of the queue document from thr DB (None will auto-succeed)
    :param item_num: The operation number (DBQ_STEP_NUM)
    :param new_value: Status of operation (QOP_COMPLETED) [yes/no]
    :return: The decoded response, the error details, and the response status code
    """
    if queuename is None:
        return

    if new_value not in [DB_YES, DB_NO]:
        raise ValueError(f"Value of QOP_COMPLETED must be '{DB_YES}' or '{DB_NO}' not {new_value}")

    old_value = DB_NO if new_value == DB_YES else DB_YES
    item_num = str(int(item_num))

    return update_field(name, f'operations.{item_num}.completed', queuename, old_value, new_value)


def mongo_auth(name):
    """ Utility function for making the authentication field for Mongo access

    Converts subsystem names to system names (only systems have username-password pairs) and makes strings lowercase to
    match the username-password pairs stored on the database

    :param name: System or subsystem logging in
    :return: A dictionary which can be used "incoming_request = {'auth': mongo_auth(name), 'request': ..."
    """
    name = mcs.rs_dir(name).lower()
    return {'port': DATABASE_PORT, 'user': f'xxxxxxx', 'password': f'xxxxxxxxx'}


def move_to_historical(q_name):
    """ Utility method for moving a queue from Queue to Historical

    :raises DatabaseRequestError: When DB socket fails

    :param q_name: Name of queue document (for 'search_term': q_name)
    :return: move_database_document() or Raises DatabaseRequestError
    """
    decoded_content, _, _ = query_document('MC', COL_QUEUE, Q_NAME, q_name)

    doc_id = decoded_content[DBG_ID]

    return move_database_document('MC', doc_id, COL_QUEUE, COL_HISTORICAL)


def check_queue_step(q_name: str, step_number: str):
    """ Utility function to check is a given step is completed or not

    :raises DatabaseRequestError: When DB socket fails
    :raises BadQueueFormatError: When the Queue document appears incorrectly formatted

    :param q_name: Name of queue document
    :param step_number: Step # being inquired
    :return: Value of completion from DB or raises DatabaseRequestError/BadQueueFormatError
    """

    db_q_document, err_details, resp_code = query_document('MC', COL_QUEUE, Q_NAME, q_name)

    try:
        operation_list = db_q_document[Q_OPERATIONS_LIST]
    except KeyError:
        raise cexc.BadQueueFormatError(f"Queue {q_name} does not have an '{Q_OPERATIONS_LIST}' field")

    try:
        operation = operation_list[step_number]
    except KeyError:
        raise cexc.BadQueueFormatError(f"Queue {q_name} does not have a step #{step_number}")

    try:
        return operation[QOP_COMPLETED]
    except KeyError:
        raise cexc.BadQueueFormatError(f"Queue {q_name} step #{step_number} does not have a '{QOP_COMPLETED}' field")


def insert_queue_steps(queue_name, index, new_steps: List[dict]):
    """ Inserts a set of ordered steps into a queue document

    Example: given a queue document {'1': step1, '2': step2} if (index, new_steps) = (1, [stepN, stemp M]) then
    the resulting queue document will be {'1': step1, '2': stepN, '3': stepM, '4': step2}

    :param queue_name: Name of the queue being edited
    :param index: (int/str) the index after which new_steps should be inserted
    :param new_steps: A list of new steps
    :return: None
    """
    index = int(index)
    original_queue_document, _, _ = query_document('MC', COL_QUEUE, Q_NAME, queue_name)
    updated_queue_document = deepcopy(original_queue_document)

    sorted_step_nums = sorted([int(step_number) for step_number in updated_queue_document[Q_OPERATIONS_LIST].keys()])
    sorted_steps = [updated_queue_document[Q_OPERATIONS_LIST][str(s)] for s in sorted_step_nums]
    pre = sorted_steps[:index]
    post = sorted_steps[index:]
    new_steps = pre + new_steps + post

    del updated_queue_document[Q_OPERATIONS_LIST]
    updated_queue_document[Q_OPERATIONS_LIST] = {str(i + 1): v for i, v in enumerate(new_steps)}

    update_document('MC', 'queue', original_queue_document, updated_queue_document)


def move(item: str, *, _from: str, _to: str):
    """ Changes a queue status value

    Name is legacy (back when queues were stored locally and moved between lists)

    :param item: Name of queue document
    :param _from: A queue status: ['idle', 'in_progress', 'error', 'complete']
    :param _to: A queue status: ['idle', 'in_progress', 'error', 'complete']
    :return: update_field("MC", "status", item, _from, _to)
    """
    if _from == "Lookup":
        doc, *_ = query_document('MC', COL_QUEUE, Q_NAME, item)
        _from = doc[DBQ_STATUS]

    valid_stata = [DBQS_IDLE, DBQS_WORK, DBQS_FAIL, DBQS_DONE]
    if _from not in valid_stata:
        raise ValueError(f"move.arg '_from' must be in {valid_stata} not '{_from}'")
    if _to not in valid_stata:
        raise ValueError(f"move.arg '_to' must be in {valid_stata} not '{_to}'")

    return update_field("MC", DBQ_STATUS, item, _from, _to)


def mark_time(queue_name, step_number, operation=None, *, start=False, stop=False, reset=False, logger=None):
    """ Used to mark times for paired operations

    :param queue_name: Name of queue documents
    :param step_number: Step of the paired task
    :param operation: Name of the operation (to find the task if the step number has changed due to insertion)
    :param start: Marks the start time (clears old end-time)
    :param stop: Marks the end time
    :param reset: Resets start and end time to None
    :param logger: A logger object
    :return:
    """
    if 1*start + 1*stop + 1*reset != 1:
        raise ValueError(f"Must specify exactly one mode for mark_time")

    try:
        step_number = str(int(step_number))
    except ValueError as ve:
        ve.args = tuple([f"Queue step number must be a number not {step_number} ({type(step_number)})", ]
                        + list(ve.args))
        raise
    _step_number = step_number

    try:
        doc, _, _ = query_document('MC', COL_QUEUE, Q_NAME, queue_name)
    except cexc.DatabaseRequestError:
        doc, _, _ = query_document('MC', COL_HISTORICAL, Q_NAME, queue_name)

    if operation is not None:
        step_number = functionals.find_shifted_step_number(doc[Q_OPERATIONS_LIST], step_number, operation, queue_name)

    if logger:
        preamble = "Marking Start"*start + "Marking Stop"*stop + "Resetting"*reset
        addendum = f" ({_step_number}->{step_number})"*(step_number != _step_number)
        logger.info(f"{preamble} time for {queue_name}:{step_number}{addendum}")

    try:
        if start:
            old_value = doc[Q_OPERATIONS_LIST][step_number][QOP_START]
            time_stamp = datetime.datetime.now().strftime(mcs.TIME_FORMAT)
            update_field('MC', f'{Q_OPERATIONS_LIST}.{step_number}.{QOP_START}', queue_name, old_value, time_stamp)
            # And wipe stop_time
            old_value = doc[Q_OPERATIONS_LIST][step_number][QOP_END]
            if old_value:
                update_field('MC', f'{Q_OPERATIONS_LIST}.{step_number}.{QOP_END }', queue_name, old_value, None)
        elif stop:
            old_value = doc[Q_OPERATIONS_LIST][step_number][QOP_END]
            time_stamp = datetime.datetime.now().strftime(mcs.TIME_FORMAT)
            update_field('MC', f'{Q_OPERATIONS_LIST}.{step_number}.{QOP_END}', queue_name, old_value, time_stamp)
        elif reset:
            old_value = doc[Q_OPERATIONS_LIST][step_number][QOP_START]
            if old_value:
                update_field('MC', f'{Q_OPERATIONS_LIST}.{step_number}.{QOP_START}', queue_name, old_value, None)
            old_value = doc[Q_OPERATIONS_LIST][step_number][QOP_END]
            if old_value:
                update_field('MC', f'{Q_OPERATIONS_LIST}.{step_number}.{QOP_END}', queue_name, old_value, None)
    except KeyError as ke:
        raise cexc.BadQueueFormatError(f"Queue '{queue_name}' step #{step_number} missing key '{ke.args[0]}'")


def get_pairedness(task):
    """ Gets the is_paired value

    :param task: A step in an operations list
    :return: the value of is_paired or 'no' if details or is_paired is missing
    """
    return task.get(QOP_DETAILS, dict()).get(QDET_PAIRED, 'no')


def get_reasonable_time(task, default=1, minimum=1):
    """ Gets filtered version of a time estimate

    :param task: A step in an operations list
    :param default: The default time
    :param minimum: The minimum time
    :return: max(task.get(QOP_TIME, default), minimum)
    """
    return max(task.get(QOP_TIME, default), minimum)


def get_time_est(queue_name, step, default):
    """ Retrieves the 'time_est' value from a queue's operation list

    :param queue_name: The name of a queue document
    :param step: The step number
    :param default: Return if database access fails
    :return: Pulls the time estimate from a queue document
    """
    try:
        doc, _, _ = query_document('MC', COL_QUEUE, Q_NAME, queue_name)
        return int(doc[Q_OPERATIONS_LIST][str(step)][QOP_TIME])
    except (cexc.DatabaseRequestError, KeyError, ValueError, TypeError):
        return int(default)


def set_time_est(queue_name, step, new_time_est, _raise=False):
    """ Updates the 'time_est' value in a queue's operation list

    :param queue_name: The name of a queue document
    :param step: The step number
    :param new_time_est: New value for the time estimate
    :param _raise: True - raise exceptions, False - ignore them
    :return: None
    """
    try:
        doc, _, _ = query_document('MC', COL_QUEUE, Q_NAME, queue_name)
        old_time_est = doc[Q_OPERATIONS_LIST][str(step)][QOP_TIME]
        update_field('MC', f"{Q_OPERATIONS_LIST}.{step}.{QOP_TIME}", queue_name, old_time_est, new_time_est)
    except Exception as e:
        if _raise:
            raise e


def get_plate_details(plate_name, user_name="MC", location_restriction=None, include_historical=False):
    """ Looks up a plate/container by name ('container_name') in multiple collections

    Searches, in order: wellplates, consumables, reagents, solvent, (historical - if specified)

    :param plate_name: Name of the plate/container
    :param user_name: (optional, for DB logon purposes)
    :param location_restriction: Only return if the location matches (matches location[0])
    :param include_historical: Adds historical to searched locations
    :return: (database document: dict, name of the collection: str)
    """
    collections = [COL_WELLPLATES, COL_CONSUMABLES, COL_REAGENTS, COL_SOLVENTS]
    if include_historical:
        collections.append(COL_HISTORICAL)
    # hit = {DBG_ID: None,
    #        DBG_CONTAINER_NAME: None,
    #        DBG_DATE_UPDATED: None,
    #        DBG_LABWARE_TYPE: None,
    #        DBG_LOCATION: None,
    #        W_CONT_CATEGORY: None,
    #        DBG_CONTENTS: None,
    #        DBG_DATE_CREATED: None,
    #        W_BARCODE: None,
    #        CON_NAME: None,
    #        CON_NUM: None,
    #        CON_TIP_LOCATIONS: None,
    #        CON_STATUS: None}
    hit = dict()
    hit_collection = None

    for collection in collections:
        try:
            doc, _, _ = query_document(user_name, collection, DBG_CONTAINER_NAME, plate_name)
        except cexc.DatabaseGeneralException:
            pass
        else:
            hit_collection = collection
            hit.update(doc)
            break
    if location_restriction is None:
        return hit, hit_collection
    location = hit.get(DBG_LOCATION, [None, ])[0]
    if location == location_restriction:
        return hit, hit_collection
    else:
        return None, hit_collection


def get_plates_by_type(_type, user_name="MC", location_restriction=None, include_historical=False,
                       sort=False, must_be_empty=True):
    """ Looks up containers by a plate_type

    Searches: wellplates, consumables, reagents, solvent, (historical - if specified)

    :param _type: The type ('labware_type') of the plate/container
    :param user_name: (Optional) for DB logon
    :param location_restriction: Filter out search hits if the location[0] field does not match the restriction
    :param include_historical: Includes historical in collections being searched
    :param sort: Will sort results by their location field (missing/type-mismatch values will be sorted as 0/'0')
    The sorting hierachy is location[0] by str value, then location[1][0, then 1, then 2] by int high-low
    :param must_be_empty: Filters out plates which are not empty (contents: {})
    :return: A list of DB documents
    """
    collections = [COL_WELLPLATES, COL_CONSUMABLES, COL_REAGENTS, COL_SOLVENTS]
    if include_historical:
        collections.append(COL_HISTORICAL)
    hits = list()
    for collection in collections:
        try:
            doc, _, _ = query_collection(user_name, collection)
        except cexc.DatabaseGeneralException:
            pass
        else:
            for _, item in doc.items():
                if item.get(DBG_LABWARE_TYPE, None) != _type:
                    continue

                if location_restriction and (item.get(DBG_LOCATION, [None, ])[0] != location_restriction):
                    continue

                if must_be_empty and item.get(DBG_CONTENTS, {}):
                    continue

                hits.append(item)
    if sort:
        hits.sort(key=lambda x: _loc_sort(x, 3, int), reverse=True)
        hits.sort(key=lambda x: _loc_sort(x, 2, int), reverse=True)
        hits.sort(key=lambda x: _loc_sort(x, 1, int), reverse=True)
        hits.sort(key=lambda x: _loc_sort(x, 0, str), reverse=False)
        # sorted by location[0] then location[1][0], then location[1][1], etc.
    return hits


def record_fault(fault: mcs.Fault, file_locator=None):
    """ Saves a copy of a Fault to the queue document

    :param fault: The Fault
    :param file_locator: The address of the log file associated (optional)
    :return: An error message, An Exception object, or None upon success
    """
    queue_name = fault.queue
    if queue_name is None:
        return "No queue name associated"
    parent_system = mcs.rs_dir(fault.location)
    try:
        doc, _, _ = query_document(parent_system, COL_QUEUE, Q_NAME, queue_name)
    except cexc.DatabaseRequestError as dre:
        return dre
    if not isinstance(doc, dict):
        return "DB did not return a dictionary when queue was requested"
    original_doc = deepcopy(doc)
    new_doc: dict = doc
    new_doc.setdefault(Q_RECORD, list())
    existing_faults = [mcs.Fault.loadd(fault_dict) for fault_dict in new_doc[Q_RECORD]]
    if fault in existing_faults:
        return "DB already has fault"
    else:
        q_fault = fault.dumpd()
        try:
            q_fault.update({'step': get_step_number(original_doc)})
        except KeyError:
            pass
        q_fault.pop('queue', None)
        if file_locator:
            try:
                q_fault['file_locator'] = file_locator[0]
            except KeyError:
                pass
        new_doc[Q_RECORD].append(q_fault)
        new_doc[Q_RECORD].sort(key=lambda md: datetime.datetime.strptime(md['timestamp'], mcs.TIME_FORMAT))
    try:
        update_document(parent_system, COL_QUEUE, original_doc, new_doc)
    except cexc.DatabaseRequestError as dre:
        return dre
    return


def catalogue_fault(fault: mcs.Fault, _return=None):
    """ Logs a fault to its very own log file with additional diagnostic data

    :param fault: the Fault
    :param _return: A Pass-By-Reference object (such as a list; must support append(str)) to which the file's name is
      passed
    :return: An error message, An Exception object, or None upon success
    """
    # Get the Queue name
    queue_name = fault.queue
    if queue_name is None:
        return "No queue name associated"
    if _return is None:
        _return = list()

    # Get the Location
    faulty_system = fault.location
    parent_system = mcs.rs_dir(fault.location)
    parent_system = "MC" if parent_system == "__" else parent_system

    # Get the Level
    fault_level = str(fault.level)

    # Grab the queue doc for more details
    try:
        doc, _, _ = query_document(parent_system, COL_QUEUE, Q_NAME, queue_name)
    except cexc.DatabaseRequestError as dre:
        return dre
    if not isinstance(doc, dict):
        return "DB did not return a dictionary when queue was requested"
    existing_faults = [mcs.Fault.loadd(fault_dict) for fault_dict in doc.get(Q_RECORD, list())]
    if fault in existing_faults:
        return "DB already has fault"

    # Get the step number
    step_num = get_step_number(doc)

    try:
        # Get the associated well plate nickname, and a copy of it
        container_nickname = doc[Q_OPERATIONS_LIST].get(step_num, {}).get(QOP_CONTAINER, '(no plate)')
        container_details = doc[Q_CONTAINERS].get(container_nickname, '(no details)')

        # Get failed operation
        failed_operation = doc[Q_OPERATIONS_LIST].get(step_num, {}).get(QOP_OPERATION, '(no operation)')
    except KeyError as ke:
        return ke

    file_name = f"{queue_name} ({step_num.rjust(2, '0')} - {faulty_system} - {fault_level}).log"
    _return.append(file_name)
    file_path = pathfinder("Logs")
    with safe_open((file_path, "Faults", file_name), 'a+') as fh:
        fh.write(str(fault) + f" - {failed_operation} on {container_nickname}\n")
        if fault_level != mcs.V_BUSY:
            yaml.dump(container_details, fh)
        fh.write(f"\n{'-'*64}"*3 + "\n")

    return


def _loc_sort(x, index, _type):
    """ Helper function for sorting, called by python's sort() function's key keyword argument

    :param x: object being evaluated
    :param index: Linearized index
    :param _type: datatype for the comparison
    :return: the value of _type(x[index])
    """
    location = x.get(DBG_LOCATION, None)
    if location is None:
        return _type(0)
    sort_values = [location[0], ]
    sort_values += (location[1:2] or [[]])[0]
    try:
        return _type(sort_values[index])
    except IndexError:
        return _type(0)
    except ValueError:
        return _type(0)


def add_container_to_queue(queue_name, container_ekename, mode=None, username='MC'):
    """ Adds a container to a queue document

    :param queue_name: Name of the queue
    :param container_ekename: Nickname of the container that the queue uses
    :param mode: None (Default) - Do nothing, (name, type, autofill, append) - See 'add_container_arg_helper()'
    :param username: Logon for DB
    :return: None - Failed to access DB, True - container already present, False - container not already present,
    str - nickname the container was added under
    """
    # mode: None, (name, type, autofill, False), (name, type, autofill, True)
    #       None - If already present, do nothing
    #       (name, type, autofill, False) - If already present, overwrite
    #       (name, type, autofill, True) - If already present, increment
    try:
        doc, _, _ = query_document(username, COL_QUEUE, Q_NAME, queue_name)
    except cexc.DatabaseRequestError:
        return None
    original_doc = deepcopy(doc)
    container_key = container_ekename
    if container_ekename in doc[Q_CONTAINERS]:
        if mode is None:
            # If already present, do nothing
            return True
        else:
            # If already present, ...
            if not mode[3]:
                # OVERWRITE
                # container_key = container_ekename
                pass
            else:
                # APPEND
                container_key = get_next_available_key(doc[Q_CONTAINERS], container_ekename)
                doc[Q_CONTAINERS][container_key] = dict()

            doc[Q_CONTAINERS][container_key][DBG_CONTAINER_NAME] = mode[0]
            doc[Q_CONTAINERS][container_key][DBG_CONTAINER_PLATETYPE] = mode[1]
            # If autofill flag is True, ...
            if mode[2]:
                # Autofill
                hit, _ = get_plate_details(mode[0], username)
                if hit:
                    doc[Q_CONTAINERS][container_key][DBG_CONTENTS] = hit[DBG_CONTENTS]
    else:
        if mode is None:
            # If not already present, do nothing
            return False
        else:
            # If not already present, add it
            doc[Q_CONTAINERS].setdefault(container_ekename, dict())
            doc[Q_CONTAINERS][container_ekename][DBG_CONTAINER_NAME] = mode[0]
            doc[Q_CONTAINERS][container_ekename][DBG_CONTAINER_PLATETYPE] = mode[1]
            doc[Q_CONTAINERS][container_ekename][DBG_CONTENTS] = dict()
            # If autofill flag is True, ...
            if mode[2]:
                # Autofill
                hit, _ = get_plate_details(mode[0], username)
                if hit:
                    doc[Q_CONTAINERS][container_ekename][DBG_CONTENTS] = hit[DBG_CONTENTS]
    try:
        update_document(username, COL_QUEUE, original_doc, doc)
    except cexc.DatabaseRequestError:
        return None
    return container_key


def add_container_arg_helper(container_name: str, container_type: str, autofill: bool, append: bool):
    """ Helper for the args of the add_container_to_queue method

    (May no longer be needed given improvements to Python IDE hints)

    :param container_name: Container name
    :param container_type: Container type
    :param autofill: True - Contents will be autofilled from collections, False - Contents will be unaltered
    :param append: True - uses get_next_available_key(), False - Overwrites existing entry
    :return: (container_name, container_type, autofill, append)
    """
    return container_name, container_type, autofill, append


def get_next_available_key(dictionary: dict, original_key: str, use_queue=()):
    """ creates an available key for a dictionary

    Designed originally for making container nicknames for queue documents

    :param dictionary: A dict()
    :param original_key: A string-type dictionary key
    :param use_queue: Supercedes dictionary, will pull a queue document from the database (q_name, username)
    :return: the original key or the original key appended with "__#"; (or None if use_queue and the lookup fails)
    """
    if use_queue:
        try:
            doc, _, _ = query_document(use_queue[1], COL_QUEUE, Q_NAME, use_queue[0])
        except cexc.DatabaseRequestError:
            return None
        except IndexError:
            return None
        else:
            dictionary = doc[Q_CONTAINERS]

    if original_key not in dictionary:
        return original_key
    else:
        suffix = 1
        while original_key + f"__{suffix}" in dictionary:
            suffix += 1
        return original_key + f"__{suffix}"


def attach_dependency(older_sibling_queue_name, younger_sibling_queue_name):
    """ Searches the database and adds a queue dependency (younger sibling) to any queue dependent
    on another queue (older sibling)

    Any queue dependent on 'older sibling' is now also dependent on 'younger sibling'

    :param older_sibling_queue_name: A queue name to be searched in the queue collection
    :param younger_sibling_queue_name: A queue name to be appended to the queue dependency list of all queues dependent
    on the older sibling
    :return: Generator (queue_name, append_dependency())
    """
    manifest = list()
    # Pull DB
    try:
        doc, _, _ = query_collection("MC", COL_QUEUE)
    except cexc.DatabaseRequestError as dre:
        manifest.append(("init", dre))
        return manifest
    else:
        manifest.append(("init", None))

    # For every queue in the queue collection
    for _, v in doc.items():
        # Grab the queue's name and dependencies
        dependencies = v.get(Q_DEPENDENCY, None)
        target_queue_name = v[Q_NAME]
        # print(f"checking {target_queue_name}'s dependencies ({dependencies}) for '{older_sibling_queue_name}'")
        # If the older sibling is in the dependency list, add the younger sibling to that dependency list
        if dependencies is None:
            continue
        elif isinstance(dependencies, str):
            if dependencies == older_sibling_queue_name:
                manifest.append(
                    (target_queue_name,
                     append_dependency(target_queue_name, younger_sibling_queue_name, dependencies))
                )
        elif isinstance(dependencies, list):
            if older_sibling_queue_name in dependencies:
                manifest.append(
                    (target_queue_name,
                     append_dependency(target_queue_name, younger_sibling_queue_name, dependencies))
                )
        else:
            manifest.append((target_queue_name, False))

    return manifest


def append_dependency(target_queue_name: str, dependant_queue_name: str,
                      pre_loaded_old_value: Union[None, str, list] = Narg()):
    """ Adds a queue (dependant) to the dependency list of a target queue

    :param target_queue_name: Name of queue to which a dependency is being added
    :param dependant_queue_name: Name of the queue being added to the dependency list
    :param pre_loaded_old_value: (Optional) Provide the old dependency value to skip redundant requests
    :return: True (success), False (failed), or caught Exception object
    """
    # Check if the user provided an old value
    if isinstance(pre_loaded_old_value, Narg):
        try:
            doc, _, _ = query_document("MC", COL_QUEUE, Q_NAME, target_queue_name)
        except cexc.DatabaseRequestError as dre:
            return dre
        else:
            pre_loaded_old_value = doc.get(Q_DEPENDENCY, None)

    # Add the dependant (type matching)
    if pre_loaded_old_value is None:
        new_dependency_list = dependant_queue_name
    elif isinstance(pre_loaded_old_value, str):
        new_dependency_list = [pre_loaded_old_value, dependant_queue_name, ]
    elif isinstance(pre_loaded_old_value, list):
        new_dependency_list = pre_loaded_old_value + [dependant_queue_name, ]
    else:
        return False

    # Update the field
    try:
        update_field("MC", Q_DEPENDENCY, target_queue_name, pre_loaded_old_value, new_dependency_list)
    except cexc.DatabaseRequestError as dre:
        return dre
    else:
        return True


def get_step_number(v):
    """ Gets the step number for the first incomplete step in a queue

    :param v: A queue document such that v[Q_OPERATIONS_LIST][step_num][QOP_COMPLETED] == 'no'/'yes'
    :return: str(step number) or "0" if complete
    """
    try:
        return str(
                min(
                    [int(step_num) for step_num in v[Q_OPERATIONS_LIST].keys()
                     if v[Q_OPERATIONS_LIST][step_num][QOP_COMPLETED] == 'no']
                )
            )
    except ValueError:
        return "0"


def pretty_format_queue(queue_doc, *, show_containers=True, trim_containers=False):
    """ Prints a queue document in a nicer format than default yaml

    :param queue_doc: A queue document (dictionary)
    :param show_containers: Boolean if containers should be shown
    :param trim_containers: Boolean if empty or non-referenced containers should be omitted
    :return: A printout of the queue document
    """
    # Header (name, status, date created)
    queue_name = queue_doc.get(Q_NAME, "N.R.")
    queue_status = queue_doc.get(Q_STATUS, "N.R.")
    queue_date = queue_doc.get(DBG_DATE_CREATED, "N.R.")
    header = f"Queue Name:   {queue_name}\n" \
             f"Status:       {queue_status}\n" \
             f"Date Created: {queue_date}\n"

    # Subheader (dependencies, next_step)
    queue_dependencies = queue_doc.get(Q_DEPENDENCY, None)
    if queue_dependencies is None:
        pass
    elif isinstance(queue_dependencies, list):
        pass
    elif isinstance(queue_dependencies, str):
        queue_dependencies = [queue_dependencies, ]
    else:
        queue_dependencies = f"<Improper Format ({type(queue_dependencies)})>"
    if queue_dependencies:
        subheader = f"Dependent on: {', '.join(queue_dependencies)}\n"
    else:
        subheader = f""
    next_step = get_step_number(queue_doc)
    if next_step == "0":
        subheader += "All steps Complete\n"
    else:
        subheader += f"Next Step:    {next_step}\n"

    # Operation list
    operations_list_repr = "Operations List:\n"
    operations_list = queue_doc[Q_OPERATIONS_LIST]
    sorted_step_nums = sorted(operations_list.keys(), key=lambda x: int(x))
    container_list = list()
    for k in sorted_step_nums:
        this_op = operations_list[k]
        agent = this_op[QOP_AGENT]
        operation = this_op[QOP_OPERATION]
        container = this_op.get(QOP_CONTAINER, None)
        if container in ['none', 'None', 'null', '', '--']:
            container = None
        if container is None:
            container = ""
        if int(k) >= int(next_step):
            flag = "*"
        else:
            flag = " "
        operations_list_repr += f" {flag}{k}:\t{agent}, {operation}({container})\n"
        op_details = deepcopy(this_op.get(QOP_DETAILS, dict()))
        if op_details:
            stream = StringIO()
            yaml.safe_dump(op_details, stream)
            op_details_repr = stream.getvalue()
            _count = op_details_repr.count("\n") - 1
            del stream
            spacer = "     \t      "

            op_details_repr = op_details_repr.replace("\n", "\n" + spacer, _count)
            op_details_repr = spacer + op_details_repr
            operations_list_repr += op_details_repr
        if container:
            if container not in container_list:
                container_list.append(container)

    # containers
    if show_containers:
        containers = deepcopy(queue_doc.get(Q_CONTAINERS, dict()))
        if trim_containers:
            for k in list(containers.keys()):
                if (not containers[k][DBG_CONTENTS]) or (k not in container_list):
                    del containers[k]
        if not containers:
            containers_repr = ""
        else:
            stream = StringIO()
            yaml.safe_dump(containers, stream)
            containers_repr = stream.getvalue()
            del stream
    else:
        containers_repr = ""

    return header + subheader + operations_list_repr + containers_repr


def redux_queue_as_clone(queue_name, new_status=DBQS_IDLE):
    """ Helper method to make a copy of a queue document that is reset

    :param queue_name: The parent queue document's name
    :param new_status: (Default: Idle)
    :return: Error message or "Success"
    """
    try:
        doc, _, _ = query_collection("MC", COL_HISTORICAL)
    except cexc.DatabaseRequestError:
        return 'DB access failed (1)'
    requested_queue = functionals.search_dict_by_value(doc, Q_NAME, queue_name)
    if requested_queue is None:
        return 'Queue not found'
    try:
        del requested_queue[DBG_ID]
        requested_queue[Q_STATUS] = new_status
        requested_queue[Q_NAME] += "_redux"
        for step in requested_queue[Q_OPERATIONS_LIST]:
            requested_queue[Q_OPERATIONS_LIST][step][QOP_COMPLETED] = DB_NO
        warning_flag = False
        for container in requested_queue[Q_CONTAINERS]:
            req_q_cont = requested_queue[Q_CONTAINERS][container].get(DBG_CONTENTS, "Empty")
            if not isinstance(req_q_cont, str):
                print(f"Warning, wellplate with contents: {container}")
                warning_flag = True
        if warning_flag:
            user_resp = input("Detected wellplate(s) with contents, a redux of the queue is likely unstable\n"
                              "Are you sure you want to make a revived version of this queue? (y/n)")
            if user_resp.lower() in ['y', 'yes', 'ok', 'affirmative', ]:
                pass
            else:
                return "Aborted by user"
    except KeyError as ke:
        return "Key error for " + ", ".join(ke.args)
    try:
        add_document("MC", COL_QUEUE, requested_queue)
    except cexc.DatabaseRequestError:
        return 'DB access failed (2)'
    return "Success"


def silent_monitor_checkpoint(checkpoint,
                              pipes, status: mcs.Status, task_key, task_ref,
                              queue=None, step=None, operation=None,
                              logger=None
                              ):
    """
    A utility to watch for tardy checkpoint updates.

    Spawn as a daemon thread when the mainloop times out waiting for a checkpoint

    :param checkpoint: The Checkpoint object being monitored
    :param pipes: Caller's pipe collection
    :param status: Caller's status object
    :param task_key: Idempotency key of message that issued the command
    :param task_ref: "(queue_name, #step_num: operation_name)"
    :param queue: task's queue
    :param step: task's step
    :param operation: task's operation
    :param logger: system_log
    :return: None
    """
    try:
        completion, queue, step, operation = checkpoint.wait_for_tardy(queue, step, operation)

        if logger:
            logger.info(f"{task_key} has (TARDY) been updated to {checkpoint.completion}")

        if completion is None:
            mark_time(queue, step, operation, stop=True, logger=logger)
            move(queue, _from="Lookup", _to=DBQS_IDLE)
            return
        else:
            if operation != 'complete_queue':
                mark_time(queue, step, operation, stop=True, logger=logger)

        if completion:
            if logger:
                logger.info(f"{task_key} is completed (based on function returns), "
                            f"tardy responses are not checked against DB")
            status.release_checkpoint(task_key)
        else:
            if checkpoint.level == mcs.V_BUSY:
                mark_time(queue, step, operation, reset=True, logger=logger)
            try:
                functionals.handle_false_ret_obj(pipes, status, task_key, task_ref)
                # syncs faults and releases checkpoint
            except cexc.TaskExecutionError:
                move(queue, _from="Lookup", _to=DBQS_FAIL)
        move(queue, _from="Lookup", _to=DBQS_IDLE)
    except:  # noqa
        if logger:
            logger.exception(f"Silent monitor of {task_ref} encountered an exception and was terminated, "
                             f"see UI for current status of checkpoint")


def test_connection():
    """ Tests connectivity to database

    :return: "Connected" or "No connection"
    """
    try:
        get_available_queue_names()
    except cexc.BadQueueFormatError:
        return "Connected"
    except Exception as e:  # noqa
        print("Platform Database Test Exception", repr(e))
        return "No connection"
    return "Connected"


if __name__ == '__main__':
    print(test_connection())
    update_database('localhost', 8000)
    print(test_connection())
