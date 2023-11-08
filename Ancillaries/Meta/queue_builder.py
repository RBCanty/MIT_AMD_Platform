"""
@author: Ben C

TODO: May remove methods that are not used in practice to clean things up
"""


from pprint import pformat
import time
from datetime import datetime
from copy import deepcopy

import functionals
from database_constants import *
import database_interface as dbi
import custom_exceptions as cexc
from constants import WELL_IDS
import yaml


def build_step(agent, operation='void', details: dict = None, container=None, time_est=-1,
               *, completed='no', start_time=None, end_time=None):
    """ Converts arguments (**dict) into a well-formatted step

    :param agent: The system which performs the task
    :param operation: The task
    :param details: Additional details required for operation
    :param container: The container nickname used by the queue (Special case: '--' --> '')
    :param time_est: The estimate for how long the operation will take to complete (in seconds)
    :param completed: Default: 'no'
    :param start_time: Default: None
    :param end_time: Default: None
    :return:
    """
    if container == '--':
        container = None
    template = {
        QOP_COMPLETED: completed,
        QOP_AGENT: agent,
        QOP_OPERATION: operation,
        QOP_CONTAINER: container,
        QOP_TIME: time_est,
        QOP_START: start_time,
        QOP_END: end_time,
        QOP_DETAILS: dict() if details is None else details,
    }
    if QDET_PAIRED in template[QOP_DETAILS]:
        if QDET_TIME not in template[QOP_DETAILS]:
            template[QOP_DETAILS][QDET_TIME] = 0
    return template


def detail_operation_agreement_checker(queue_step):
    """ Checks that all required details are defined for a given operation

    :param queue_step:
    :return:
    """
    operation = queue_step.get(QOP_OPERATION, None)
    detail_info = OPERATIONS_DETAIL_REQUIREMENTS.get(operation, None)
    if detail_info is None:
        return
    step_details = queue_step.get(QOP_DETAILS, dict())
    for req_key in detail_info.get('required', list()):
        if req_key not in step_details:
            print(f"Warning: {operation} missing detail {req_key}")
    for d_key, d_value in detail_info.get('dependant', dict()).items():
        if d_key in step_details:
            for dependent in d_value:
                if dependent not in step_details:
                    print(f"Warning: {operation} missing detail {dependent} required by {d_key}")

    pass


def build_operations_list(operations):
    """ Constructor for a queue operations list

    Calls build_step() for each step
    If operations is dict, the keys will be reformed as sorted (1,2,4-->1,2,3)

    :param operations: A list, dictionary, tuple, (or NoneType)
    :return: A dictionary of queue steps with step numbers as keys
    """
    if not operations:
        return {'1': {'agent': 'MC', 'operation': 'complete_queue'}}
    # Convert form
    if isinstance(operations, dict):
        sorted_step_nums = sorted([int(step_number) for step_number in operations.keys()])
        operations = [operations[str(s)] for s in sorted_step_nums]
    elif isinstance(operations, list):
        pass
    elif isinstance(operations, tuple):
        operations = list(operations)
    else:
        raise ValueError(f"Unable to convert operations into usable form")
    # Add complete queue step is missing
    if operations[-1].get(QOP_OPERATION, "void") != 'complete_queue':
        operations.append({'agent': 'MC', 'operation': 'complete_queue'})
    # Build queue list
    op_list = {str(i+1): build_step(**op) for i, op in enumerate(operations)}
    return op_list


def verify_container_names(queue_doc):
    """ Determines if all containers referenced in operations are declared by queue

    :param queue_doc: a queue document
    :return: a list of queues referenced barring queues declared ([] if all referenced queues are declared)
    """
    operations_list = queue_doc.get(Q_OPERATIONS_LIST, None)
    if operations_list is None:
        raise ValueError("Queue contains no operations")

    def _grab(op_list):
        for op in op_list:
            container: str = op_list[op].get(QOP_CONTAINER, None)
            if container is None:
                return
            elif container.lower() == 'none':
                return
            else:
                yield container

    containers_referenced_by_operations = list(_grab(operations_list))
    containers_defined_by_queue = list(queue_doc.get(Q_CONTAINERS, dict()).keys())

    return functionals.list_subtract(containers_referenced_by_operations, containers_defined_by_queue)


def verify_container_identities(queue_doc):
    """ Determines if all containers referenced in the queue.containers section which are not None/'none' are in the
    database

    :param queue_doc: A queue document with a Q_CONTAINERS field
    :return: A list of containers which are referenced by the queue but not defined in the database
    """
    containers_referenced = list()
    for container_nickname, container_details in queue_doc.get(Q_CONTAINERS, dict()).items():
        container_name: str = container_details.get(DBG_CONTAINER_NAME, None)
        if container_name is None:
            pass
        elif container_name.lower() == 'none':
            pass
        else:
            containers_referenced.append(container_name)
    containers_defined = get_wellplate_containers()
    return functionals.list_subtract(containers_referenced, containers_defined)


def get_reagent_containers():
    """ Fetches a list of container names from the reagents collection

    Includes "--" at index 0

    :return: A list of container names
    """
    try:
        doc, _, _ = dbi.query_collection('MC', COL_REAGENTS)
    except cexc.DatabaseRequestError:
        return ['--', ]
    container_names = [v[DBG_CONTAINER_NAME] for _, v in doc.items()]
    container_names.insert(0, "--")
    return container_names


def get_solvent_containers():
    """ Fetches a list of container names from the solvents collection

        Includes "--" at index 0

        :return: A list of container names
        """
    try:
        doc, _, _ = dbi.query_collection('MC', COL_SOLVENTS)
    except cexc.DatabaseRequestError:
        return ['--', ]
    container_names = [v[DBG_CONTAINER_NAME] for _, v in doc.items()]
    container_names.insert(0, "--")
    return container_names


def get_wellplate_containers():
    """ Fetches a list of container names from the wellplates collection

        Includes "--" at index 0

        :return: A list of container names
        """
    try:
        doc, _, _ = dbi.query_collection('MC', COL_WELLPLATES)
    except cexc.DatabaseRequestError:
        return ['--', ]
    container_names = [v[DBG_CONTAINER_NAME] for _, v in doc.items()]
    container_names.insert(0, "--")
    return container_names


def get_queues():
    """ Retrieves all queues in the database (queue collection and historical)

    :return: the list-sum of get_available_queue_names and  get_historical_queue_names from database_interface.py
    """
    try:
        live_q = dbi.get_available_queue_names()
    except cexc.DatabaseRequestError:
        return list()
    try:
        hist_q = dbi.get_historical_queue_names()
    except cexc.DatabaseRequestError:
        return live_q
    return live_q + hist_q


def find_in_database(name, collection, cue=(DBG_CHEM_NAME, DBG_CHEM_SMILES)):
    """ Produces the document containing the item if present

    :param name: Term being searched for
    :param collection: Name of DB collection to search (solvents and reagents are the only ones implemented)
    :param cue: Iterable of terms to be searched (default is chemical_name and chemical_smiles)
    :return: DB doc containing the item or dict() if not found
    """
    if collection == COL_SOLVENTS:
        search_field = [f"contents.XX.{term}" for term in cue]
    elif collection == COL_REAGENTS:
        search_field = [f"contents.{well_id}.{term}" for well_id in WELL_IDS for term in cue]
    else:
        raise ValueError(f"Collection '{collection}' not implemented")

    try:
        doc, _, _ = dbi.query_collection('MC', collection)
    except cexc.DatabaseRequestError:
        return dict()

    for k in doc.keys():
        for search in search_field:
            f1, f2, f3 = search.split(".")
            try:
                if doc[k][f1][f2][f3] == name:
                    return doc[k]
            except KeyError:
                pass
    return dict()


def read_in_queue_file(filepath):
    """ Converts a yaml file containing queues (keys) and wellplates (keys under 'wellplates' key into
    documents which are uploaded to the database.  Will skip existing documents.

    :param filepath: A filepath which can be opened with python's open(..., 'r') functions
    :return: (Number of documents which were uploaded, Number of documents in the file)
    """
    document = None
    error_msg = None
    n_doc_uploaded = 0
    with open(filepath, 'r') as file_handle:
        try:
            document = yaml.safe_load(file_handle)
        except yaml.YAMLError as ye:
            raise ye
    if document is None:
        return n_doc_uploaded, 0
    n_docs_possible = len(document.keys()) + len(document.get(COL_WELLPLATES, dict()).keys())
    if COL_WELLPLATES in document:
        n_docs_possible = n_docs_possible - 1

    # Load in any wellplates if they're here
    wellplates = list()
    if COL_WELLPLATES in document:
        for k, v in document[COL_WELLPLATES].items():
            v.update({DBG_CONTAINER_NAME: k})
            wellplates.append(v)
    extant_wellplates = get_wellplate_containers()
    for wellplate in wellplates:
        if wellplate[DBG_CONTAINER_NAME] in extant_wellplates:
            print(f"Wellplate '{wellplate[DBG_CONTAINER_NAME]}' already exists")
            continue
        for _ in range(0, 4):
            time.sleep(0.1)
            try:
                dbi.add_document('MC', COL_WELLPLATES, wellplate)
            except cexc.DatabaseRequestError as dre:
                error_msg = dre
                continue
            else:
                n_doc_uploaded += 1
                break
        else:
            print(f"Failed to upload wellplate '{wellplate[DBG_CONTAINER_NAME]}'")
            print(error_msg)
            return n_doc_uploaded, n_docs_possible
    if COL_WELLPLATES in document:
        del document[COL_WELLPLATES]

    # Load in the queues
    queue_docs = list()
    for k, v in document.items():
        v.update({Q_NAME: k})
        queue_docs.append(v)
    extant_queues = get_queues()
    for queue_doc in queue_docs:
        if queue_doc[Q_NAME] in extant_queues:
            print(f"Queue '{queue_doc[Q_NAME]}' already exists")
            continue
        # TODO: Make sure queue has default values for other fields (status/date_created)
        print(f"Debug: Verifying container names and identities for '{queue_doc[Q_NAME]}':")
        invalid_container_names = verify_container_names(queue_doc)
        invalid_container_identities = verify_container_identities(queue_doc)
        print("All names clear" if not invalid_container_names else invalid_container_names)
        print("All identities clear" if not invalid_container_identities else invalid_container_identities)
        print("Debug: Rectifying steps...")
        treat_queue_doc(queue_doc)
        print("... rectified!")
        for _ in range(0, 4):
            time.sleep(0.1)
            try:
                dbi.add_document('MC', COL_QUEUE, queue_doc)
            except cexc.DatabaseRequestError as dre:
                error_msg = dre
                continue
            else:
                n_doc_uploaded += 1
                break
        else:
            print(f"Failed to upload queue '{queue_doc[Q_NAME]}'")
            print(error_msg)
            return n_doc_uploaded, n_docs_possible
    return n_doc_uploaded, n_docs_possible


def treat_queue_doc(queue_doc: dict):
    """ Runs the operations of a queue document through a formatter to ensure all necessary keys are present

    :param queue_doc: A queue document
    :return: None
    """
    queue_doc[Q_OPERATIONS_LIST] = build_operations_list(queue_doc.get(Q_OPERATIONS_LIST, None))
    queue_doc.setdefault(DBG_DATE_CREATED, datetime.now().strftime(DBQ_DATE_FORMAT))
    queue_doc.setdefault(Q_STATUS, DBQS_IDLE)
    queue_doc.setdefault(Q_CONTAINERS, dict())
    queue_doc.setdefault(Q_DEPENDENCY, None)
    well_contents_errors = list(validate_well_contents(queue_doc))
    if well_contents_errors:
        raise RuntimeError(f"Errors in wellplate contents:\n{pformat(well_contents_errors)}")


def make_reaction_well(well_id: str, reagents=None, solvent=None, target_product: str = "", total_volume=0):
    """ Verifies the contents of a well {well_id: {details}} for schema and types and returns a formatted version

    :param well_id: Wellplate well ID (e.g. "A1", "B12")
    :param reagents: A list of reagents [[Name, concentration, smiles], ...]
    :param solvent: A list of solvents [[Name, concentration], ...]
    :param target_product: A SMILES string of the target product
    :param total_volume: A total volume in uL
    :return: {well_id: {details}}
    """
    # Load args
    if solvent is None:
        solvent = list(list())
    if reagents is None:
        reagents = list(list())
    if not target_product:
        print(f"Warning, no target product for {well_id}")

    # Verify values
    well_id = str(well_id).strip("\"\'")
    if well_id not in WELL_IDS:
        raise ValueError(f"Well ID '{well_id}' is not a valid well ID")
    for index, reagent in enumerate(reagents):
        name, concentration, smiles = reagent
        reagents[index] = [str(name), float(concentration), str(smiles)]
    for index, solv in enumerate(solvent):
        name, ratio = solv
        solvent[index] = [str(name), float(ratio)]
    target_product = str(target_product)
    total_volume = float(total_volume)

    # Give resultant well
    return {
        well_id: {
            DBG_CONFIRM: None,
            DBG_REAGENTS: reagents,
            DBG_SOLVENT: solvent,
            DBG_TARGET: target_product,
            DBG_TOTAL_VOLUME: total_volume
        }
    }


def validate_well_contents(queue_doc: dict):
    """ Verifies properly formatted well contents for containers with contents

    :param queue_doc: A queue document
    :return: Generator Object[all exceptions raised by the well contents checker: make_reaction_well(dict)]
    """
    for container in queue_doc.get(Q_CONTAINERS, []):
        contents_field = container.get(DBG_CONTENTS, dict())
        if isinstance(contents_field, str):
            contents_field = dict()
        for well, well_contents in contents_field.items():
            try:
                contents_field[well] = make_reaction_well(**well_contents)
            except Exception as e:
                yield [well, e]


def reset_queue(queue_name):
    """ Resets completion, start & end times, and queue status to 'no', None, and 'idle' (respectively)

    :param queue_name: Name of a queue document
    :return: (queue_name, "Success"/"Failed")
    """
    try:
        queue_doc, _, _ = dbi.query_document('MC', COL_QUEUE, Q_NAME, queue_name)
    except cexc.DatabaseRequestError:
        return queue_name, "Failed"
    edited_doc = deepcopy(queue_doc)
    edited_doc[DBQ_STATUS] = DBQS_IDLE
    operations_list = edited_doc.get(Q_OPERATIONS_LIST, dict())
    for _, operation in operations_list.items():
        operation[QOP_START] = None
        operation[QOP_END] = None
        operation[QOP_COMPLETED] = 'no'
    try:
        dbi.update_document("MC", COL_QUEUE, queue_doc, edited_doc)
    except cexc.DatabaseRequestError:
        return queue_name, "Failed"
    else:
        return queue_name, "Success"


if __name__ == '__main__':
    # path = r"C:\Users\chemegrad2018\Desktop\pairedness_test_queues.txt"
    # print(read_in_queue_file(path))
    print(reset_queue("pairedness_test10"))
    print(reset_queue("pairedness_test11"))
    print(reset_queue("pairedness_test12"))

"""
list smiles
list of properties (for each: mean, width, tolerance/weighting)
how many molecules you wanna target
how many reactions are you willing to do (max # wellplates)
Monetary budget
"""
