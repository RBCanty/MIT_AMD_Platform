""" Functionals
Collection of methods used by Systems

Decorators
  - standard_exception()  : Puts a functional in a try-except with a general return_object

Macros
  - unpack_kwargs()     : Extracts original_message, callers_state_object, & child_com from kwargs
  - unpack_db_kwargs()  : Extracts queue name, step number, plate id, & agent from kwargs

Debug -> RetObj
  - echo()                 : Prints args[0] to the console
  - exception_echo()       : Raises an exception before RetObj can be returned
  - wait_debug()           : Waits 10 seconds
  - nop()                  : Does nothing
  - db_nop()               : Does nothing (operates from DB)
  - print_from_database()  : Prints (pprint) a document requested form the DB

Core Operation -> RetObj
  - execute_message()             : Takes a message and executes it according to the command list; however,
    it then handles ibid processing based on the return codes of the operation it executed.  For example
    the message may say "Read me as an operation" and so being read as an operation generates an operation,
    but we still need  to do the operation.  As such there are many "read"/"execute" pairs
  - read_as_string()              : Returns the message as a string
  - read_as_file()                : Saves the message to a file
  - read_as_operation()           : Returns the message as an Operation object
  - execute_operation()           : Executes an Operation object
  - _execute_database_operation() : execute_operation() with additional Database functionality
  - run_from_db()                 : (Copy of read_as_operation) used for book-keeping
  - proxy_run()                   : Handles the execution of an operation from a system other than the MC
  - confirmation_of_action()      : Returns the message as an update on a task
  - update_checklist_status()     : updates the caller's status object to reflect the task update
  - handle_false_ret_obj()        : Converts checkpoints into useful Return objects
  - complete_queue()              : Marks a queue as complete and handle migrating to Historical

Stata -> RetObj
  - send_state_to()                : Enqueues a status report message to the outbox
  - read_status()                  : (MC only) Updates the caller's status object to reflect the status update
  - resuscitate_queue_operation()  : Recreates the monitor for an operation
  - _handle_queue_checkpoint()     : Handles the waiting for & post-processing of a Checkpoint
  - _handle_proxy_checkpoint()     : Handles the waiting for & post-processing of a Checkpoint for a Proxy operation
  - resolve_fault()                : Removes a fault
  - add_fault()                    : Adds a fault

Network -> RetObj
  - note_network_update()       : (MC only) Notes a disconnect from the network
  - update_database()           : Changes DB configuration without restart
  - retry_update_mongo_queue()  : Attempts to update a Queue operation to compelte

User Interfacing
  - wait_for_user() -> RetObj  : Prompts the user to press OK (blocking)
  - quick_gui() -> Any         : Prompts the user for button/entry-based input
  - quick_select() -> str      : Prompts the user to select form a list of options

Helper methods
  - realtime_wait() -> bool         : Used to make while loops depend on real time
  - search_dict_by_value()          : Searches a dictionary (first-level) for a key based on a value
  - dictionary_recursive_search()   : Searches a dictionary (multi-level) for a value based on a key
  - dictionary_direct_access()      : Accesses a dictionary using a period-delimitered string of keys
  - list_subtract() -> list         : Creates a list that is the "difference" of the two constituents
  - make_interval_and_wait_times()  : Generates scheduling/monitoring times based on a time estimate
  - void() -> None                  : pass
  - sum_if()                        : Sums elements in an Iterable provided the element meets a criterion
  - assign_list()                   : Emulates pass-by-reference using a container with __set_item__
  - queue_name_to_campaign_name()   : Extracts the Campaign name from a Queue name
  - find_shifted_step_number()      : Given the original step number and operation, attempts to find it again in an
    edited queue document and return the new location of the step.
  - semantic_join()                 : Joins elements in an Iterable with special treatment for the last element

@author: Ben C
"""

import datetime
import json
import re
import time
import tkinter as tk
from ui_database_prompt import DatabasePrompt
from functools import reduce, wraps
from pprint import pformat
from typing import Union

import custom_exceptions as cexc
import data_repository_interface as dri
import database_interface as dbi
import mcn_status as mcs
import message as msgm
import operations as oprtn
import synchronization as sync
import ui_file_save
import ui_quick_dialog
import ui_user_wait
from args_and_kwargs import *
from constants import *
from custom_classes import Narg
from database_constants import *
from mcn_logging_manager import system_log
from mcn_queues import MCNPriorityQueue


def standard_exception(func):
    """ Decorator, used to provide general bad-return objects from a function call

    If the function decorated raises an exception that isn't handled in-place,
    this decorator will generate a return object with the location of the error
    determined by AGENT or RECIP keywords (in that order)--otherwise reporting the default
    value of "Unknown"--and repr(the exception).

    :param func: Function to be wrapped
    :return: A wrapped function
    """
    @wraps(func)
    def wrapper_function(*args, **kwargs) -> mcs.RetObj:
        try:
            return func(*args, **kwargs)
        except Exception as e:
            system_log.exception("Exception from standard_exception wrapper")
            try:
                agent = kwargs.get(oprtn.AGENT,
                                   kwargs.get(ORIGINAL_MSG,
                                              {msgm.RECIP: "Unknown"}).get(msgm.RECIP,
                                                                           "Unknown caller of message.py"))
                exc_type = type(e).__name__  # str(e.__class__).strip("<>").replace("class ", "").strip("'")
                if e.args:
                    exc_det = str(e.args[0])  # .strip("'\"")
                else:
                    exc_det = "No details provided"
                return mcs.RetObj.incomplete(agent, mcs.V_PROBLEM, f"{exc_type}:{exc_det}")
            except Exception as e:  # noqa
                system_log.exception("Wrapper 'standard_exception' encountered an exception")
                return mcs.RetObj.incomplete("__", mcs.V_PROBLEM, f"@standard_exception fault--see logs\n{repr(e)}")
    return wrapper_function


def unpack_kwargs(kwargs, d_msg=None, d_cso=None, d_cc=None):
    """ Unpacks the basal elements from kwargs

    :param kwargs: The kwargs a functionals or _module method should receive when called by the MCN
    :param d_msg: Default value for the message
    :param d_cso: Default value for the caller's state object
    :param d_cc: Default value for the child communicator
    :return: original_message, callers_state_object, child_com
    """
    try:
        if d_msg is None:
            original_message: msgm.Message = kwargs[ORIGINAL_MSG]
        else:
            original_message: msgm.Message = kwargs.get(ORIGINAL_MSG, d_msg)

        if d_cso is None:
            callers_state_object: mcs.Status = kwargs[INTERNAL]
        else:
            callers_state_object: mcs.Status = kwargs.get(INTERNAL, d_cso)

        if d_cc is None:
            child_com: dict = kwargs[CHILD_COM]
        else:
            child_com: dict = kwargs.get(CHILD_COM, d_cc)
    except KeyError as ke:
        raise RuntimeError(f"Kwargs not properly specified: {repr(ke)}")

    return original_message, callers_state_object, child_com


def unpack_db_kwargs(kwargs, d_qn: str = Narg, d_sn: str = Narg, d_pid: str = Narg, d_a: str = Narg):
    """ Unpacks DB elements of kwargs (only for functions called from DB)

    :param kwargs: The kwargs a functionals or _module method should receive when called by the MCN via DB
    :param d_qn: Default value for the queue_name
    :param d_sn: Default value for the step number
    :param d_pid: Default value for plate id
    :param d_a: Default value for agent
    :return: queue name, step number, plate id, agent
    """
    if d_qn is Narg:
        queue_name = kwargs[Q_NAME]
    else:
        queue_name = kwargs.get(Q_NAME, d_qn)

    if d_sn is Narg:
        step_num = kwargs[DBQ_STEP_NUM]
    else:
        step_num = kwargs.get(DBQ_STEP_NUM, d_sn)

    if d_pid is Narg:
        plate_id = kwargs[DBQ_PLATE_ID]
    else:
        plate_id = kwargs.get(DBQ_PLATE_ID, d_pid)

    if d_a is Narg:
        agent = kwargs[oprtn.AGENT]
    else:
        agent = kwargs.get(oprtn.AGENT, d_a)

    try:
        return queue_name, step_num, plate_id, agent
    except KeyError as ke:
        raise RuntimeError(f"Kwargs not properly specified for a database operation {repr(ke)}")


@standard_exception
def echo(*args, **__) -> mcs.RetObj:
    """ Debugging, prints args[0] to console

    :param args: args[0] should be a printable object (like a string)
    :return: Standard code: (bool - success, string - details)
    """
    print(f"Debug.Echo: {args[0]}")
    return mcs.RetObj.complete("Echo Complete")


@standard_exception
def exception_echo(*args, **__):
    """ Debugging, prints args[0] to console

    :param args: args[0] should be a printable object (like a string)
    :return: Standard code: (bool - success, string - details)
    """
    print(f"Debug.Echo: {args[0]}")
    raise Exception(f"This is a test of exception handling")


@standard_exception
def wait_debug(*_, **__) -> mcs.RetObj:
    """ Debugging, waits 10 seconds

    :return: Standard code: (bool - success, string - details)
    """
    time.sleep(10)
    return mcs.RetObj.complete("Waited 10 seconds")


@standard_exception
def nop(*_, **__) -> mcs.RetObj:
    """ No Operation (an empty field)

    :return: Success
    """
    return mcs.RetObj.complete("Did nothing")


@standard_exception
def db_nop(*_, **kwargs) -> mcs.RetObj:
    """ No Operation (an empty field)

    :return: Success
    """
    q_name, q_step, _, agent = unpack_db_kwargs(kwargs, d_pid=None)  # noqa
    details = kwargs.get(QOP_DETAILS, dict())
    sleep_time = details.get('wait', 0)
    time.sleep(sleep_time)
    # try:
    #     # system_log.info(f"Void: marking {q_name} step {q_step} as complete")
    #     dbi.update_mongo_queue(mcs.rs_dir(agent), q_name, q_step, DB_YES)
    # except cexc.DatabaseRequestError:
    #     system_log.exception(f"Failed to update DB for completion of ({q_name}, {q_step}) by '{mcs.rs_dir(agent)}'")
    #     return mcs.RetObj.incomplete(agent, mcs.V_PROBLEM, 'db_nop failed to update DB for completion')
    return mcs.RetObj.complete("Did nothing")


@standard_exception
def print_from_database(*_, **kwargs) -> mcs.RetObj:
    """ Calls pretty print on a database document from the MCN control network

    The collection is specified using the 'container' keyword

    If no details are provided, it will query the collection.  If details are present, it will expect a 'field' and
    'term' keyword for querying a document.

    :param kwargs: Contains database lookup information
    :return: A return object
    """
    internals = kwargs[INTERNAL]
    q_name, q_step, _, _ = unpack_db_kwargs(kwargs, d_qn=None, d_sn=None, d_pid=None, d_a=None)  # noqa
    name = internals.name
    collection = kwargs.get('container', None)
    if collection is None:
        return mcs.RetObj.complete("missing 'collection' kwarg")
    details = kwargs.get(QOP_DETAILS, None)
    if details:
        field = details.get('field', None)
        term = details.get('term', None)
        if field is None or term is None:
            return mcs.RetObj.complete("missing 'field' or 'term' kwargs")
        try:
            doc, _, _ = dbi.query_document(name, collection, field, term)
        except cexc.DatabaseRequestError:
            system_log.exception("Failed to make request to DB")
            return mcs.RetObj.incomplete(name, mcs.V_PROBLEM, "Failed to connect to DB")
        else:
            system_log.info(f"Readout from DB:\n{pformat(doc)}")
    else:
        try:
            doc, _, _ = dbi.query_collection(name, collection)
        except cexc.DatabaseRequestError:
            system_log.exception("Failed to make request to DB")
            return mcs.RetObj.incomplete(name, mcs.V_PROBLEM, "Failed to connect to DB")
        else:
            system_log.info(f"Readout from DB:\n{pformat(doc)}")

    return mcs.RetObj.complete("Print from DB successful")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


def execute_message(message: msgm.Message, data_pipes: dict, *, cmd_list: dict, status: mcs.Status):
    """ Used for interpreting Messages as commands

    **Note**: Most Messages are, in a sense, interpreted twice (only Read_String and Read_File are interpreted once).
    In general, a Message is used to create an Operation---this Operation is stored in the data attribute of the Return
    Object (RetObj) yielded by the interpretation of the Message.  This Operation is then executed and produces its own
    RetObj.

    So the workflows are:
     - (A) Message -> Dispatcher -> RetObj
     - (B) Message -> Dispatcher -> temporary RetObj containing an Operation -> Dispatcher -> RetObj
     - (C) Message -> Dispatcher -> temporary RetObj containing an Operation -> Further in-place processing -> RetObj

    :param message: The message being interpreted
    :param data_pipes: The data pipes of the caller
    :param cmd_list: The functions available to the caller (dictionary: "string name" : function handle)
    :param status: The state object of the caller
    :return: None (Called in Thread, cannot provide useful return)
    """
    system_log.info(pformat(message.prints()))

    # Determine the intent of the Message
    message_type = message[msgm.CONTENTS]
    function = cmd_list[message_type]

    # Generate default Keyword Arguments
    positional_arguments = tuple()
    keyword_arguments = dict()
    keyword_arguments[ORIGINAL_MSG] = message
    keyword_arguments[INTERNAL] = status
    keyword_arguments[MSG_Q_OUT] = data_pipes[MSG_Q_OUT]
    keyword_arguments[MSG_Q_IN] = data_pipes[MSG_Q_IN]
    keyword_arguments[CHILD_COM] = data_pipes[CHILD_COM]

    # Execute the Function indicated by the Message's Contents
    ret_obj: mcs.RetObj = function(*positional_arguments, **keyword_arguments)

    # Examine the result of this function call...
    if ret_obj.completion:
        # ...If the functional call was successful, then (A) we're done, (B) the Message's Data was processed into
        #   an Operation object which needs to be executed, or (C) the Message's Data was processed into the input
        #   for further processes which can be done here [CrossRef: A, B, and C with the doc-string's workflows].
        msg_result_data = ret_obj.data
        # Flow A
        if message_type == 'ib.Read_String':
            pass
        # Flow A
        elif message_type == 'ib.Read_File':
            pass
        # Flow B
        elif message_type == 'ib.Read_Operation':
            keyword_arguments.update(msg_result_data.__dict__())
            ret_obj = execute_operation(msg_result_data, cmd_list, **keyword_arguments)
        # Flow B
        elif message_type == '__.Run_From_DB':
            keyword_arguments.update(msg_result_data.__dict__())
            ret_obj = _execute_database_operation(status.name, msg_result_data, keyword_arguments, data_pipes, cmd_list)
        # Flow C
        elif message_type == 'ib.Resuscitate':
            ret_obj = resuscitate_queue_operation(**keyword_arguments)
        # Flow C
        elif message_type == 'ib.Confirmation':
            ret_obj = update_checklist_status(msg_result_data, status)
        # Flow C* (*the processing is done on another system)
        elif message_type == 'ib.Run_Local':
            try:
                ret_obj = proxy_run(_pipes=data_pipes,
                                    _name=status.name,
                                    _recip=msg_result_data[oprtn.AGENT],
                                    _type="ib.Read_Operation",
                                    _operation=msg_result_data,
                                    _wait=30,
                                    _timeout=3*60)
            except cexc.ConfirmationTimeout as cte:
                system_log.exception(f"Failed to proxy an operation ({cte.get_task_key()}) within time limit")
                ret_obj = mcs.RetObj.incomplete(msg_result_data[oprtn.AGENT], mcs.V_FATAL, "Proxy run timeout")
        else:
            pass
    else:
        # ...If the functional call was unsuccessful, then log the failure
        system_log.exception(f"{status.name} failed to execute message '{message[msgm.IDKEY]}'."
                             f"Encountered exception: {ret_obj.data}")

    system_log.debug(f"Executed '{message.prints()}' with return object: {str(ret_obj)}")

    if message[msgm.PRIORITY] == msgm.P_CMD_RR:
        # This message requires a response code to be written
        response = msgm.Message.confirm_ro(message[msgm.RECIP], message[msgm.SENDER], message[msgm.IDKEY],
                                           was_successful=ret_obj.completion,
                                           location=ret_obj.location,
                                           level=ret_obj.level,
                                           info=ret_obj.data)
        data_pipes[MSG_Q_OUT].enqueue(response)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


@standard_exception
def read_as_string(*_, **kwargs) -> mcs.RetObj:
    """ Message is a string

    :param kwargs: Contains the original message and the caller's state object
    :return: Standard code (see return_object function)
    """
    message, _, _ = unpack_kwargs(kwargs, d_cso="", d_cc="")
    msg_data = message.read_data()
    system_log.info(f"Received message: {msg_data}")
    return mcs.RetObj.complete(msg_data)


@standard_exception
def read_as_file(*_, **kwargs) -> mcs.RetObj:
    """ Message is a file

    :param _: args[0] optional location for filename if not in kwargs (sender's filename)
    :param kwargs: Contains the original message, the caller's state object, and the desired filename
    :return: Standard code: Standard code (see return_object function)
    """
    message, _, _ = unpack_kwargs(kwargs, d_cso="", d_cc="")
    # filename = kwargs.get(KWA_FILENAME, args[0])

    prompt = ui_file_save.FileSaver(tk.Tk())
    # filename = prompt.run() if prompt.run() else filename
    filename = prompt.run()

    if not filename:
        system_log.warning("User failed to specify a file destination")
        return mcs.RetObj.incomplete("_U", mcs.V_PROBLEM, "Failed to provide file destination")

    successful_ftp = message.read_data(filename=filename)
    if successful_ftp:
        system_log.info(f"Received file\n"
                        f"\t\tfrom '{message[msgm.SENDER]}' and saved to\n"
                        f"\t\t{filename}.")
    else:
        system_log.warning(f"Failed to receive file from '{message[msgm.SENDER]}'.")
        return mcs.RetObj.incomplete(mcs.rs_dir(message[msgm.RECIP]), mcs.V_PROBLEM, "Failed to receive a file")
    return mcs.RetObj.complete(successful_ftp)


@standard_exception
def read_as_operation(*_, **kwargs) -> mcs.RetObj:
    """ Message is an operation

    :param kwargs: Contains the original message and the caller's state object
    :return: Standard code: Standard code (see return_object function)
    """
    message, _, _ = unpack_kwargs(kwargs, d_cso="", d_cc="")
    msg_data = message.read_data()
    operation = oprtn.Operation.build_from_json(json.loads(msg_data))
    system_log.info(f"Received operation: {operation.prints()}")
    return mcs.RetObj.complete(operation)


@standard_exception
def execute_operation(operation, cmd_list: dict, **kwargs) -> mcs.RetObj:
    """ Runs an operation

    :param operation: Operation object to be run
    :param cmd_list: Functions available to the caller
    :param kwargs: Contains the original message and the caller's state object
    :return: Standard code: Standard code (see return_object function)
    """
    function = cmd_list[operation[oprtn.FUNC]]
    positional_arguments = operation[oprtn.ARGS]
    kwargs.update(operation[oprtn.KWARGS])
    ret_obj: mcs.RetObj = function(*positional_arguments, **kwargs)
    # Operations requiring DB support get additional keywords:
    #   {'queue_name': q_name, 'plate_id': p_id, 'line': task_no}
    operation[oprtn.RETVAL] = ret_obj.completion  # unfortunate
    system_log.info(f"Recent call on '{function.__name__}' returned {str(ret_obj)}")
    return ret_obj


@standard_exception
def _execute_database_operation(operator_name: str, directive: oprtn.Operation, keyword_arguments: dict,
                                data_pipes: dict, cmd_list: dict) -> mcs.RetObj:
    """ Wrapper for 'execute_operation' which handles interfacing with the database

    :param operator_name:  The name of the System
    :param directive: An Operative derived from a Message
    :param keyword_arguments: kwargs for the Operation
    :param data_pipes: Connects execution to the System across Threads
    :param cmd_list: Allowed operations
    :return: Return Object
    """
    source_queue_name = keyword_arguments[oprtn.KWARGS][Q_NAME]
    queue_line = keyword_arguments[oprtn.KWARGS][DBQ_STEP_NUM]
    plate_id = keyword_arguments[oprtn.KWARGS][DBQ_PLATE_ID]
    agent = keyword_arguments[oprtn.AGENT]
    function_name = directive[oprtn.FUNC]

    system_log.debug(f"{plate_id} assigned to {agent}")
    # these 'details' and 'contents' are from when the queue document was imported, some processes may have
    # changed these values, so we want to grab the most recent values before performing the action
    try:
        refreshed_document, _, _ = dbi.query_document(operator_name, 'queue', 'queue_name', source_queue_name)
        updated_details = refreshed_document[Q_OPERATIONS_LIST][str(queue_line)].get(QOP_DETAILS, {})
        keyword_arguments.update({QOP_DETAILS: updated_details})
    except cexc.DatabaseGeneralException:
        system_log.exception(f"Method executor unable to update queue task information before execution")

    dbi.mark_time(source_queue_name, str(queue_line), directive[oprtn.FUNC], start=True, logger=system_log)
    ret_obj: mcs.RetObj = execute_operation(directive, cmd_list, **keyword_arguments)
    fault_timestamp = None

    # # # # start: Allow directed Recovery of the DB-update-completed Error \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    if ret_obj.completion and (function_name != "complete_queue"):
        try:
            dbi.update_mongo_queue(agent, source_queue_name, queue_line, dbi.DB_YES)
        except cexc.DatabaseRequestError:
            db_err_msg = f"DB request failed - Failed to mark ({source_queue_name}, #{queue_line}) complete"
            system_log.exception(db_err_msg)
            db_error_obj = mcs.RetObj.incomplete(agent,
                                                 mcs.V_FATAL,
                                                 f"{db_err_msg}\n----original return----\n{ret_obj.data}")
            temp_fault = mcs.Fault.from_retobj(db_error_obj, queue=source_queue_name)
            fault_timestamp = temp_fault.timestamp
            sync.synchronized_fault_add(data_pipes, temp_fault)
            response = [False, ]
            while True:
                uie = quick_gui(title="User Resolution Requested",
                                dialog=f"{db_err_msg}\n\n"
                                       f"Please mark complete once able\n"
                                       f"Use 'Abort' to redirect to Fatal Return Object",
                                buttons={'OK': lambda *_: assign_list(response, 0, True),
                                         'Abort': lambda *_: assign_list(response, 0, None),
                                         'Retry': lambda _a=agent, _q=source_queue_name, _s=queue_line:
                                         retry_update_mongo_queue(_a, _q, _s, response)
                                         },
                                ret_if_ok="OK")

                if isinstance(uie, mcs.RetObj):  # The UI encountered an unhandled exception, abort
                    ret_obj = uie
                    break

                if response[0] is False:  # User clicked RETRY but it failed
                    continue  # try again
                elif response[0]:  # User clicked OKAY or Retry was Successful
                    sync.synchronized_remove_fault(data_pipes, temp_fault)
                elif response[0] is None:  # User clicked ABORT
                    ret_obj = db_error_obj
                break
    # # # # end:   Allow directed Recovery of the DB-update-completed Error ///////////////////////////////////////////

    # If fault hasn't been locally reported yet, report it (unless it was just busy)
    if ret_obj.completion or (ret_obj.level == mcs.V_BUSY):
        pass
    else:
        new_fault = mcs.Fault.from_retobj(ret_obj, queue=source_queue_name, timestamp=fault_timestamp)
        sync.synchronized_fault_add(data_pipes, new_fault)

    return ret_obj


@standard_exception
def resuscitate_queue_operation(*_, **kwargs) -> mcs.RetObj:
    """ Recreates the monitoring of an operation

    :param kwargs: Contains System pipes and DB information
    :return: A Return Object
    """
    og_msg, status, child_com = unpack_kwargs(kwargs)
    task_key = kwargs['checkpoint_key']
    _pipes = {INTERNAL: status,
              MSG_Q_OUT: kwargs[MSG_Q_OUT],
              MSG_Q_IN: kwargs[MSG_Q_IN],
              CHILD_COM: child_com}

    queue_name = kwargs.get('queue_name', None)
    is_queue_based = True
    if (queue_name == "Internal") or (queue_name is None):
        is_queue_based = False

    # Case 1: The checkpoint exists but is deactivated, just need to reactivate
    # Case 2: The checkpoint does not exist, need to recreate it (was Queue-based)
    # Case 3: The checkpoint does not exist, need to recreate it (was Proxy-based)

    # Retrieve or Recreate Checkpoint:
    checkpoint = status.get_checkpoint(task_key)
    if checkpoint is None:
        checkpoint = mcs.Checkpoint(completion=None,
                                    location=og_msg.get('recipient', "__"),
                                    level=None,
                                    data=None,
                                    queue=queue_name,
                                    step=kwargs['queue_step'],
                                    operation=None,
                                    task_key=task_key)
        status.add_checkpoint(task_key, checkpoint)
    else:
        checkpoint.clear()  # clears the Event not the checkpoint itself
    task_ref = f"({checkpoint.queue}, #{checkpoint.step}): {checkpoint.operation})"

    if is_queue_based:
        return _handle_queue_checkpoint(checkpoint, task_key, task_ref, queue_name, status, _pipes)
    else:  # Was based on Proxy-Run
        return _handle_proxy_checkpoint(checkpoint, task_key, task_ref, kwargs, status, _pipes)


def _handle_queue_checkpoint(checkpoint: mcs.Checkpoint, task_key, task_ref: str,
                             queue_name, status: mcs.Status, _pipes: dict) -> mcs.RetObj:
    """ Handles the waiting for and post-processing of a checkpoint for a Queue-based operation

    :param checkpoint: The Checkpoint being watched
    :param task_key: The task-key of the Checkpoint
    :param task_ref: A string representation of the task
    :param queue_name: The name of the associated queue
    :param status: The Status object of the calling System
    :param _pipes: Pipes for communication with the calling System
    :return: A Return Object
    """
    # Step 1 # Mark Queue as in operation (if not already)
    try:
        dbi.move(queue_name, _from="Lookup", _to=DBQS_WORK)
    except (ValueError, cexc.DatabaseGeneralException):
        return mcs.RetObj.incomplete("__", mcs.V_PROBLEM, f"Failed to move queue '{queue_name}' to '{DBQS_WORK}'")

    # Step 2 # Await checkpoint update from executor
    checkpoint_update = checkpoint.wait_for_update(
        True,
        get_time_est=lambda: max(dbi.get_time_est(queue_name, checkpoint.step, 60), 15),
        get_inw_times=lambda x: make_interval_and_wait_times(x),
    )
    system_log.info(f"Checkpoint updated: {str(checkpoint)}")

    # Step 3 # Analyze result of checkpoint
    # User-void
    if checkpoint_update == -1:  # -> exit
        checkpoint.data = "User has suspended the checkpoint"
        return handle_false_ret_obj(_pipes, status, task_key, task_ref, _raise=False, _release=False)
    # Timeout
    elif checkpoint_update == 0:  # -> exit
        checkpoint.set()
        raise cexc.ConfirmationTimeout(f"Task '{task_ref}' has timed out", task_key=task_key)
    # checkpoint was explicitly a success or failure (opposed to user-void or timeout)
    else:
        if checkpoint.operation != 'complete_queue':
            dbi.mark_time(queue_name, checkpoint.step, checkpoint.operation, stop=True, logger=system_log)

    # Handle Successful and Failed Checkpoints
    # Failure
    if checkpoint_update != 1:  # -> exit
        try:
            if dbi.check_queue_step(queue_name, checkpoint.step) == 'yes':
                return mcs.RetObj.incomplete("_U", mcs.V_PROBLEM,
                                             f"Function returns indicated {task_ref} had failed, "
                                             f"but DB has this task marked as complete!")
        except cexc.DatabaseGeneralException:
            system_log.exception(f"Could not check {task_ref} against DB to confirm task failure")

        if checkpoint.level == mcs.V_BUSY:
            dbi.mark_time(queue_name, checkpoint.step, checkpoint.operation, reset=True, logger=system_log)
        # completion was false, boo, now handle the response
        try:
            return handle_false_ret_obj(_pipes, status, task_key, task_ref)
        except cexc.TaskExecutionError:
            try:
                dbi.move(queue_name, _from="Lookup", _to=DBQS_FAIL)
            except (ValueError, cexc.DatabaseGeneralException):
                return mcs.RetObj.incomplete('__', mcs.V_PROBLEM, f"Failed to move '{queue_name}' to {DBQS_FAIL}")
            else:
                return mcs.RetObj.incomplete(checkpoint.location, mcs.V_PROBLEM,
                                             f"{checkpoint.location} encountered a fault")

    # checkpoint_update == 1 (Success)
    system_log.info(f"{task_key} is completed (based on function returns), asking DB to confirm...")
    if checkpoint.operation == 'complete_queue':
        for _ in range(0, 4):
            if queue_name in dbi.get_historical_queue_names():
                break
            else:
                time.sleep(5)
        system_log.info(f'A complete_queue operation returned success '
                        f'but after exhausting the maximum number of DB access attempts, '
                        f'its queue document is still not in historical')
        return mcs.RetObj.incomplete("__", mcs.V_PROBLEM,
                                     f"Task {task_ref} failed to move {queue_name} to historical")
    # May need to wait for DB
    for _ in range(0, 4):
        try:
            if dbi.check_queue_step(queue_name, checkpoint.step) == DB_YES:
                # DB confirmed, exit
                system_log.info(f"DB confirmed completion of {task_ref}!")
                break
            else:
                time.sleep(5)
        except cexc.DatabaseGeneralException:
            return mcs.RetObj.incomplete("__", mcs.V_PROBLEM, 'DatabaseGeneralException in confirmation step')
    else:
        system_log.warning(f'Maximum number of attempts to to check {queue_name} for '
                           f'completion of step #{checkpoint.step} has been reached.')
        return mcs.RetObj.incomplete("__", mcs.V_PROBLEM, f"Exhausted maximum number of DB access attempts"
                                                          f" when trying to confirm step completion")

    # Update te checkpoint and queue on operation completion
    status.release_checkpoint(task_key)
    try:
        dbi.move(queue_name, _from=DBQS_WORK, _to=DBQS_IDLE)
    except (ValueError, cexc.DatabaseGeneralException):
        return mcs.RetObj.incomplete("__", mcs.V_PROBLEM,
                                     f"Failed to move queue '{queue_name}' to '{DBQS_IDLE}'")

    return mcs.RetObj.complete(f"Resuscitation of '{task_ref}' complete")


def _handle_proxy_checkpoint(checkpoint: mcs.Checkpoint, task_key, task_ref: str,
                             kwargs: dict, status: mcs.Status, _pipes: dict) -> mcs.RetObj:
    """  Handles the waiting for and post-processing of a checkpoint for a Proxy operation

    :param checkpoint: The Checkpoint being watched
    :param task_key: The task-key of the Checkpoint
    :param task_ref: A string representation of the task
    :param kwargs: Specify the waiting parameters
    :param status: The Status object of the calling System
    :param _pipes: Pipes for communication with the calling System
    :return: A Return Object
    """
    checkpoint_update = checkpoint.wait_for_update(
        True,
        get_inw_times=lambda x: (int(kwargs.get('update_interval', 15)), int(kwargs.get('timeout', 60)))
    )
    system_log.info(f"Checkpoint updated: {str(checkpoint)}")

    if checkpoint_update == -1:
        checkpoint.data = "User has voided the checkpoint"
        return mcs.RetObj.complete(f"User has voided the checkpoint {task_key}")
    elif checkpoint_update == 0:
        checkpoint.set()
        system_log.warning(f"Task '{task_key}' timed out")
        return mcs.RetObj.incomplete(checkpoint.location, mcs.V_PROBLEM,
                                     "Timeout error - from functionals.resuscitate()")
    elif checkpoint_update == 1:
        status.release_checkpoint(task_key)
        return mcs.RetObj.complete(f"Resuscitation of '{task_ref}' complete")
    else:
        return handle_false_ret_obj(_pipes, status, task_key, task_ref, _raise=False)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


@standard_exception
def send_state_to(*_, **kwargs):
    """ Sends the caller's state to a recipient

    :param kwargs: Contains the original message, the caller's state object, and the callers outbox
    :return: Standard code: Standard code (see return_object function)
    """
    message, status, _ = unpack_kwargs(kwargs, d_cc="")
    recip = kwargs.get(KWA_REQUESTER, message[msgm.SENDER])
    rsp_msg = msgm.Message.status_report(status.get_short_copy(), sender=message[msgm.RECIP], recipient=recip)
    kwargs[MSG_Q_OUT].enqueue(rsp_msg)
    return mcs.RetObj(True)


@standard_exception
def read_status(*_, **kwargs):
    """ Incorporates a status update into the caller's status

    :param kwargs: Contains the status report and the caller's state object
    :return: Standard code: Standard code (see return_object function)
    """
    requested_status = json.loads(kwargs[KWA_STS_REPORT])

    for fi, fv in enumerate(requested_status[mcs.S_FAULTS]):
        if not fv:
            continue
        try:
            requested_status[mcs.S_FAULTS][fi] = mcs.Fault.cast2fault(fv)
        except ValueError:
            system_log.exception(f"Problem loading Fault ({fi}, {fv})")
    requested_status[mcs.S_FAULTS] = [f for f in requested_status[mcs.S_FAULTS] if f]

    _, self_status, _ = unpack_kwargs(kwargs, d_msg="", d_cc="")

    self_status.short_update(requested_status)

    return mcs.RetObj.complete()


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


@standard_exception
def resolve_fault(*_, **kwargs):
    """ Removes a fault from the status object

    Note: Returns True if the error could not be found

    :param kwargs: A dict containing the caller's status object,
                   and either a fault value 'fault_value'
    :return: Standard code: Standard code (see return_object function)
    """
    _, internals, _ = unpack_kwargs(kwargs, d_msg="", d_cc="")
    fault_value = kwargs.get(KWA_FAULT_VAL, None)

    if fault_value is None:
        return mcs.RetObj.complete(f"Call to resolve_fault missing fault kwargs '{KWA_FAULT_VAL}'")

    internals.remove_fault(fault_value)

    kwargs[MSG_Q_OUT].enqueue(msgm.Message.status_report(internals.get_short_copy(),
                                                         sender=internals.name,
                                                         recipient="MC"))
    return mcs.RetObj.complete()


@standard_exception
def add_fault(*_, **kwargs):
    """ Used to manually add Faults to a system (Do not use to programmatically report faults)

    :param kwargs: A dict containing the caller's status object, and a valid Fault representation
    :return: Standard code (see return_object function)
    """
    _, internals, _ = unpack_kwargs(kwargs, d_msg="", d_cc="")
    fault_value = kwargs.get(KWA_FAULT_VAL, None)

    if fault_value is None:
        return mcs.RetObj.complete(f"Call to add_fault missing fault kwargs '{KWA_FAULT_VAL}'")

    internals.add_fault(fault_value)

    kwargs[MSG_Q_OUT].enqueue(msgm.Message.status_report(internals.get_short_copy(),
                                                         sender=internals.name,
                                                         recipient="MC"))
    return mcs.RetObj.complete()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


@standard_exception
def confirmation_of_action(*_, **kwargs):
    """ Reads Message as a confirmation receipt

    :param kwargs: Contains the original message and the caller's state object
    :return: Standard code: Standard code (see return_object function)
    """
    message, _, _ = unpack_kwargs(kwargs, d_cso="", d_cc="")
    msg_data = json.loads(message.read_data())

    system_log.info(f"Confirmed receipt and operation of {message[msgm.IDKEY]}\n"
                    f"{type(msg_data)}\n"
                    f"{msg_data}")
    return mcs.RetObj.complete(msg_data)


@standard_exception
def update_checklist_status(receipt, status: mcs.Status):
    """ Updates the caller's (MC) checklist with the confirmation given by confirmation_of_action()

    Split from confirmation_of_action() for traceback purposes

    :param receipt: The idemp key of the message being confirmed
    :param status: The state object of the caller (state being updated with new info)
    :return: Standard code: Standard code (see return_object function)
    """
    task_key = receipt[KWA_TASK_KEY]
    value = receipt[KWA_TASK_RESULT]

    if isinstance(value, mcs.RetObj):
        value = value.dumpd()

    system_log.info(f"Updating checkpoint {task_key}: " + ", ".join([f"{k} --> {v}" for k, v, in value.items()]))
    status.modify_checkpoint(task_key, **value)

    return mcs.RetObj.complete()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


@standard_exception
def note_network_update(*_, **kwargs):
    """ (MC only) Updates caller's state object that a dependant has left the network

    :param kwargs: Contains the caller's state object and the name of the disconnected system
    :return: Standard code: Standard code (see return_object function)
    """
    _, internals, _ = unpack_kwargs(kwargs, d_msg="", d_cc="")

    current_network_members = kwargs.get(KWA_NET_MEMBERS, None)

    if current_network_members is None:
        return mcs.RetObj.complete(f"note_network_update missing kwarg '{KWA_NET_MEMBERS}'")

    current_network_members = list_subtract(current_network_members, ["_S"])
    internals.update_network_members(current_network_members)
    return mcs.RetObj.complete()


@standard_exception
def update_database(*_, **kwargs):
    """ Used to change database settings without requiring restart

    :param kwargs: Contains the caller's state object and the name of the disconnected system
    :keyword ID: Expects a string, either "Platform" or "Data", for the given database
    :return: Standard code (see return_object function)
    """
    database_id = kwargs.get('ID', 'Not Specified')

    new_database_settings = DatabasePrompt(
        tk.Tk(),
        defaults=[dbi.DATABASE_URL, dbi.DATABASE_PORT]
    ).run()

    if database_id == "Platform":
        dbi.update_database(*new_database_settings)  # NOTE: uses "DBI"
    elif database_id == "Data":
        dri.update_database(*new_database_settings)  # NOTE: uses "DRI"
    else:
        system_log.info(f"Unable to find database '{database_id}', "
                        f"unable to apply network settings: {new_database_settings}")

    return mcs.RetObj.complete(data=f"Changed '{database_id}' settings to: {new_database_settings}")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


@standard_exception
def run_from_db(*_, **kwargs):
    """ Clone of execute_operation() (for now)

    :param kwargs: Contains the original message and the caller's state object
    :return: Standard code: Standard code (see return_object function)
    """
    message, _, _ = unpack_kwargs(kwargs, d_cso="", d_cc="")
    msg_data = message.read_data()
    operation = oprtn.Operation.build_from_json(json.loads(msg_data))
    system_log.info(f"Received operation: {operation.prints()} --> DB")
    return mcs.RetObj.complete(operation)


@standard_exception
def proxy_run(*, _pipes: dict,
              _name: str, _recip: str, _type: str, _operation: oprtn.Operation,
              _wait: int = 0, _timeout: int = -1):
    """ Allows a Client to issue a command to another client via the MC

    The operation object loaded should be written for the intended client not the intermediate

    The operation object's function should not be an ibid (ib.*) type message

    :raises ConfirmationTimeout: When timeout is exceeded
    :raises ValueError: When _type is not "ib.Read_Operation" or "ib.Run_Local"

    :param _pipes: Internal data pipes of the calling system
    :param _name: The sender of the message
    :param _recip: The recipient of the message (MC if client using proxy_run, AH/SP/LC/MC if MC using proxy_run)
    :param _type: If client->MC: "ib.Run_Local", if MC->client: "ib.Read_Operation"
    :param _operation: If client->MC: "ib.Read_Operation", if MC->client: the name of the operation in command_list
    :param _wait: How long (blocking) the system should wait before checking for completion
    :param _timeout: How long the system should wait after it starts checking for completion
    :return: A return object
    """
    if _type not in ["ib.Read_Operation", "ib.Run_Local"]:
        raise ValueError(f"A proxy run message type must be 'ib.Run_Local' or 'ib.Read_Operation', not '{_type}'.")
    task_key = msgm.Message.gen_idemp_key(_name)
    directive = msgm.Message.build_from_args(
        "CMD", _name, _recip, task_key, _type, _operation.package(), msgm.P_CMD_RR)
    outbox: MCNPriorityQueue = _pipes[MSG_Q_OUT]
    outbox.enqueue(directive)
    internals: mcs.Status = _pipes[INTERNAL]
    checkpoint = mcs.Checkpoint(completion=None,
                                location=_recip,
                                level=None,
                                data=_operation[oprtn.FUNC],
                                queue="Internal",
                                step="Proxy",
                                operation=_operation[oprtn.FUNC],
                                task_key=task_key)
    internals.add_checkpoint(task_key, checkpoint)
    time.sleep(_wait)

    checkpoint_update = checkpoint.wait_for_update(_timeout > 0, get_inw_times=lambda x: (1, _timeout))
    system_log.info(f"Checkpoint updated: {str(checkpoint)}")

    if checkpoint_update == -1:
        checkpoint.data = "User has suspended the checkpoint"
        return handle_false_ret_obj(_pipes, internals, task_key, directive.prints(), _raise=False, _release=False)
    elif checkpoint_update == 0:
        checkpoint.set()
        raise cexc.ConfirmationTimeout(f"Task '{directive.prints()}' did not clear within {_timeout} seconds",
                                       task_key=task_key)
    elif checkpoint_update == 1:
        internals.release_checkpoint(task_key)
        return mcs.RetObj.complete(task_key)
    else:
        return handle_false_ret_obj(_pipes, internals, task_key, directive.prints(), _raise=False)


def handle_false_ret_obj(pipes, state, task_key, ekename, _raise=True, _release=True):
    """ Manages the checkpoint and Fault reporting of Incomplete Return Objects

    :raises TaskExecutionError: On a FATAL response code (needed by scheduler)

    :param pipes: Internals for synchronous updates
    :param state: Caller's state object for logging and checkpoints (redundant with pipes)
    :param task_key: The key used to make the message and build the checkpoint
    :param ekename: How the logger should refer to the task that is being checked
    :param _raise: True (raise TaskExecutionError), False (Log FATAL fault)
    :param _release: Should the checkpoint be released?
    :return: A return object built from the checkpoint value
    """
    with state as s:
        comp, loc, level, info, q_name, *_ = s[mcs.S_SYSTEM_][mcs.S_CHECKPOINTS][task_key].unpack()
        if comp:
            ret_val = mcs.RetObj.complete(info)
            # Isn't used in any current implementations, but nice to be general
        else:
            ret_val = mcs.RetObj.incomplete(loc, level, info)

    if _release:
        state.release_checkpoint(task_key)

    # task_ref = f"({queue_name}, #{queue_step}): {queue_operation})"
    backup = re.match(r"\A\((\S*), (\S*): (\S*)\)\Z", ekename)
    if backup:
        try:
            _q, _s, _o = backup.groups()
        except ValueError:
            pass
        else:
            q_name = q_name if q_name else _q

    system_log.warning(f'Checkpoint for {ekename} failed')

    p_fault = mcs.Fault(location=loc, level=level, data=info, queue=q_name)

    if level == mcs.V_BUSY:
        system_log.info(f'Agent busy, will retry {ekename}')
    elif level == mcs.V_PROBLEM:
        sync.synchronized_fault_add(pipes, p_fault)
        system_log.info(f'Failure due to problem with the agent (problem), will retry {ekename}')
    elif level == mcs.V_FATAL:
        sync.synchronized_fault_add(pipes, p_fault)
        if _raise:
            raise cexc.TaskExecutionError(f"{ekename} failed")
        system_log.info(f'Failure due to problem with the agent (fatal): {ekename}')
    else:
        system_log.warning(f"Unrecognized fault 'level' field for a failed process: '{level}', "
                           f"will retry {ekename} later.")
        sync.synchronized_fault_add(pipes, p_fault)
    return ret_val

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


@standard_exception
def complete_queue(*_, **kwargs):
    """
    Completes a queue.
      1. Changes status from working to done.
      2. Marks the final step complete.
      3. Marks the end time of the final step.
      4. Moves the queue to Historical.

    :param _: Used to consume args passed by manager
    :param kwargs: Used for queue name and step number
    :return: mcs.RetObj.complete()
    """
    q_name, q_step, _, _ = unpack_db_kwargs(kwargs, d_pid=None, d_a=None)  # noqa
    if q_name is None:
        system_log.info(f"complete_queue() missing kwarg: '{Q_NAME}'")
        return mcs.RetObj.incomplete("MC", mcs.V_PROBLEM, "No queue_name kwarg given")

    model_details = kwargs.get(QOP_DETAILS, {}).get('model_details', None)
    completion_tasks = [
        ("Marking Queue as complete", dbi.update_field,
         ("MC", 'status', q_name, DBQS_WORK, DBQS_DONE), {},
         "status     -> complete"),

        ("Marking queue Step as complete", dbi.update_mongo_queue,
         ('MC', q_name, q_step, DB_YES), {},
         "final_step -> complete"),

        ("Marking completion time", dbi.mark_time,
         (q_name, q_step, 'complete_queue'), {'stop': True, 'logger': system_log},
         "end_time   -> now"),

        ("Posting update to data DB",  dri.complete_queue,
         (queue_name_to_campaign_name(q_name), model_details), {},
         "Notified Data Database"),

        ("Moving to Historical", dbi.move_to_historical,
         (q_name,),  {},
         "location   -> historical"),
    ]
    # If no model details, do not call step 3: 'dri.complete_queue'
    if not model_details:
        del completion_tasks[3]
    n_tasks = len(completion_tasks)
    n_attempts = 4

    start_from = kwargs.get(QOP_DETAILS, {}).get('start_from', 0)
    completion_tasks = completion_tasks[start_from:]
    caught_exceptions = []
    success_returns = []

    system_log.info(f"Performing completion of '{q_name}'...")
    for i, (_verb, _func, _args, _kwargs, _msg) in enumerate(completion_tasks, start=1):
        for attempt in range(n_attempts):
            try:
                system_log.info(_verb)
                _ret = _func(*_args, **_kwargs)
                success_returns.append((_verb, _ret))
                system_log.info(f"{q_name}: {_msg:25}({i}/{n_tasks})")
            except (cexc.DatabaseGeneralException, AttributeError, IndexError) as e:
                caught_exceptions.append(repr(e))
                msg_tag = "aborting" if attempt == (n_attempts - 1) else "will retry"
                system_log.warning(f"{q_name}: error on ({i}: {_verb}) attempt {attempt+1} of 4, {msg_tag}")
                time.sleep(1)
            else:
                break  # Success, continue to next task
        else:  # Ran out of attempts
            bad_ret_obj = mcs.RetObj.incomplete("MC", mcs.V_FATAL,
                                                f"Failed to complete queue '{q_name}' ({i}/{n_tasks}: "
                                                f"failed on '{_verb}')")
            uie = quick_gui(title="User Resolution Requested",
                            dialog=f"Problem completing queue '{q_name}'\n"
                                   + "\n".join([v[0] for v in completion_tasks[:i]])
                                   + f"\n\nPlease resolve (OK) or abort operation (Abort).",
                            buttons={'OK': lambda *_: None,
                                     'Abort': lambda *_: None},
                            ret_if_ok="OK")
            if isinstance(uie, mcs.RetObj):  # The UI encountered an unhandled exception, abort
                bad_ret_obj.data = bad_ret_obj.data + f"\n\nGUI Recovery Error:\n{uie.data}" \
                                                      f"\n\nDatabase(s) Error History:\n{caught_exceptions}"
                return bad_ret_obj
            elif uie == "OK":  # User clicked OK, trust that they resolved the issue
                pass
            else:  # User clicked Abort or X'd out, assume not resolved
                return bad_ret_obj
    system_log.info(f"... '{q_name}' complete!")

    return mcs.RetObj.complete()


def retry_update_mongo_queue(agent, source_queue_name, queue_line, output: list):
    """
    Used to invoke dbi.update_mongo_queue and report success via a list return-argument

    :param agent: The name of the system
    :param source_queue_name: The name of the queue
    :param queue_line: The step number
    :param output: A list with at least 1 element
    :return: output[0] = True (success), False (failure)
    """
    try:
        dbi.update_mongo_queue(agent, source_queue_name, queue_line, dbi.DB_YES)
    except cexc.DatabaseRequestError:
        output[0] = False
    else:
        output[0] = True

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   GUI and Utility Functions                                                                                         #
#                                                                                                                     #
# Which, honestly, belong in their own module* since the rest of functionals.py concerns itself with operations that  #
# a System can perform and with the actual execution of operations.                                                   #
# *wait_for_user is a queue-document-level command, so it should stay here.                                           #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


@standard_exception
def wait_for_user(*_, **kwargs):
    """ A simple wait dialog (Continue/Halt)

    Continue gives Return Object (True, value="User told MCN to continue")

    Halt gives Return Object (False, location="_U", level=FATAL, info="User issued Halt")

    :param kwargs: Used for SYSTEM_DATA -> logging
    :return: A Return Object (or None if X-ed out)
    """
    info = kwargs.get('info', "Waiting for User")

    prompt = ui_user_wait.Waiter(tk.Tk(), info)
    ret_val = prompt.run()

    system_log.info(f"User Wait returned: {ret_val.data}")

    return ret_val


def quick_gui(*, title: str = None, dialog: str = None, buttons: dict = None, has_entry=False, ret_if_ok=""):
    """
    Creates a simple UI without needing to know Tkinter

    By default, all UI will have an OK button (you may overwrite it by adding your own {"OK": func} to the buttons
    input argument).  When has_entry is True, an OK is replaced by Submit button and the return of run() will be
    the text in the entry field (the return is "" otherwise).  When has_entry is True, ALL functions in the buttons
    dictionary will receive the text of the entry widget as their args[0].  All buttons will close the message box
    after their referenced functions execute.

    Is decorated by @standard_exception, so if quick_gui encounters an exception, instead of returning like normal
    it will return a return_object(False, ..., info=repr(Exception))

    :param title: The text shown on the header bar of the message box
    :param dialog: The text shown in the message box
    :param buttons: A dictionary of {"Button text": Python function handle}
    :param has_entry: (Default: False) True - adds a text entry field
    :param ret_if_ok: (Default: "") If using the default OK function, then run() will return ret_if_ok if the user
      clicks OK but "" if any other button is clicked
    :return: (has_entry: False) an empty string, (has_entry=True) the text stored in the textbox
    """
    prompt = ui_quick_dialog.QuickUI(tk.Tk(),
                                     title=title,
                                     dialog=dialog,
                                     buttons=buttons,
                                     has_entry=has_entry,
                                     ret_if_ok=ret_if_ok)
    ret_val = prompt.run()
    return ret_val


def quick_select(*, title: str = None, dialog: str = None, options: list = None, default=None) -> str:
    """ A quick UI for having the user select one of a given set of options.

    Returns the user's selection (string) or the default value (string, 'None' unless set by caller)
    The default value must be handled by the caller as an option (even if that's just raising an exception)

    Example Use::
        initial_column = quick_select(
                                      title="Column Selector",
                                      dialog="Please specify what column is currently installed",
                                      options=["Analytic", "Semi-Prep"],
                                      default="Neither")

    :param title: The title of the window
    :param dialog: A prompt for the user
    :param options: A list of options (str) for the user to choose
    :param default: The default value if nothing is selected
    :return: The user's selection (or default) as a string
    """
    prompt = ui_quick_dialog.QuickSelectUI(tk.Tk(),
                                           title=title,
                                           dialog=dialog,
                                           options=options,
                                           default=default)
    ret_val = prompt.run()
    return ret_val

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


def realtime_wait(start: datetime.datetime, *, hours: int = 0, minutes: int = 0, seconds: int = 0) -> bool:
    """ logical shorthand for while (elapsed time < specified time): do stuff

    Use: "timer = datetime.now()" then "while realtime_wait(timer, [kwargs]):"

    :param start: A datetime object which remembers the start time
    :param hours: An integer representing the number of hours to wait (Default = 0)
    :param minutes: An integer representing the number of minutes to wait in addition to hours (Default = 0)
    :param seconds: An integer representing the number of seconds to wait in addition to hours and minutes (Default = 0)
    :return: True - elapsed time <= specified time / False - elapsed time > specified time
    """
    delta = datetime.timedelta(hours=hours, minutes=minutes, seconds=seconds)
    return datetime.datetime.now() - start <= delta


def search_dict_by_value(library, field, value):
    """ Utility for finding dict[k] such that dict[k][field] == value

    :param library: A parent dictionary
    :param field: The term being matched
    :param value: The value of the term being matched
    :return: The child of the parent which contains the field that matched value (or None)
    """
    for k in library.keys():
        try:
            if library[k][field] == value:
                return library[k]
        except KeyError:
            continue
    return None


def _dictionary_recursive_search(library, searched_value, depth, sub_search):
    """ recursive underbelly of the dictionary search

    :param library: A dictionary to be searched
    :param searched_value: The value being searched for
    :param depth: (used to track search depth)
    :param sub_search: Search in non-dictionary elements?
    :return: (True/False - was found, List - dictionary keys)
    """
    if depth == 0:
        return False, list()
    if type(searched_value) == type(library):
        if searched_value == library:
            return True, None
    if sub_search:
        for i, dk_lv in enumerate(library):
            if isinstance(library, dict):
                sub_key = dk_lv
            else:
                sub_key = i
            if sub_key == 0 and dk_lv == library:
                if dk_lv == searched_value:
                    return True, [sub_key, ]
                else:
                    return False, None
            is_found, path_raw = _dictionary_recursive_search(library[sub_key], searched_value, depth-1, sub_search)
            if is_found:
                if path_raw is None:
                    return True, [sub_key, ]
                else:
                    return True, [sub_key, ] + path_raw
    else:
        if isinstance(library, dict):
            for k in library:
                is_found, path_raw = _dictionary_recursive_search(library[k], searched_value, depth-1, sub_search)
                if is_found:
                    if path_raw is None:
                        return True, [k, ]
                    else:
                        return True, [k, ] + path_raw
        else:
            return False, None
    return False, list()


def dictionary_recursive_search(library: dict, searched_value, depth: int = -1, sub_search: bool = False):
    """ Utility function for finding a value within a dictionary

    :param library: A dictionary to be searched
    :param searched_value: The value being searched for
    :param depth: (-1) Search all levels, otherwise specifies maximum number of levels to search
    :param sub_search: False (Default) - Search exact values only, True - Search in non-dictionary elements
    :return: None if not found, or a list of keys to the value
    """
    is_found, path = _dictionary_recursive_search(library, searched_value, depth, sub_search)
    if is_found:
        return path
    else:
        return None


def dictionary_direct_access(dictionary: dict, key_path: str, delimiter="."):
    """
    Allows a dictionary to be accessed with a flat key string (e.g. "key_1.key_2.key_3")

    :param dictionary: The dictionary being searched
    :param key_path: A string of keys separated with a delimiter
    :param delimiter: (Default: ".") The delimiter of key_path
    :return: The value in dictionary addressed by key_path
    :raises KeyError: If any key in key_path is not present in dictionary
    """
    _key_path = key_path.split(delimiter)
    temp = reduce(dict.get, _key_path[:-1], dictionary)
    last_key = _key_path[-1]
    try:
        return temp[last_key]
    except TypeError:
        raise KeyError(_key_path) from None


def list_subtract(minuend: list, subtrahend: list):
    """ Utility function for subtracting one list from another

    :param minuend: The list providing values
    :param subtrahend: The list providing exceptions
    :return: minuend less the elements in subtrahend
    """
    return [x for x in minuend if x not in subtrahend]


def sum_if(iterable, *, start=0, condition=lambda x: bool(x), quantifier=lambda x: x):
    """ Utility function for Sum-If functionality

    :param iterable: The object being summed over
    :param start: (from sum) value to initialize sum
    :param condition: An expression which is evaluated (True = include in sum) on each element of the iterable
    :param quantifier: Used to convert each element of the iterable into a number
    :return: sum([quantifier(i) for i in iterable if condition(i)], start)
    """
    return sum([quantifier(i) for i in iterable if condition(i)], start)


def make_interval_and_wait_times(time_estimate):
    """ Helper method to convert a time_estimate into an interval and wait time.

    :param time_estimate: A duration in seconds
    :return: An appropriate interval (seconds), An appropriate wait time (seconds)
    """
    if time_estimate < 60:
        interval = 5
    elif time_estimate < 5*60:
        interval = 15
    elif time_estimate < 15*60:
        interval = 30
    else:
        interval = 60

    wait_time = max(2 * time_estimate + interval, 60) + 15

    return interval, wait_time


def assign_list(pass_by_reference: list, index: int, value):
    """  Helper method for when a return-less function needs to use a return-argument to pass a response

    :param pass_by_reference: A list or __set_item__ object
    :param index: The key provided in  pass_by_reference[index]
    :param value: The value of the keyed field
    :return: None (pass_by_reference[index] = value)
    """
    pass_by_reference[index] = value


def queue_name_to_campaign_name(queue_name: str):
    """ Converts a queue_name to a campaign_name

    All queues are named f"{campaign_name}_{index}_{date}"

    :param queue_name: The string value of a queue name
    :return: a campaign name
    """
    campaign_name, *_ = queue_name.rsplit("_", 2)
    return campaign_name


def find_shifted_step_number(operations_list: dict,
                             step_number: Union[int, str],
                             operation_name: str = None,
                             q_name: str = "<not specified>") -> str:
    """
    Searches through an operations list for the first instance of an operation on or after a given step number.
    Used to protect against when queues add steps for recovery.

    :param operations_list: A dictionary taken from the queue document.
    :param step_number: the step number to start the search from
    :param operation_name: the name of the operation being searched for
    :param q_name: The name of the queue being searched--for logging
    :return: A string of the step number that was found.
    """
    # The default behavior of these searches is that if no operation name is specified for searching, then no search is
    #   performed and the unaltered step number is used.
    if not operation_name:
        return str(step_number)

    _step = int(step_number)
    max_iter = max(int(s) for s in operations_list.keys())
    while _step <= max_iter:
        if operations_list[str(_step)][QOP_OPERATION] == operation_name:
            return str(_step)
        _step += 1
    raise cexc.BadQueueFormatError(f"Queue '{q_name}' is "
                                   f"missing a step for operation '{operation_name}' after {step_number}")


def void(*_, **__):
    """ Does nothing (on purpose)

    :return: None
    """
    pass


def semantic_join(iterable, separator: str, last_separator: str = None,
                  replace_separator_from_last_if_dual=True,
                  replace_with: str = None):
    """ Joins elements in a list with special treatment for the last join and for lists of 2 elements

    Examples:
     - A, B, ..., C, and D (separator=", ", last_separator=", and ")
     - A, B, ..., C, D (separator=", ")
     - A or B (separator=", ", last_separator=", or ")
     - A, or B (separator=", ", last_separator=", or ", replace_separator_from_last_if_dual=False)

    Based on 'format_authors' by stackoverflow user Shashank
    (May 6, 2015: https://stackoverflow.com/a/30084397/17421441)

    :param iterable: An iterable object
    :param separator: The default delimiter of items
    :param last_separator: The delimiter for the last two items (defaults to separator)
    :param replace_separator_from_last_if_dual: Replaces separator from last_separator with 'replace_with' if there
      are only two items in the iterable.  Defaults to True which will replace the separator in last_separator with " ".
    :param replace_with: A string to replace the separator in the last separator if there are only two items.  Only
      applies when both (a) replace_separator_from_last_if_dual is True, (b) separator is a substring of last_separator,
      and (c) there are two items in the iterable.  (Default value is a whitespace, " ")
    :return: A string of the list concatenated by the provided delimiters
    """
    if last_separator is None:
        last_separator = separator
    if replace_with is None:
        replace_with = " "
    if not iterable:
        return ''
    n = len(iterable)
    if n > 2:
        return (f'{{}}{separator}'*(n-2) + f'{{}}{last_separator}{{}}').format(*iterable)
    elif n == 2:
        if replace_separator_from_last_if_dual:
            last_separator = last_separator.replace(separator, replace_with)
        return (f'{{}}{separator}'*(n-2) + f'{{}}{last_separator}{{}}').format(*iterable)
    elif n == 1:
        return f'{iterable[0]}'
    else:
        return ''


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


if __name__ == '__main__':
    pass

    """
    # Example use of proxy_run()
    # Here SP is asking the Robotic arm (Ra) to put the lid on the thermo-reactor plate
    try:
        rtrnobj = proxy_run(_pipes=dict(),
                            # This is the dictionary of system.data_pipes,
                            # The thing that holds CHILD_COM, INTERNALS, MSG_Q_IN, and MSG_Q_OUT
                            _name="SP",
                            # The SP system is making this request
                            _recip="MC",
                            # This says MC because the message is routed through the MC,
                            # not because the agent (Ra) is a member of MC!  If this example
                            # were SP asking LC to do something, _recip would still be "MC".
                            _type="ib.Run_Local",
                            # This is "ib.Run_Local" in all use cases you _should_ encounter, the MC
                            # changes this to "ib.Read_Operation" when it forwards (which is why it's even an argument)
                            _operation=oprtn.Operation.build_from_args(func="load_Th_plate_lid",
                                                                       agent="Ra",
                                                                       # kwargs={},
                                                                       ),
                            # The function name "load_Th_plate_lid" is a placeholder
                            # The agent is the agent of the operation, like normal
                            # kwargs=dict() where you can specify the arguments that load_Th_plate_lid may need
                            _wait=35,
                            # I want SP to wait 35 seconds before checking if the MC said Ra finished
                            # MC waits 30 seconds by default since I don't know how to give it a good
                            # time estimate to use
                            _timeout=5 * 60
                            # After waiting _wait seconds, give the MC 5 minutes to respond to SP
                            # The MC gives the Ra 3 minutes to respond to MC.
                            # When you specify your timeout, two things that can happen here:
                            # If _timeout > 3 minutes, then you should receive a message from MC that Ra failed
                            # to respond in time and you can handle the timeout by inspecting the ret_obj given
                            # by proxy_run.  Note, this will clear your checkpoint.
                            # If _timeout < 3 minutes, then you can timeout before MC, this means you go to
                            # the following 'except ConfirmationTimeout' clause where you can handle the
                            # exception as you see fit, but this will not clear the checkpoint.  If you want
                            # to manually clear the checkpoint, you will need to get the task_key via the
                            # getter method demonstrated below.  Note: The MC may send a response (if Ra
                            # does eventually respond within 5 minutes), and the checkpoint will be updated;
                            # however, nothing will be waiting for it, so it will remain in State[S_CHECKPOINT]
                            # with updated values but no one listening to it--unless you use the task_key to
                            # look it up.
                            )
    except cexc.ConfirmationTimeout as cte:
        print(f"The timed-out task key was {cte.get_task_key()}")
    """
