""" Aceso

Library of recovery options for the MCN

(Aceso, sister of Panacea and goddess of the healing process)

Author: Ben C
"""
import custom_exceptions
import mcn_status
from queue_builder import build_step
from database_constants import *
import database_interface as dbi
import json
from constants import INTERNAL, CHILD_COM, MSG_Q_OUT, MSG_Q_IN
from threading import Thread


def resource_request(queue_name, step, target_destination: str, is_new: bool, name_basis: str, plate_type=None,
                     generate=False, directly_specify: tuple = None):
    """ Inserts steps into a queue doc's operations list to request a resource (such as a wellplate or tips)

    :param queue_name: Name of the queue document
    :param step: Step number of step who needs the request (steps get inserted prior to this index)
    :param target_destination: Where the resource is needed
    :param is_new: True - An entry is needed in the queue/ False - The entry exists in the queue already
    :param name_basis: The queue's nickname for a plate (or the basis of a series of names)
    :param plate_type: The labware type (can be None [param ignored] is is_new=False)
    :param generate: False - nickname = name_basis / True nickname = f(name_basis): f(n) = "n__##"
    :param directly_specify: Don't Use (str - "True Name of Plate", bool - True to autofill contents, False to not)
    :return: True/False (Success), The nickname it wound up using or an Error Message, Exception Object or None if n.a.
    """
    # Very special case:
    if directly_specify:
        true_name = directly_specify[0]
        autofill = directly_specify[1]
    else:
        true_name = None
        autofill = False
    # Do we need to make a new container in the queue?
    if is_new:
        # Check
        if plate_type is None:
            return False, 'Cannot create a new wellplate of type None', None
        # Make a container
        nickname_used = \
            dbi.add_container_to_queue(queue_name,
                                       name_basis,
                                       mode=dbi.add_container_arg_helper(
                                           true_name,  # Plate does not have a true name yet
                                           plate_type,  # We know the plate type
                                           autofill,  # Search wellplates for the true name and autofill contents
                                           generate  # Autogenerate a nickname if user requests (eg. tips, FC plates)
                                       )
                                       )
    else:
        nickname_used = name_basis

    try:
        queue_doc, _, _ = dbi.query_document('MC', COL_QUEUE, Q_NAME, queue_name)
    except custom_exceptions.DatabaseRequestError as dre:
        return False, "Failed to access DB", dre
    if nickname_used not in queue_doc[Q_CONTAINERS].keys():
        return False, f"Reference '{nickname_used}' not found in queue '{queue_name}'", None

    # Cleanup
    if target_destination.lower() == 'lh':
        target_destination = 'liquid_handler'
    elif target_destination.lower() == 'fc':
        target_destination = 'fraction_collector'
    step = int(step) - 1

    # BUILD
    step_1 = build_step('Ss',
                        operation='Ss.request',
                        details={QDET_PAIRED: 'yes'},
                        container=nickname_used,
                        time_est=60)
    step_3 = build_step('Ra',
                        operation='move_wellplate',
                        details={QDET_PAIRED: 'yes',
                                 'target_destination': target_destination},
                        container=nickname_used,
                        time_est=30)
    if target_destination == 'liquid_handler':
        step_2 = build_step('Lh',
                            operation='transfer_wellplate',
                            details={QDET_PAIRED: 'yes',
                                     'target_destination': 'transfer_prep'},
                            container=nickname_used,
                            time_est=30)
        step_4 = build_step('Lh',
                            operation='transfer_wellplate',
                            details={QDET_PAIRED: 'yes',
                                     'target_destination': 'transfer_cleanup'},
                            container=nickname_used,
                            time_est=80)
    elif target_destination == 'fraction_collector':
        step_2 = None
        step_4 = None
    else:
        return False, f"Resource Request not established for '{target_destination}'", None
    queue_steps = [step for step in [step_1, step_2, step_3, step_4] if step is not None]
    try:
        dbi.insert_queue_steps(queue_name, step, queue_steps)
    except custom_exceptions.DatabaseGeneralException as dge:
        return False, "Failed to insert steps", dge

    return True, nickname_used, None


def tip_request(queue_name, step: str):
    """ Convenience method for requesting tips onto the liquid handler

    :param queue_name: Name of queue
    :param step: Step number (str) of the operation which failed due to not having tips
    --> tip fetching sequence will be inserted before this step number
    :return: None
    """
    return resource_request(queue_name, step, 'liquid_handler', True, 'tips_from_lpx', 'DiTi SBS Waste', True, 'no')


def request_optimization(queue_name, step):
    raise NotImplementedError  # TODO: generalize code from Matt's LC code


def request_sample_adjustment(queue_name, step: str, container: str, cleanup_key: str, iter_num: int, n_changes: int):
    """
    For the Plate Reader to submit a plate to the Liquid Handler to dilute or concentrate wells for optimal signal

    :param queue_name: The name of the queue
    :param step: The current step (additional steps will be inserted after this step)
    :param container: The Queue nickname for the plate
    :param cleanup_key: A key for cleanup_characterization_plate's target_cleanup detail
    :param iter_num: The current iteration number (the newly added Pr/check step will use int(iter_num - 1) as its iter)
    :param n_changes: Used for estimating times
    :return: A triple (0 - True/False for success, 1 - A message, 2 - A caught exception object)
    """
    try:
        queue_doc, _, _ = dbi.query_document('MC', COL_QUEUE, Q_NAME, queue_name)
    except custom_exceptions.DatabaseRequestError as dre:
        return False, "Failed to access DB", dre

    # Move the plate back to the Lh, perform reductions & augmentations, move plate back to host
    step_1 = build_step('Pr',
                        operation='Pr.prepare_to_send_or_receive',
                        container=container,
                        details={QDET_PAIRED: 'yes',
                                 'schedule_time': 0,
                                 'direction': 'send'},
                        time_est=5)
    step_2 = build_step('Lh',
                        operation='transfer_wellplate',
                        container=container,
                        details={QDET_PAIRED: 'yes',
                                 'schedule_time': 0,
                                 'target_destination': 'bed_position'},
                        time_est=60)
    step_3 = build_step('Pr',
                        operation='Pr.go_to_initial_state',
                        container=container,
                        time_est=30)
    step_4 = build_step('Lh',
                        operation='cleanup_characterization_plate',
                        container=container,
                        details={'target_cleanup': cleanup_key},
                        time_est=60 + 120*n_changes)
    step_5 = build_step('Pr',
                        operation='Pr.prepare_to_send_or_receive',
                        container=container,
                        details={QDET_PAIRED: 'yes',
                                 'schedule_time': 0,
                                 'direction': 'receive'},
                        time_est=5)
    step_6 = build_step('Lh',
                        operation='transfer_wellplate',
                        container=container,
                        details={QDET_PAIRED: 'yes',
                                 'schedule_time': 0,
                                 'target_destination': 'spark'},
                        time_est=60)
    step_7 = queue_doc[Q_OPERATIONS_LIST][step]
    step_7.update({'completed': 'no',
                   'end_time': None,
                   'start_time': None})
    step_7.setdefault(QOP_DETAILS, dict())
    step_7[QOP_DETAILS]['iter'] = iter_num - 1

    queue_steps = [step_1, step_2, step_3, step_4, step_5, step_6, step_7]

    if iter_num == 0:
        del queue_steps[6]

    try:
        dbi.insert_queue_steps(queue_name, step, queue_steps)
    except custom_exceptions.DatabaseGeneralException as dge:
        return False, "Failed to insert steps", dge

    return True, "Steps added", None


def anabasis(file_path: str, _data_pipes: dict, resuscitate_function=None):
    """ Used to recover checkpoint data from a save file upon startup

    anabasis - an ascent, often from the underworld (system startup)

    :param file_path: Location of the save file
    :param _data_pipes: The calling system's self.data_pipes field
    :param resuscitate_function: The function used to handle resuscitation (funcationals.resuscitate_queue_operation)
    :return: A list of child threads spawned
    """
    status_object = _data_pipes[INTERNAL]
    child_threads = list()
    try:
        with open(file_path, 'r') as fh:
            try:
                past_life = json.load(fh)
            except json.decoder.JSONDecodeError:
                past_life = {}

        for k, v in past_life.items():
            checkpoint = mcn_status.Checkpoint(**v)
            if checkpoint.task_key is None:
                checkpoint.task_key = k
            status_object.add_checkpoint(k, checkpoint)
            if resuscitate_function is not None:
                _kwargs = {INTERNAL: status_object,
                           CHILD_COM: _data_pipes[CHILD_COM],
                           MSG_Q_IN: _data_pipes[MSG_Q_IN],
                           MSG_Q_OUT: _data_pipes[MSG_Q_OUT],
                           'checkpoint_key': k,
                           'queue_name': checkpoint.queue,
                           'queue_step': checkpoint.step,
                           'update_interval': 15,
                           'timeout': 5*60,
                           }
                child = Thread(target=resuscitate_function,
                               kwargs=_kwargs,
                               daemon=True)
                child.start()
                child_threads.append(child)
        with open(file_path, 'w') as fh:  # noqa
            pass  # wipes old data
    except FileNotFoundError:
        pass
    finally:
        return child_threads


def katabasis(file_path, status_object):
    """ Used to save checkpoint data upon shutdown

    katabasis - a descent, often into the underworld (the system shutting down)

    :param file_path: Location of the save file
    :param status_object: Calling system's status object
    :return: None
    """
    try:
        with open(file_path, 'w+') as fh:
            temp = status_object.get_checkpoints(_all=True)
            json.dump(temp, fh, default=lambda o: o.__dict__())
    except:  # noqa
        pass
    # Unknown
