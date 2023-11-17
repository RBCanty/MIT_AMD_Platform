""" synchronization.py

Synchronizes the Faults between an Agent and the Master Controller

Uses a change-request format.  The MC holds the "official" version of the Faut list.  Agents submit reports to the MC,
and the MC records them.

Diagram::

| Fault update --> Local Agent --(if not MC)--report--> MC --> Official Fault Record
|                  V
|                  Local Fault Record
| or
| Fault update --> Agent --signal--> Local Agent --(if not MC)--report--> MC --> Official Fault Record
|                                    V
|                                    Local Fault Record

@author: Ben C
"""

from pprint import pformat

import database_interface as dbi
import mcn_status as mcs
import message as msgm
from constants import *
from mcn_logging_manager import system_log
from slacker import post_fault_to_slack


# Protect the private methods
__all__ = ["synchronized_fault_add", "synchronized_remove_fault"]


def synchronized_fault_add(data_pipes, new_fault: mcs.Fault):
    """ Report for adding a Fault

    :param data_pipes: Pipes to the local fault record
    :param new_fault: The Fault
    :return: None
    """
    system_log.info(pformat(new_fault))
    internals: mcs.Status = data_pipes[INTERNAL]
    if new_fault.location in mcs.s_dir(internals.name):
        _local_add_fault(data_pipes, new_fault)
    else:
        _remote_add_fault(data_pipes, new_fault)


def synchronized_remove_fault(data_pipes, cleared_fault: mcs.Fault):
    """ Report for removing a Fault

    :param data_pipes: Pipes to the local fault record
    :param cleared_fault: The Fault
    :return: None
    """
    system_log.info(pformat(cleared_fault))
    internals: mcs.Status = data_pipes[INTERNAL]
    if cleared_fault.location in mcs.s_dir(internals.name):
        _local_remove_fault(data_pipes, cleared_fault)
    else:
        _remote_remove_fault(data_pipes, cleared_fault)


def _local_add_fault(data_pipes, new_fault: mcs.Fault):
    """ (Private) Adds a fault at the local level

    :param data_pipes: Pipes to the local fault record
    :param new_fault: The Fault
    :return: None
    """
    internals: mcs.Status = data_pipes[INTERNAL]
    outbox = data_pipes[MSG_Q_OUT]
    my_name = internals.name
    # Catalogue Busy, Problem, and Fatal faults for data analytics
    file_locator = list()
    if new_fault.level in [mcs.V_BUSY, mcs.V_PROBLEM, mcs.V_FATAL]:
        catalogue_error = dbi.catalogue_fault(new_fault, file_locator)
        if catalogue_error:
            system_log.warning(f"Queue '{new_fault.queue}' was not catalogued: '{repr(catalogue_error)}'")
    # Only log Problem and Fatal faults to Queue doc (Busy Faults are not important here)
    if new_fault.level in [mcs.V_PROBLEM, mcs.V_FATAL]:
        record_error = dbi.record_fault(new_fault, file_locator)
        if record_error:
            system_log.warning(f"Queue '{new_fault.queue}' was not posted to DB: '{repr(record_error)}'")
    if new_fault not in internals.get_faults(_all=True):
        internals.add_fault(new_fault)
        post_fault_to_slack(new_fault)
    if my_name != "MC":
        outbox.enqueue(msgm.Message.status_report(internals.get_short_copy(), sender=my_name, recipient="MC"))


def _local_remove_fault(data_pipes, cleared_fault: mcs.Fault):
    """ (Private) Removes a fault at the local level

    :param data_pipes: Pipes to the local fault record
    :param new_fault: The Fault
    :return: None
    """
    internals: mcs.Status = data_pipes[INTERNAL]
    outbox = data_pipes[MSG_Q_OUT]
    my_name = internals.name
    internals.remove_fault(cleared_fault)
    if my_name != "MC":
        outbox.enqueue(msgm.Message.status_report(internals.get_short_copy(), sender=my_name, recipient="MC"))


def _remote_add_fault(data_pipes, new_fault: mcs.Fault):
    """ (Private) Adds a Fault remotely

    :param data_pipes: For accessing the outbox and for creating a Checkpoint on the addition of a new Fault to the
      affected system
    :param new_fault: The Fault
    :return: "timed out", "success", "no", "void"
    """
    internals: mcs.Status = data_pipes[INTERNAL]
    outbox = data_pipes[MSG_Q_OUT]
    my_name = internals.name
    my_msg = msgm.Message.add_fault(sender=my_name, system_with_fault=new_fault.location, fault=new_fault.dumps())
    idemp_key = my_msg[msgm.IDKEY]
    my_msg[msgm.PRIORITY] = msgm.P_CMD_RR
    outbox.enqueue(my_msg)
    checkpoint = mcs.Checkpoint(completion=None,
                                location=new_fault.location,
                                level=None,
                                data=f"Trying to add Fault: {new_fault}",
                                queue="Internal",
                                step="Sync",
                                operation="Sync",
                                task_key=idemp_key)
    internals.add_checkpoint(idemp_key, checkpoint)

    checkpoint_update = checkpoint.wait(120.0)
    internals.release_checkpoint(idemp_key)
    if not checkpoint_update:
        system_log.warning(f"remote add fault operation timed out ({str(new_fault)})")
        return "timed out"  # ConfirmationTimeout
    if checkpoint.completion is True:
        system_log.info("remote add fault success")
        return "success"
    elif checkpoint.completion is False:
        system_log.warning(f"remote add fault failed ({str(new_fault)})")
        return "no"
    else:
        system_log.warning(f"remote add fault ignored ({str(new_fault)})")
        return "void"


def _remote_remove_fault(data_pipes, old_fault: mcs.Fault):
    """ (Private) Remove a Fault remotely

    :param data_pipes: For accessing the outbox and for creating a Checkpoint on the removal of an old Fault from the
      affected system
    :param new_fault: The Fault
    :return: "timed out", "success", "no", "void"
    """
    internals: mcs.Status = data_pipes[INTERNAL]
    outbox = data_pipes[MSG_Q_OUT]
    my_name = internals.name
    my_msg = msgm.Message.remove_fault(sender=my_name,
                                       system_removing_fault=old_fault.location,
                                       fault=old_fault.dumps())
    idemp_key = my_msg[msgm.IDKEY]
    my_msg[msgm.PRIORITY] = msgm.P_CMD_RR
    outbox.enqueue(my_msg)
    checkpoint = mcs.Checkpoint(completion=None,
                                location=old_fault.location,
                                level=None,
                                data=f"Trying to remove Fault: {old_fault}",
                                queue="Internal",
                                step="Sync",
                                operation="Sync",
                                task_key=idemp_key)
    internals.add_checkpoint(idemp_key, checkpoint)

    checkpoint_update = checkpoint.wait(60.0)
    internals.release_checkpoint(idemp_key)
    if not checkpoint_update:
        system_log.warning(f"remote remove fault operation timed out ({str(old_fault)})")
        return "timed out"  # ConfirmationTimeout
    if checkpoint.completion is True:
        system_log.info("remote remove fault success")
        return "success"
    elif checkpoint.completion is False:
        system_log.warning(f"remote remove fault failed ({str(old_fault)})")
        return "no"
    else:
        system_log.warning(f"remote remove fault ignored ({str(old_fault)})")
        return "void"
