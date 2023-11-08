""" Functionality of the thermoreactor
@authors: Ben C, Brent K
"""

from pprint import pformat

import custom_exceptions as cexc
import database_interface as dbi
import functionals
import mcn_status as mcs
from constants import *
from database_constants import *
from mcn_logging_manager import system_log
from operations import Operation
from ui_thermoreactor import ThermoreactorGUI
import tkinter as tk


@functionals.standard_exception
def thermo_prepare(*_, **kwargs):
    """ Prepares thermoreactor for a transfer

    :param kwargs: Standard kwargs from MCN control
    :return: a RetObj
    """
    queue_name, operation_number, _, _ = functionals.unpack_db_kwargs(kwargs, d_pid="", d_a="")
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    with child_com as cc:
        thermo_controller = cc['Th']
    with thermo_controller as my_controller:
        response = my_controller.prepare_for_transfer()
    if len(response) == 1:
        response += ["(Status queue empty)"]
    quick_code = response[0]
    if quick_code == 'success':
        pass
    elif quick_code == 'busy':
        system_log.info(f"Th reporting busy:\n{pformat(response)}")
        return mcs.RetObj.incomplete("Th", mcs.V_BUSY, "Th reporting busy")
    elif quick_code == 'error':
        system_log.info(f"Th reporting error:\n{pformat(response)}")
        return mcs.RetObj.incomplete("Th", mcs.V_FATAL, "Th reporting error")
    else:
        system_log.info(f"Th reporting unknown state:\n{pformat(response)}")
        return mcs.RetObj.incomplete("Th", mcs.V_FATAL, "Th state unknown")

    try:
        doc, err_details, resp_code = dbi.query_document('MC', COL_CONSUMABLES, CON_NAME, 'PFA_films')
    except cexc.DatabaseRequestError:
        system_log.exception(f"Th DB request failed\n"
                             f"RECORD: Th status:\n"
                             f"{pformat(response)}")
        return mcs.RetObj.incomplete('SP', mcs.V_FATAL, "Th DB request failed")

    film_state = doc['location'][1]
    if film_state == 'occupied':
        _func = 'unload_pfa_film'
    elif film_state == 'empty':
        _func = 'load_pfa_film'
    else:
        system_log.info(f"Th.thermo_prepare failed due to invalid PFA film state: {film_state}\n"
                        f"RECORD: Th status:\n"
                        f"{pformat(response)}")
        return mcs.RetObj.incomplete('Th', mcs.V_PROBLEM, f"Invalid pfa film state: {film_state}")

    with thermo_controller as my_controller:
        response = my_controller.prepare_for_film_transfer()

    try:
        ret_obj = functionals.proxy_run(_pipes={MSG_Q_IN: kwargs[MSG_Q_IN],
                                                MSG_Q_OUT: kwargs[MSG_Q_OUT],
                                                CHILD_COM: kwargs[CHILD_COM],
                                                INTERNAL: kwargs[INTERNAL]},
                                        _name="SP",
                                        _recip="MC",
                                        _type="ib.Run_Local",
                                        _operation=Operation.build_from_args(func=_func, agent="Ra"),
                                        _wait=20,
                                        _timeout=5*60
                                        )
    except cexc.ConfirmationTimeout:
        system_log.warning(f"Operation {_func} has timed out\n"
                           f"RECORD: Th status:\n"
                           f"{pformat(response)}")
        return mcs.RetObj.incomplete('SP', mcs.V_PROBLEM, f"{_func} has timed out")

    if not ret_obj.completion:
        ret_obj.data = pformat({"Thermo Response": response,
                                "Error Details": ret_obj.data})
        return ret_obj

    with thermo_controller as my_controller:
        response = my_controller.prepare_for_transfer()

    return mcs.RetObj.complete(response)


@functionals.standard_exception
def thermo_run(*_, **kwargs):
    """ Executes a thermoreactor tun

    :param kwargs: Standard kwargs from MCN control
    :return: a RetObj
    """
    details = kwargs[QOP_DETAILS]
    queue_name, operation_number, _, _ = functionals.unpack_db_kwargs(kwargs, d_pid="", d_a="")
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    with child_com as cc:
        thermo_controller = cc['Th']
    with thermo_controller as my_controller:
        response = my_controller.thermo_reactor_profile_selector({'heating_profile': details['heating_profile'],
                                                                  'queue_name': queue_name})
    if response[0] == 'Finished operation':
        pass
    elif response[0] == 'busy':
        system_log.info(f"Th is busy: {pformat(my_controller.check_for_status())}")
        return mcs.RetObj.incomplete('Th', mcs.V_BUSY, "Th busy")
    elif response == 'Thermoreactor encountered an error':
        system_log.warning(f"Th encountered an error: {pformat(my_controller.check_for_status())}")
        return mcs.RetObj.incomplete('Th', mcs.V_FATAL, "Th error")
    else:
        system_log.warning(f"Th reporting unknown status: {pformat(response)}")
        return mcs.RetObj.incomplete('Th', mcs.V_FATAL, "Th state unknown")

    return mcs.RetObj.complete(response)


@functionals.standard_exception
def thermo_check(*_, **kwargs):
    """ Inquired for the state of the thermoreactor

    :param kwargs: Standard kwargs from MCN control
    :return: a RetObj (complete if 'idle', incomplete otherwise)
    """
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    with child_com as cc:
        thermo_controller = cc['Th']
    with thermo_controller as my_controller:
        response = my_controller.check_for_status()

    if response[0] == 'idle':
        # return [GOOD, response[1]]
        return mcs.RetObj.complete(response)
    elif response[0] == 'busy':
        # return [BUSY, response[1]]
        system_log.info(f"Th reporting busy:\n{pformat(response)}")
        return mcs.RetObj.incomplete("Th", mcs.V_BUSY, "Th reporting busy")
    elif response[0] == 'error':
        # return [FATAL, response[1]]
        system_log.info(f"Th reporting error:\n{pformat(response)}")
        return mcs.RetObj.incomplete("Th", mcs.V_FATAL, "Th reporting error")
    else:
        # return [PROBLEM, response[1]]
        system_log.info(f"Th reporting unknown state:\n{pformat(response)}")
        return mcs.RetObj.incomplete("Th", mcs.V_PROBLEM, "Th state unknown")


@functionals.standard_exception
def thermo_initial_state(*_, **kwargs):
    """ Attempts to move the thermoreactor into a default state

    :param kwargs: Standard kwargs from MCN control
    :return: a RetObj
    """
    queue_name, operation_number, _, _ = functionals.unpack_db_kwargs(kwargs, d_pid="", d_a="")
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    with child_com as cc:
        thermo_controller = cc['Th']
    with thermo_controller as my_controller:
        response = my_controller.go_to_initial_state()
    if len(response) == 1:
        response += ["(Status queue empty)"]
    quick_code = response[0]
    if quick_code == 'success':
        pass
    elif quick_code == 'error':
        system_log.info(f"Th reporting error:\n{pformat(response)}")
        return mcs.RetObj.incomplete('Th', mcs.V_FATAL, "Th reporting error")
    elif quick_code == 'busy':
        system_log.info(f"Th reporting busy:\n{pformat(response)}")
        return mcs.RetObj.incomplete('Th', mcs.V_BUSY, "Th reporting busy")
    else:
        system_log.info(f"Th reporting unknown state:\n{pformat(response)}")
        return mcs.RetObj.incomplete('Th', mcs.V_FATAL, "Th state unknown")

    return mcs.RetObj.complete(response)


@functionals.standard_exception
def thermo_door_open(*_, **kwargs):
    """ Opens door

    :param kwargs: Standard kwargs from MCN control
    :return: a RetObj
    """
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    with child_com as cc:
        thermo_controller = cc['Th']
    with thermo_controller as my_controller:
        response = my_controller.pcb_controller.open_door()
    system_log.debug(str(response))
    return mcs.RetObj.complete(response)


@functionals.standard_exception
def thermo_door_close(*_, **kwargs):
    """ closes door

    :param kwargs: Standard kwargs from MCN control
    :return: a RetObj
    """
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    with child_com as cc:
        thermo_controller = cc['Th']
    with thermo_controller as my_controller:
        response = my_controller.pcb_controller.close_door()
    system_log.debug(str(response))
    return mcs.RetObj.complete(response)


@functionals.standard_exception
def thermo_piston_up(*_, **kwargs):
    """ Moves the piston up (release)

    :param kwargs: Standard kwargs from MCN control
    :return: a RetObj
    """
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    with child_com as cc:
        thermo_controller = cc['Th']
    with thermo_controller as my_controller:
        response = my_controller.pcb_controller.release()
    system_log.debug(str(response))
    return mcs.RetObj.complete(response)


@functionals.standard_exception
def thermo_piston_down(*_, **kwargs):
    """ Moves the piston down (seal)

    :param kwargs: Standard kwargs from MCN control
    :return: a RetObj
    """
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    with child_com as cc:
        thermo_controller = cc['Th']
    with thermo_controller as my_controller:
        response = my_controller.pcb_controller.press()
    system_log.debug(str(response))
    return mcs.RetObj.complete(response)


@functionals.standard_exception
def pfa_film_load(*_, **kwargs):
    """ Loads a PFA film into the reactor

    :param kwargs: Standard kwargs from MCN control
    :return: a RetObj
    """
    _func = 'load_pfa_film'
    try:
        return functionals.proxy_run(_pipes={MSG_Q_IN: kwargs[MSG_Q_IN],
                                             MSG_Q_OUT: kwargs[MSG_Q_OUT],
                                             CHILD_COM: kwargs[CHILD_COM],
                                             INTERNAL: kwargs[INTERNAL]},
                                     _name="SP",
                                     _recip="MC",
                                     _type="ib.Run_Local",
                                     _operation=Operation.build_from_args(func=_func, agent="Ra"),
                                     _wait=20,
                                     _timeout=5 * 60
                                     )
    except cexc.ConfirmationTimeout:
        system_log.exception("PFA Film Loading Timeout")
        return mcs.RetObj.incomplete('SP', mcs.V_PROBLEM, f"Error Details: {_func} has timed out")


@functionals.standard_exception
def pfa_film_unload(*_, **kwargs):
    """ Removes a PFA film from the reactor

    :param kwargs: Standard kwargs from MCN control
    :return: a RetObj
    """
    _func = 'unload_pfa_film'
    try:
        return functionals.proxy_run(_pipes={MSG_Q_IN: kwargs[MSG_Q_IN],
                                             MSG_Q_OUT: kwargs[MSG_Q_OUT],
                                             CHILD_COM: kwargs[CHILD_COM],
                                             INTERNAL: kwargs[INTERNAL]},
                                     _name="SP",
                                     _recip="MC",
                                     _type="ib.Run_Local",
                                     _operation=Operation.build_from_args(func=_func, agent="Ra"),
                                     _wait=20,
                                     _timeout=5 * 60
                                     )
    except cexc.ConfirmationTimeout:
        system_log.exception("PFA Film Unloading Timeout")
        return mcs.RetObj.incomplete('SP', mcs.V_PROBLEM, f"Error Details: {_func} has timed out")


@functionals.standard_exception
def run_thermoreactor_gui(*_, **kwargs):
    """ Spawns an instance of the thermoreactor GUI

    :param kwargs: Standard kwargs from MCN control
    :return: a RetObj
    """
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    with child_com as cc:
        thermo_controller = cc['Th']

    thermo_gui_core = tk.Tk()
    ThermoreactorGUI(thermo_gui_core, thermo_controller)
    thermo_gui_core.mainloop()

    return mcs.RetObj.complete("Thermoreactor GUI terminated")
