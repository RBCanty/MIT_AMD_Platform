""" Functionality for the Storage Carousel (LPX)
@author: Ben C
"""

import time
from pprint import pformat

import serial

import custom_exceptions as cexc
import database_interface as dbi
import functionals
import mcn_status as mcs
from AMD_LPX_Control import LpxController, TRANSFER_STATION, LPX_UPDATING
from database_constants import *
from mcn_logging_manager import system_log
from ui_exception import LPXPopupDialog

TRANSFER_LOCATION = ['lpx', TRANSFER_STATION]
UPDATING_LOCATION = ['lpx', LPX_UPDATING]


@functionals.standard_exception
def lpx_request(*_, **kwargs):
    """ Moves a plate from the LPX to the transfer station

    :param kwargs: Standard kwargs from MCN control
    :return: a RetObj
    """
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    q_name, step_num, _, _ = functionals.unpack_db_kwargs(kwargs, d_pid="", d_a="")
    with child_com as cc:
        my_lpx_controller = cc['Ss']

    # STEP ONE: Get some information:
    try:
        queue_doc, _, _ = dbi.query_document('SP', COL_QUEUE, Q_NAME, q_name)
    except cexc.DatabaseRequestError:
        error_message = f"Failed to pull queue document '{q_name}' for details"
        system_log.exception(error_message)
        return mcs.RetObj.incomplete('Ss', mcs.V_PROBLEM, error_message)
    this_operation = queue_doc[Q_OPERATIONS_LIST][str(step_num)]
    nickname_used = this_operation.get(QOP_CONTAINER, None)
    if nickname_used is None:
        return mcs.RetObj.incomplete('Ss', mcs.V_PROBLEM, f"Request made for undefined plate ({q_name})")
    this_plate = queue_doc[Q_CONTAINERS].get(nickname_used, None)
    if this_plate is None:
        return mcs.RetObj.incomplete('Ss',
                                     mcs.V_PROBLEM,
                                     f"Request made for unrecognized plate ({q_name}: {nickname_used})")
    plate_name = this_plate[DBG_CONTAINER_NAME]
    plate_type = this_plate[DBG_CONTAINER_PLATETYPE]

    # STEP TWO: If no plate name, find one
    if plate_name is None:
        with my_lpx_controller as ctrl:
            while True:
                found_location = ctrl.get_last_slot_of_type_from_db(plate_type)
                if found_location:
                    hotel_num, position_num = found_location
                    print("DEBUG(lpx_request):", hotel_num, position_num)  # DEBUG
                    break
                else:
                    system_log.warning(f"LPX lacking a plate of type {plate_type}, requesting reload form User")
                    loading_attempted = ctrl.load_onto_carousel()
                    if not loading_attempted:
                        system_log.info(f"User was unable to reload plates")
                        return mcs.RetObj.incomplete('Ss', mcs.V_FATAL, f"No plates of '{plate_type}' on LPX")
                time.sleep(1)
        try:
            doc, _, _ = dbi.query_location("SP", location_tag='lpx', sublocation=[str(hotel_num), int(position_num)])
        except cexc.DatabaseRequestError:
            system_log.warning(f"Database Request Error when looking up ['lpx', ['{hotel_num}', {position_num}]]")
            plate_name = None
        else:
            for _match in doc.get('full_matches', {}).values():
                plate_name = _match[dbi.DBG_CONTAINER_NAME]
                break
            else:
                system_log.warning(f"Database returned no plate in ['lpx', ['{hotel_num}', {position_num}]]"
                                   f" despite the plate there being promised")
        if plate_name is None:
            system_log.info(f"Either DB could not be accessed or "
                            f"no plates were found of type '{plate_type}' on the LPX")
            return mcs.RetObj.incomplete('Ss', mcs.V_FATAL, f"No plates of '{plate_type}' on LPX")
        try:
            dbi.update_field('SP', f"{Q_CONTAINERS}.{nickname_used}.{DBG_CONTAINER_NAME}", q_name, None, plate_name)
        except cexc.DatabaseRequestError:
            return mcs.RetObj.incomplete('Ss', mcs.V_PROBLEM, f"Failed to update plate name in DB (queue, ref, name):\n"
                                                              f"({q_name}, {nickname_used}, {plate_name})")
    else:
        try:
            deetz, _ = dbi.get_plate_details(plate_name, location_restriction='lpx')
        except TypeError:
            system_log.info(f"Plate '{plate_name}' location data unreachable")
            return mcs.RetObj.incomplete('Ss', mcs.V_FATAL, f"No plate location for '{plate_name}' found in DB")
        if not deetz:
            system_log.info(f"Plate '{plate_name}' either not found in DB or contains no information")
            return mcs.RetObj.incomplete('Ss', mcs.V_FATAL, f"No plate details for '{plate_name}' found in DB")

        _, [hotel_num, position_num] = deetz.get(DBG_LOCATION, [None, [None, None]])
        if hotel_num is None or position_num is None:
            system_log.info(f"Plate '{plate_name}' location data unreachable")
            return mcs.RetObj.incomplete('Ss', mcs.V_FATAL, f"No plate location for '{plate_name}' found in DB")

    # STEP THREE: Fetch the plate
    with my_lpx_controller as ctrl:
        # Check the Transfer Station
        ts_code = ctrl.is_transfer_station_occupied()  # (-1 error; 0 clear; 1 occupied)
        if ts_code == 1:
            system_log.info(f"LPX could not complete ({plate_name}, {plate_type}) request: Transfer Station occupied")
            return mcs.RetObj.incomplete("Ss", mcs.V_PROBLEM, "Transfer station occupied")
        elif ts_code == 0:
            pass
        else:
            system_log.warning(f"LPX encountered exception during operation ({plate_name}, {plate_type})")
            return mcs.RetObj.incomplete("Ss", mcs.V_PROBLEM, "Ss Exception")
        # Check the Database
        db_code = ctrl.check_location(*TRANSFER_STATION)
        if db_code:
            system_log.warning(f"Conflict on ({plate_name}, {plate_type}) request: "
                               f"DB reports Transfer Station occupied but LPX detects no plate there")
            return mcs.RetObj.incomplete("Ss", mcs.V_PROBLEM, "DB-LPX Conflict: Occupation of Transfer Station")

        try:
            sequence = ctrl.generate_transfer_sequence(hotel_num, position_num, -1)
        except LPXPopupDialog as popup:
            system_log.exception("Plate source location is not locatable, requesting human recovery")
            if not popup.run(plate_name, plate_type, -1):
                raise RuntimeError(f"Failed moving plate {plate_name} ({plate_type})")
            sequence = []

        manifest = list()
        try:
            ctrl.update_db_location(plate_name, UPDATING_LOCATION, source=['lpx', [hotel_num, position_num]])
        except cexc.DatabaseRequestError as dre:
            manifest.append(["DB-Pre", "Make DB Location Transient", repr(dre)])
            system_log.exception(f"Failed to make DB location transient:\n{pformat(manifest)}")
            return mcs.RetObj.incomplete("Ss", mcs.V_PROBLEM, "Failed to update DB location")
        for index, step in enumerate(sequence):
            result = _update_manifest(manifest, index, step, ctrl.execute(step))
            if result:
                system_log.info(f"LPX failed step ({index}, {step}) [{result}]:\n{pformat(manifest)}")
                return mcs.RetObj.incomplete("Ss", mcs.V_PROBLEM, f"LPX failed step ({index}, {step})")

        ts_code = ctrl.is_transfer_station_occupied()  # (-1 error; 0 clear; 1 occupied)
        if ts_code == 1:
            pass
        elif ts_code == 0:
            system_log.info(f"LPX could not complete ({plate_name}, {plate_type}) request: "
                            f"Transfer Station missing plate")
            return mcs.RetObj.incomplete("Ss", mcs.V_PROBLEM, "Transfer station missing plate post request")
        else:
            system_log.warning(f"LPX encountered exception during operation ({plate_name}, {plate_type})")
            return mcs.RetObj.incomplete("Ss", mcs.V_PROBLEM, "Ss Exception")

        try:
            ctrl.update_db_location(plate_name, TRANSFER_LOCATION)
        except cexc.DatabaseRequestError as dre:
            manifest.append(["DB-Post", "Make DB Location Real", repr(dre)])
            system_log.exception(f"LPX failed to make DB location real:\n{pformat(manifest)}")
            return mcs.RetObj.incomplete("Ss", mcs.V_PROBLEM, "Failed to update DB location")

    return mcs.RetObj.complete(manifest)


@functionals.standard_exception
def lpx_stow(*_, **kwargs):
    """ Moves a plate from the transfer station into the LPX

    :param kwargs: Standard kwargs from MCN control
    :return: a RetObj
    """
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    q_name, step_num, _, _ = functionals.unpack_db_kwargs(kwargs, d_pid="", d_a="")

    plate_name = kwargs[QOP_DETAILS].get('plate_name', None)
    if plate_name is None:
        try:
            queue_doc, _, _ = dbi.query_document('SP', COL_QUEUE, Q_NAME, q_name)
        except cexc.DatabaseRequestError:
            error_message = f"Failed to pull queue document '{q_name}' for details"
            system_log.exception(error_message)
            return mcs.RetObj.incomplete('Ss', mcs.V_PROBLEM, error_message)
        nickname = queue_doc[Q_OPERATIONS_LIST][str(step_num)][QOP_CONTAINER]
        plate_name = queue_doc[Q_CONTAINERS][nickname][DBG_CONTAINER_NAME]
        plate_type = queue_doc[Q_CONTAINERS][nickname].get(DBG_CONTAINER_PLATETYPE,
                                                           kwargs[QOP_DETAILS].get(DBG_CONTAINER_PLATETYPE, None))
        plate_oloc = queue_doc[Q_CONTAINERS][nickname].get(DBG_LOCATION,
                                                           kwargs[QOP_DETAILS].get(DBG_LOCATION,
                                                                                   TRANSFER_LOCATION)
                                                           )
    else:
        plate_type = kwargs[QOP_DETAILS].get('plate_type', None)
        plate_oloc = kwargs[QOP_DETAILS].get(DBG_LOCATION, TRANSFER_LOCATION)

    if plate_type is None:
        deetz, _ = dbi.get_plate_details(plate_name, location_restriction='lpx')
        if deetz is None:
            system_log.info(f"No plate '{plate_name}' found in DB")
            return mcs.RetObj.incomplete('Ss', mcs.V_FATAL, f"No plate '{plate_name}' found in DB")
        plate_type = deetz.get(DBG_LABWARE_TYPE, None)
        if plate_type is None:
            system_log.info(f"No plate '{plate_name}' found in DB")
            return mcs.RetObj.incomplete('Ss', mcs.V_FATAL, f"No plate '{plate_name}' found in DB")
        plate_oloc = deetz.get(DBG_LOCATION, TRANSFER_LOCATION)

    with child_com as cc:
        my_lpx_controller: LpxController = cc['Ss']
    with my_lpx_controller as ctrl:
        # Check the Transfer Station
        ts_code = ctrl.is_transfer_station_occupied()  # (-1 error; 0 clear; 1 occupied)
        if ts_code == 1:
            pass
        elif ts_code == 0:
            system_log.info(f"LPX could not complete ({plate_name}, {plate_type}) stow: No plate in Transfer Station")
            return mcs.RetObj.incomplete("Ss", mcs.V_PROBLEM, "No plate in Transfer Station")
        else:
            system_log.warning(f"LPX encountered exception during operation ({plate_name}, {plate_type})")
            return mcs.RetObj.incomplete("Ss", mcs.V_PROBLEM, "Ss Exception")
        # Check the Database
        db_code = ctrl.check_location(*TRANSFER_STATION)
        if not db_code:
            system_log.warning(f"Conflict on ({plate_name}, {plate_type}) stow: "
                               f"LPX reports Transfer Station occupied but DB reports no plate there")
            return mcs.RetObj.incomplete("Ss", mcs.V_PROBLEM, "DB-LPX Conflict: Occupation of Transfer Station")

        manifest = list()
        while True:
            try:
                found_location = ctrl.get_first_available_slot_from_db(plate_type)
                if found_location:
                    hotel_num, position_num = found_location
                    break
                else:
                    _reorganize(ctrl, manifest, True)
                    # This means we've checked all plates in the LPX and none could be moved to make a different
                    system_log.warning(f"LPX lacking a space for a {plate_type}, requesting unload form User")
                    ctrl.load_onto_carousel()
            except cexc.BreakOuterLoop:
                pass
            finally:
                time.sleep(1)  # There are a lot of DB accesses here, so a bit of a cool-down

        # One way or another, there is now a place for the plate
        sequence = ctrl.generate_transfer_sequence(hotel_num, position_num, 1)
        try:
            ctrl.update_db_location(plate_name, UPDATING_LOCATION, source=plate_oloc)
        except cexc.DatabaseRequestError as dre:
            manifest.append(["DB-Pre", "Make DB Location Transient", repr(dre)])
            system_log.exception(f"Failed to make DB location transient:\n{pformat(manifest)}")
            return mcs.RetObj.incomplete("Ss", mcs.V_PROBLEM, "Failed to update DB location")
        for index, step in enumerate(sequence):
            result = _update_manifest(manifest, index, step, ctrl.execute(step))
            if result:
                system_log.info(f"LPX failed step ({index}, {step}) [{result}]:\n{pformat(manifest)}")
                return mcs.RetObj.incomplete("Ss", mcs.V_PROBLEM, f"LPX failed step ({index}, {step})")

        ts_code = ctrl.is_transfer_station_occupied()  # (-1 error; 0 clear; 1 occupied)
        if ts_code == 1:
            system_log.info(f"LPX could not complete ({plate_name}, {plate_type}) stow: "
                            f"Plate still in Transfer Station")
            return mcs.RetObj.incomplete("Ss", mcs.V_PROBLEM, "Plate still in Transfer Station")
        elif ts_code == 0:
            pass
        else:
            system_log.warning(f"LPX encountered exception during operation ({plate_name}, {plate_type})")
            return mcs.RetObj.incomplete("Ss", mcs.V_PROBLEM, "Ss Exception")

        try:
            ctrl.update_db_location(plate_name, ['lpx', [hotel_num, position_num]])
        except cexc.DatabaseRequestError as dre:
            manifest.append(["DB-Post", "Make DB Location Real", repr(dre)])
            system_log.exception(f"LPX failed to make DB location real:\n{pformat(manifest)}")
            return mcs.RetObj.incomplete("Ss", mcs.V_PROBLEM, "Failed to update DB location")

    return mcs.RetObj.complete(manifest)


def lpx_direct(*_, **kwargs):
    """ For issuing a direct, low-level command to the LPX

    :param kwargs: Standard kwargs from MCN control
    :keyword step: A string passed to the LPX execute method (default is 'RD 1915')
    :return: a RetObj
    """
    step = kwargs.get('step', 'RD 1915')
    step = step + chr(13)

    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    with child_com as cc:
        my_lpx_controller = cc['Ss']
    with my_lpx_controller as ctrl:
        ret_val = ctrl.execute(step)

    return mcs.RetObj.complete(f"'{step}' returned '{ret_val}'")


def open_lpx_ui(*_, **kwargs):
    """ Spawns an instance of the plate storage GUI

    :param kwargs: Standard kwargs from MCN control
    :return: a RetObj
    """
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    with child_com as cc:
        my_lpx_controller: LpxController = cc['Ss']
    with my_lpx_controller as ctrl:
        ctrl.load_onto_carousel()

    return mcs.RetObj.complete("LPX GUI terminated")


@functionals.standard_exception
def lpx_locations(*_, **kwargs):
    """ Debug: provides the last occupied and fist available slots for a given plate type

    :param kwargs: Standard kwargs from MCN control
    :keyword req_plate_type: The plate type being searched (default: "96 Well DeepWell")
    :return: a RetObj with the requested info
    """
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    plate_type = kwargs.get("req_plate_type", "96 Well DeepWell")
    with child_com as cc:
        my_lpx_controller: LpxController = cc['Ss']
    with my_lpx_controller as ctrl:
        last = ctrl.get_last_slot_of_type_from_db(plate_type)
        first = ctrl.get_first_available_slot_from_db(plate_type)

    return mcs.RetObj.complete(f"For a plate_type '{plate_type}', pull from {last} and place onto {first}")


@functionals.standard_exception
def reorganize(*_, **kwargs):
    """ Moves plates around internally to make more efficient use of space

    :param kwargs: Standard kwargs from MCN control
    :return: a RetObj
    """
    manifest = list()

    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    with child_com as cc:
        my_lpx_controller = cc['Ss']

    with my_lpx_controller as ctrl:
        try:
            plate_types = ctrl.configuration['wellplates'] + ctrl.configuration['consumables']

            for plate_type in plate_types:
                found_location = ctrl.get_first_available_slot_from_db(plate_type)
                if found_location:
                    continue
                _reorganize(ctrl, manifest)
        except:  # noqa
            system_log.exception(f"Ss encountered an exception while trying to reorganize itself:\n"
                                 f"{pformat(manifest)}")
            return mcs.RetObj.incomplete("Ss", mcs.V_PROBLEM, "LPX Reogranization failed")

    return mcs.RetObj.complete(manifest)


@functionals.standard_exception
def hard_stop(*_, **kwargs):
    """ Calls the shutdown method for the LPX controller

    :param kwargs: Standard kwargs from MCN control
    :return: a RetObj
    """
    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    with child_com as cc:
        my_lpx_controller = cc['Ss']
    with my_lpx_controller as ctrl:
        manifest = ctrl.emergency_shutdown()

    return mcs.RetObj.complete(manifest)


@functionals.standard_exception
def reconnect_lpx(*_, **kwargs):
    """ Attempts to reconnec the LPX

    :param kwargs: Standard kwargs from MCN control
    :keyword port: The com port of the LPX (default: "COM10")
    :return: a RetObj
    """
    port = kwargs.get("port", "COM10")

    _, _, child_com = functionals.unpack_kwargs(kwargs, d_msg="", d_cso="")
    with child_com as cc:
        my_lpx_controller = cc['Ss']
    with my_lpx_controller as ctrl:
        try:
            ctrl.lpx_port.close()
        except:  # noqa
            pass
        ctrl.lpx_serial_settings = {'port': port, 'baudrate': 9600, 'bytesize': 8,
                                    'stopbits': serial.STOPBITS_ONE, 'parity': serial.PARITY_EVEN}
        ctrl.lpx_port = serial.Serial(**ctrl.lpx_serial_settings)
        ctrl.lpx_port.write(bytes("CR\r", 'Ascii'))
        time.sleep(1)
        while ctrl.lpx_port.in_waiting > 0:
            ctrl.last_response = "".join(ctrl.lpx_port.readline().decode('Ascii').split())
        ret_val = ctrl.last_response
        system_log.info(f"LPX reconnect received: {ret_val}")
        ctrl.initialize()

    return mcs.RetObj.complete(ret_val)


def _reorganize(ctrl, manifest, raises=False):
    # No locations found, so try to fix it yourself
    partial_matches = ctrl.get_plates_from_db('partial_matches')
    # Builds a reference for moving plates around
    plate_types = ctrl.configuration['wellplates'] + ctrl.configuration['consumables']
    limits = {_type: min(ctrl.configuration[_type]) for _type in plate_types}
    # Check all plates in LPX
    for plate_id in partial_matches:
        # Get a plate's type, current location, & the hotels that it shouldn't be in if it can be helped
        c_type = partial_matches[plate_id]['labware_type']
        c_hotel, c_position = partial_matches[plate_id]['location'][1]
        c_limit = min([x for x in limits if x > limits[c_type]])
        plate_name = partial_matches[plate_id]['container_name']
        # If it is in a hotel of a less restrictive plate type, move it to a more restrictive, compatible plate type
        if c_hotel > c_limit:
            # What is the most restrictive hotel this plate could fit on
            found_relocation = ctrl.get_first_available_slot_from_db(c_type)
            if not (found_relocation is None):
                n_hotel, n_position = found_relocation
                # Only bother if the plate would change hotel in the move
                if n_hotel < c_hotel:
                    try:
                        # Try to generate the full sequence...
                        sequence = ctrl.generate_transfer_sequence(c_hotel, c_position, 0)
                    except LPXPopupDialog as popup:
                        system_log.exception("Plate source location is not locatable, requesting human recovery")
                        # ...if it fails, then ask the user to do first part (moving to transfer station)...
                        if not popup.run(plate_name, c_type, -1):
                            raise RuntimeError(f"Failed moving plate {plate_name} ({c_type})")
                        # ...then automatically do the second part (stowing back on LPX)
                        sequence = ctrl.generate_transfer_sequence(n_hotel, n_position, 1)
                    manifest.append([f"Attempting to relocate a {c_type} to make room",
                                     c_hotel, c_position])
                    try:
                        ctrl.update_db_location(plate_name, UPDATING_LOCATION, source=['lpx', [c_hotel, c_position]])
                    except cexc.DatabaseRequestError as dre:
                        manifest.append(["DB-Pre", "Make DB Location Transient", repr(dre)])
                        system_log.exception(f"LPX exception when making location transient:\n"
                                             f"{pformat(manifest)}")
                        return mcs.RetObj.incomplete("Ss",
                                                     mcs.V_PROBLEM,
                                                     "LPX exception when making location transient")
                    for index, step in enumerate(sequence):
                        result = _update_manifest(manifest, index, step, ctrl.execute(step))
                        if result:
                            system_log.info(f"LPX failed on ({index}, {step}) [{result}]:\n"
                                            f"{pformat(manifest)}")
                            return mcs.RetObj.incomplete("Ss",
                                                         mcs.V_PROBLEM,
                                                         f"LPX failed on ({index}, {step}) [{result}]")
                    manifest.append([f"Relocate a {c_type}.", n_hotel, n_position])
                    try:
                        ctrl.update_db_location(plate_name, ['lpx', [n_hotel, n_position]])
                    except cexc.DatabaseRequestError as dre:
                        manifest.append(["DB-Post", "Make DB Location Real", repr(dre)])
                        system_log.exception(f"LPX failed to make location real in DB:\n"
                                             f"{pformat(manifest)}")
                        return mcs.RetObj.incomplete("Ss", mcs.V_PROBLEM, "LPX failed to update DB")
                    if raises:
                        raise cexc.BreakOuterLoop(memory=manifest)
    return manifest


def _update_manifest(manifest, index, step, result):
    if result == 0:
        manifest.append([index, str(step), "Success"])
    elif (step == b'RD 1915\r') and (result == 1):
        manifest.append([index, str(step), "Success"])
    elif result == 1:
        manifest.append([index, str(step), "Timeout"])
        return f"LPX transfer step {step} failed due to timeout"
    elif result == 2:
        manifest.append([index, str(step), "Fault"])
        return f"LPX transfer step {step} failed due to fault"
    else:
        system_log.warning(f"A step in the transfer sequence ({index}, {step}) was Forced, true state unknown.")
        manifest.append([index, str(step), "Forced"])
    return ""
