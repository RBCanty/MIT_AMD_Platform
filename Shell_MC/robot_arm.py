#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 13:34:23 2021

@author: Matt McDonald
@editor: Ben Canty
"""

import socket
from copy import deepcopy
from threading import Lock

import custom_exceptions as cexc
import database_interface as dbi
import mcn_status as mcs
from database_constants import *
from mcn_logging_manager import system_log

ENCODING = 'Ascii'
INIT_SIGNAL = b'Initialize'
STOP_SIGNAL = b'Stop'
LOCATION_CODEX = {'liquid_handler': ([19, 5, 84], [19, 6, 84],), # [19, 7, 84], [19, 8, 84]), #temp access to 7, 8
                  'liquid_handler_other': ([24, 1, 329],), #, [28, 2, 329]),  # for non-96 well normal plates
                  'fraction_collector': ([1], [2], [3]),
                  'thermal_reactor': ('thermal_reactor',),
                  'autosampler': ('autosampler',),
                  'ir_spectrometer': ('ir_spectrometer',),
                  'lpx': (['0', 1],)}
STD_PLATE_TYPES = ('96 Well Microplate', )  #,'96 Well Microplate Half Area'
SCARA_IP = "123.456.7.8"
SCARA_PORT = 8082
COM_ERROR = mcs.RetObj.incomplete("Ra", mcs.V_PROBLEM, "Robot COM error")


class RoboticArmController:
    """ Controller for the Robotic Arm """
    def __init__(self, child_com, ra_socket=None):
        """ Creates a controller for the Robotic Arm

        Note that Guidance Software should be running alongside this

        :param child_com: (Needed for Barcode Reading)
        :param ra_socket: Allows the IP and PORT to be specified, default values taken from *.py file
        """
        if ra_socket is None:
            ra_socket = [SCARA_IP, SCARA_PORT]
        self._lock = Lock()
        self.child_com = child_com
        self.make_socket = lambda: socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.sock_addr = tuple(ra_socket)
        self.holding_nmr_gripper = False

    def __enter__(self):
        self._lock.acquire()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type is not None:
            system_log.exception(f"Robotic Arm Controller:: "
                                 f"\n\t- type = {exc_type}"
                                 f"\n\t- value = {exc_val}",
                                 exc_info=exc_tb)
        self._lock.release()

    def acquire(self, *args, **kwargs):
        return self._lock.acquire(*args, **kwargs)

    def release(self):
        return self._lock.release()

    def __del__(self):
        try:
            system_log.info("Ra Shutdown called")
            self.stop()
        except:
            system_log.exception("Ra encountered exception during shutdown")

    def initialize(self):
        with self.make_socket() as scara_sock:
            try:
                scara_sock.connect(self.sock_addr)
                sent = scara_sock.sendall(INIT_SIGNAL)
                if sent == 0:
                    system_log.warning(f"Socket failed to transmit initialization signal: {sent}")
                    return COM_ERROR
                data = scara_sock.recv(1024).decode(ENCODING)
            except Exception:  # ConnectionRefusedError
                system_log.exception(f"Socket failed to transmit initialization signal")
                return mcs.RetObj.incomplete("Ra", mcs.V_PROBLEM, f"Please manually start Ra Guidance Program")
        if data == 'Initialized':
            system_log.info("Robot initialized")
            return mcs.RetObj.complete("Robot initialized")
        else:
            system_log.warning(f"Robot failed to initialize: {data}")
            return mcs.RetObj.incomplete("Ra", mcs.V_FATAL, "Robot failed to initialize")

    def move_plate(self, q_id, q_step, override: dict = None):
        """ Moves a wellplate via the robotic arm

        finds the wellplate source by searching for the wellplate in the plate database, then
        searches for available destination at the correct agent (agent is specified in the operation details)

        :raises DatabaseRequestError: when DB is not connected
        :raises ValueError: If override is improperly specified or if specified on top of run-by-db fields

        :param q_id: a queue name to be looked up
        :param q_step: the step in queue q_id
        :param override: Keys = [container_name, destination]
        :return: Return object
        """
        # Verify we're not running by DB and Override
        if override and any([q_id, q_step]):
            raise ValueError(f"Can't specify Run by Override and provide Run by DB arguments")

        # Load in data from override or database
        if override:
            # Like with DB lookup, we build the container name and destination first
            _container_name = override.get('container_name', None)
            _destination = override.get('destination', None)
            # ?? Maybe make container_name an object_name... to reload films, etc.?
            if not all([_container_name, _destination, ]):
                raise ValueError(f"If override is specified, it requires all arguments, "
                                 f"missing: {[k for k,v in override if v is None]}")
                
            # We then look up the plate type and source location
            try:
                _p_doc, _p_collection = dbi.get_plate_details(_container_name)
            except cexc.DatabaseRequestError:
                exc_msg = "Ra DB request failed - Plate details"
                system_log.exception(exc_msg)
                return mcs.RetObj.incomplete('Ra', mcs.V_PROBLEM, exc_msg)
            try:
                _source = _p_doc['location']
                if 'labware_type' in _p_doc.keys():
                    _plate_type = _p_doc['labware_type']
                elif 'plate_type' in _p_doc.keys():
                    _plate_type = _p_doc['plate_type']
                else:
                    return mcs.RetObj.incomplete('Ra', mcs.V_FATAL, 'plate/labware_type key error')
            except KeyError:
                exc_msg = "Ra DB request failed - likely plate not given container_name"
                system_log.exception(exc_msg)
                return mcs.RetObj.incomplete('Ra', mcs.V_PROBLEM, exc_msg)

            q_step = 'override'
        else:
            # For DB lookup, we use the queue name and step number to find the container name and destination
            # But we also get the plate type for free
            try:
                q_doc, q_err, rc = dbi.query_document('MC', 'queue', 'queue_name', q_id)
            except cexc.DatabaseRequestError:
                exc_msg = "Ra DB request failed - Destination location"
                system_log.exception(exc_msg)
                return mcs.RetObj.incomplete('Ra', mcs.V_PROBLEM, exc_msg)

            q_step = str(int(q_step))

            plate_name = q_doc[Q_OPERATIONS_LIST][q_step][QOP_CONTAINER]
            _container_name = q_doc[Q_CONTAINERS][plate_name]['container_name']
            _destination = q_doc[Q_OPERATIONS_LIST][q_step][QOP_DETAILS]['target_destination']
            if 'labware_type' in q_doc[Q_CONTAINERS][plate_name].keys():
                _plate_type = q_doc[Q_CONTAINERS][plate_name]['labware_type']
            elif 'plate_type' in q_doc[Q_CONTAINERS][plate_name].keys():
                _plate_type = q_doc[Q_CONTAINERS][plate_name]['plate_type']
            else:
                return mcs.RetObj.incomplete('Ra', mcs.V_FATAL, 'plate/labware_type key error')

            # lookup the plate details (search all collections)
            try:
                _p_doc, _p_collection = dbi.get_plate_details(_container_name)
            except cexc.DatabaseRequestError:
                exc_msg = "Ra DB request failed - Plate details"
                system_log.exception(exc_msg)
                return mcs.RetObj.incomplete('Ra', mcs.V_PROBLEM, exc_msg)
            try:
                _source = _p_doc['location']
            except KeyError:
                exc_msg = "Ra DB request failed - likely plate not given container_name"
                system_log.exception(exc_msg)
                return mcs.RetObj.incomplete('Ra', mcs.V_PROBLEM, exc_msg)

        # Either via Override or Queue-Step lookup, we now have:
        source = _source
        destination = _destination
        collection = _p_collection
        plate_type = _plate_type
        plate_document = _p_doc

        # Validate the source location
        valid, scara_source = self.location_translator(source)
        if not valid:
            system_log.info(f"Ra source invalid: {scara_source}")
            return mcs.RetObj.incomplete('Ra', mcs.V_PROBLEM, "Ra source invalid")

        # Find a sublocation at the destination
        if destination == 'liquid_handler':
            if plate_type not in STD_PLATE_TYPES:
                lib_dest = 'liquid_handler_other'
            else:
                lib_dest = 'liquid_handler'
        else:
            lib_dest = destination

        found_available_location = False
        for check_location in LOCATION_CODEX[lib_dest]:
            try:
                location_return, _, _ = dbi.query_location('MC', destination, check_location)
            except cexc.DatabaseRequestError:
                exc_msg = "Ra DB request failed - plate location lookup"
                system_log.exception(exc_msg)
                continue

            if not location_return['full_matches']:
                lib_dest = [lib_dest, check_location]
                destination = [destination, check_location]
                found_available_location = True
                break

        if not found_available_location:
            system_log.info(f"No available locations at {destination} for Ra to complete move_wellplate")
            return mcs.RetObj.incomplete('Ra', mcs.V_BUSY, 'No available locations')

        system_log.debug(f"Found location:\n{lib_dest}")

        # Validate the Destination
        valid, scara_destination = self.location_translator(destination)
        if not valid:
            system_log.info(f"Ra destination invalid: {scara_destination}")
            return mcs.RetObj.incomplete('Ra', mcs.V_PROBLEM, "Ra destination invalid")

        command_string = scara_source + " to " + scara_destination
        with self.make_socket() as scara_sock:
            scara_sock.connect(self.sock_addr)
            sent = scara_sock.sendall(command_string.encode(ENCODING))
            if sent == 0:
                system_log.warning(f"Socket failed to transmit signal: {sent}")
                return COM_ERROR
            data = scara_sock.recv(1024).decode(ENCODING)
        if data[0] != 'S':
            system_log.info(f"Robot arm reports a failed transfer: {data}")
            return mcs.RetObj.incomplete("Ra", "Problem", "Ra plate move failed")

        # If the function has made it to this point we can go ahead and update the library location
        updated_container = deepcopy(plate_document)
        updated_container['location'] = destination
        try:
            req_ret_0, req_ret_1, req_code = dbi.update_document('MC', collection, plate_document, updated_container)
            system_log.debug(f"Robotic Arm's request ({source} -> {destination}) "
                             f"to update a document was: {req_ret_0}, {req_ret_1}, {req_code}")
        except cexc.DatabaseRequestError:
            exc_msg = "Ra DB request failed - plate location update"
            system_log.exception(exc_msg)
            return mcs.RetObj.incomplete('Ra', mcs.V_PROBLEM, exc_msg)

        return mcs.RetObj.complete(data)

    def move_pfa_film(self, load_clean: bool = True):
        """ Move a pfa film from clean storage (load_clean == True) of move a 
        used PFA film to dirty
        """
        PFA_FILM_LOCATIONS = ['0', '1', '2', '3']  # not sure if this is an exhaustive list
        try:
            film_doc, film_details, rc = dbi.query_document('MC', 'consumables', 'consumable_name', 'PFA_films')
        except cexc.DatabaseRequestError:
            exc_msg = "Ra DB request failed - PFA film location"
            system_log.exception(exc_msg)
            return  mcs.RetObj.incomplete('Ra', mcs.V_PROBLEM, exc_msg)
        updated_films = deepcopy(film_doc)
        if load_clean:
            if film_doc['location'][1] == 'occupied':
                exc_msg = "PFA film already in reactor"
                system_log.info(exc_msg)
                return mcs.RetObj.incomplete('Ra', mcs.V_PROBLEM, exc_msg)
            if not film_doc['clean_locations']:
                exc_msg = "No clean PFA films available"
                system_log.info(exc_msg)
                return mcs.RetObj.incomplete('Ra', mcs.V_PROBLEM, exc_msg)
            else:
                command_string = "ts" + film_doc['clean_locations'][0] + " to tr"
                updated_films['clean_locations'].pop(0)
                updated_films['location'][1] = 'occupied'
        else:
            if film_doc['location'][1] == 'empty':
                exc_msg = "No PFA film in reactor"
                system_log.info(exc_msg)
                return mcs.RetObj.incomplete('Ra', mcs.V_PROBLEM, exc_msg)
            empty_locations = list(set(PFA_FILM_LOCATIONS) - set(film_doc['clean_locations']+film_doc['dirty_locations']))
            empty_locations.sort()
            if not empty_locations:
                exc_msg = "No available PFA film storage locations"
                system_log.info(exc_msg)
                return mcs.RetObj.incomplete('Ra', mcs.V_PROBLEM, exc_msg)
            else:
                command_string = "tr to ts" + empty_locations[0]
                updated_films['dirty_locations'].append(empty_locations[0])
                updated_films['location'][1] = 'empty'
                
        with self.make_socket() as scara_sock:
            scara_sock.connect(self.sock_addr)
            sent = scara_sock.sendall(command_string.encode(ENCODING))
            if sent == 0:
                system_log.warning(f"Socket failed to transmit signal: {sent}")
                return COM_ERROR
            data = scara_sock.recv(1024).decode(ENCODING)
        if data[0] != 'S':
            system_log.info(f"Robot arm reports a failed pfa film transfer: {data}")
            return mcs.RetObj.incomplete('Ra', mcs.V_PROBLEM, "Ra film move failed")

        try:
            dbi.update_document('MC', 'consumables', film_doc, updated_films)
        except cexc.DatabaseRequestError:
            exc_msg = "Ra DB request failed - pfa film location update"
            system_log.exception(exc_msg)
            return mcs.RetObj.incomplete('Ra', mcs.V_PROBLEM, exc_msg)
        
        return mcs.RetObj.complete(data)

    def stop(self, *_, **__):
        """ Robotic Arm shutdown method

        TODO: Implement in Ra code (needs to power down robot)

        :return: Return Object
        """
        with self.make_socket() as scara_sock:
            scara_sock.connect(self.sock_addr)
            sent = scara_sock.sendall(STOP_SIGNAL)
            if sent == 0:
                system_log.warning(f"Socket failed to transmit stop signal: {sent}")
                return COM_ERROR
            data = scara_sock.recv(1024).decode(ENCODING)
        if data == 'Stopped':
            system_log.info("Robot entered safe mode")
            return mcs.RetObj.complete("Robot entered safe mode")
        else:
            system_log.warning(f"Robot failed enter safe mode: {data}")
            return mcs.RetObj.incomplete("Ra", mcs.V_FATAL, "Ra failed to enter safe mode")

    @staticmethod
    def location_translator(library_location):
        """ Utility to translate DB location names into Ra library names

        :param library_location: The location as specified by the DB
        :return: [True, str(location as specified by Ra)] or [False, "Robot-inaccessible plate location"]
        """
        if library_location[0] == 'liquid_handler':
            scara_location = 'Lh'
            if library_location[1][0] == 19:
                scara_location += str(library_location[1][1] - 5)
            elif library_location[1][0] == 24:
                scara_location += str(library_location[1][1] + 4)
            elif library_location[1][0] == 30:
                scara_location += 7
            else:
                return False, "Robot-inaccessible plate location"
        elif library_location[0] == 'fraction_collector':
            scara_location = 'fc'
            scara_location += str(library_location[1][0] - 1)   
        elif library_location[0] == 'thermal_reactor':
            scara_location = 'Th'
        elif library_location[0] == 'autosampler':
            scara_location = 'as'
        elif library_location[0] == 'lpx':
            scara_location = 'Ss'
        elif library_location[0] == 'ir_spectrometer':
            scara_location = 'Fs'
        elif library_location[0] == 'pfa_film_storage':
            scara_location = 'ts'
            scara_location += str(library_location[1])
        elif library_location[0] == 'pfa_film_reactor':
            scara_location = 'tr'
        else:
            return False, str(library_location) + " is not valid location"
        return True, scara_location
