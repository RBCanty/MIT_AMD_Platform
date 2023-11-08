""" LPX Controller
Deployment of the control system for the Storage Carousel, the LPX
Defaults
 "6": {"pitch":  785, "size": 22}, -->  785
 "8": {"pitch": 1705, "size": 10}, --> 1718
"17": {"pitch": 2447, "size":  7}, --> 2460
@author: Ben Canty
"""

import datetime
import json
import re
import sys
import time
import tkinter
import traceback
from threading import Lock

import serial

import database_interface as dbi
import ui_storage_shell
from constants import DATE_FORMAT
from custom_exceptions import DatabaseRequestError
from dummy_serial import DumSerial
from mcn_logging_manager import system_log
from ui_exception import LPXPopupDialog

ENCODING = 'Ascii'
STARTUP_SEQUENCE = [bytes(step, ENCODING) for step in [
    "ST 1900\r",        # Clear internalized values (unnecessary for cold start, but nice for restart)
    "RD DM109\r",       # LPX Serial Number
    # Set Platform Configuration ---------------------------------------------------------- Set Platform Configuration #
    "WR DM29 2\r",      # Set number of levels (2)
    "WR DM46 10\r",     # Set number of hotels per level (10)
    "ST 1604\r",        # Vertical Numbering Mode
    # Z and dZ Specifications ---------------------------------------------------------------- Z and dZ Specifications #
    "RD DM20\r",        # Handler z-Offset (default = 600) (is at 40)
    "WR DM47 19650\r",  # Set z-Offset for upper level (default = 12400)
    "WR DM22 26510\r",  # Set handler in-transfer z position (1730-original low position) (prev: 26510)
    "WR DM24 26500\r",  # Set handler out-transfer z position (1730-original low position)
    "WR DM21 177\r",    # Set handler dz pick and place movement (180)
    "WR DM26 800\r",    # Set handler dz pick and plate movement at in-transfer station (500)
    "WR DM28 800\r",    # Set handler dz pick and plate movement at out-transfer station (500)
    # Theta specifications ---------------------------------------------------------------------- Theta specifications #
    "WR DM80 60\r",     # Set left radial handler turn position (60)
    "WR DM82 3065\r",   # Set radial handler turn position at transfer station (3550) (prev: 3065)
    # Not used at the moment ------------------------------------------------------------------ Not used at the moment #
    # # "RD DM27\r",      # BCR z-Lift Read Position offset (default = 200) (is at 300)
    # # "RS 1611\r",      # Disable plate tracing
    # # "RS 1612\r",      # Disable plate waiting
    # # "RS 1613\r",      # Disable plate verification
    # Save ------------------------------------------------------------------------------------------------------ Save #
    "ST 1900\r",        # Clear internalized values (does not clear the DM values)
    "ST 1801\r",        # Initialize handling system (commits DM values to internal memory)
    "RD 1915\r"         # Signal Ready
]]
SHUTDOWN_SEQUENCE = [bytes(step, ENCODING) for step in [
    "ST 1903\r", "ST 1903\r", "ST 1903\r",
    "ST 1800\r", "ST 1800\r", "ST 1800\r",
    "ST 1900\r"
]]
CHECK_READY_SIGNAL = bytes("RD 1915\r", ENCODING)
CHECK_ERROR_SIGNAL = bytes("RD DM200\r", ENCODING)
TRANSFER_STATION = ["0", 1]
LPX_UPDATING = ['updating', -1]
LPX_INACCESSIBLE = ["0", LPX_UPDATING[0]]
RCP = 'response_code_parser'


class LpxController:
    """ Creates a controller for the LPX storage carousel """
    def __init__(self, lpx_port, config_file, real_mode=True):
        """ A controller for the LPX

        :param lpx_port: The COM port for the instrument (e.g. "COM13")
        :param config_file: A JSON dictionary containing information on allowed wellplate types and the spacing
          and number of slots on the device, as well as some error translation
        :param real_mode: If False, skips initialization
        """
        self.config_file = config_file
        self._lock = Lock()
        self.lpx_port = DumSerial()
        with open(self.config_file, "r") as fh:
            self.configuration = json.load(fh)
        self.lpx_serial_settings = {'port': lpx_port, 'baudrate': 9600, 'bytesize': 8,
                                    'stopbits': serial.STOPBITS_ONE,
                                    'parity': serial.PARITY_EVEN}
        self.last_response = ""
        #
        system_log.info("Connecting to LPX")
        if real_mode:
            self.initialize()
            system_log.info("LPX: Startup sequence complete")
        else:
            system_log.debug("DEBUG MODE (not connected to LPX)")

    def __enter__(self):
        """
        Reserves the LPX Controller for the caller

        :return: self
        """
        self._lock.acquire()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Releases the LPX from the caller

        :return: None
        """
        if exc_type is not None:
            print(f"LpxController:: \n\t- type = {exc_type}\n\t- value = {exc_val}\n\t- traceback = {exc_tb}")
            traceback.print_tb(exc_tb, file=sys.stdout)
        self._lock.release()

    def __del__(self):
        print("LPX Shutdown called")
        if self.lpx_port:
            self.disconnect(False)

    def initialize(self, reconn=False):
        """ Connects the serial device and issues a startup sequence

        :param reconn: Boolean flag for if call to initialization involves a reconnection to the instrument.  If True,
          it will attempt to close the Serial connection before creating a new one.
        :return: None
        """
        if reconn:
            print(f"reconn is {reconn}")
            if self.lpx_port:
                self.lpx_port.write(bytes("CQ\r", ENCODING))
                self.lpx_port.close()
                self.lpx_port = serial.Serial(**self.lpx_serial_settings)
            self.lpx_port.write(bytes("CR\r", ENCODING))
        self.lpx_port = serial.Serial(**self.lpx_serial_settings)
        self.lpx_port.write(bytes("CR\r", ENCODING))
        system_log.info(f"Connecting to LPX over {self.lpx_port.name}")
        time.sleep(0.250)
        while self.lpx_port.in_waiting > 0:
            self.last_response = "".join(self.lpx_port.readline().decode(ENCODING).split())
        system_log.info(f"LPX: {self.last_response}")
        system_log.info("LPX: Starting startup sequence")
        for step in STARTUP_SEQUENCE:
            print(step.decode(ENCODING)[:-1])
            self.execute(step, n_attempts=1000, wait_time=1)

    def _tune_value(self, dm, value):
        """ Used for tuning values, sets a value then updates LPX

        :param dm: A string of a number or a string of a list of numbers (e.g. "22,24")
        :param value: A string of a number to which all dm are set
        :return: None
        """
        print("starting...")
        if "," in dm:
            dm = dm.split(",")
            dm = [x.strip() for x in dm]
            sequence = [f"WR DM{x} {value}\r" for x in dm]
            sequence = sequence + ["ST 1900\r", "ST 1801\r"]
        else:
            sequence = [f"WR DM{dm} {value}\r", "ST 1900\r", "ST 1801\r"]

        for step in sequence:
            step = bytes(step, ENCODING)
            self.execute(step, n_attempts=1000, wait_time=1)
        print("...complete")
        time.sleep(1)

    def _exc_seq(self, seq):
        """ Used to execute a sequence of steps in one go

        :param seq: a list or a json-decode-able string of a list
        :return: None
        """
        if isinstance(seq, str):
            seq = seq.replace("\r", "\\r")
            print(seq)
            seq = json.loads(seq)
        for step in seq:
            print(f"Running: '{step[:-1]}'")
            step = bytes(step, ENCODING)
            self.execute(step, n_attempts=200, wait_time=1)
        print("...complete")
        time.sleep(1)

    @staticmethod
    def get_plates_from_db(_key=None):
        """ Pulls the DB for all plates with the lpx location tag

        doc, _, _ = functionals.query_location("SP", "lpx")

        :param _key: For convenience, allows the document to be directly accessed
        :return: doc[_key] or False if DB/Request error
        """
        try:
            document, err_details, resp_code = dbi.query_location("SP", 'lpx')
        except DatabaseRequestError:
            exc_msg = "Ss DB request failed - plate location lookup"
            system_log.exception(exc_msg)
            return False

        if _key is None:
            return document
        else:
            try:
                return document[_key]
            except KeyError:
                system_log.warning(f'DB did not return Error but gave a dictionary without a "{_key}" field')
                return False

    def validate_location(self, hotel, position):
        """ Used to check locations

        Transfer_Station ["0", 1] and Updating ['updating', -1] are bypassed

        :raises ValueError: If the location is not valid

        :param hotel: str/int of the hotel number
        :param position: str/int of the position number
        :return: [hotel: str, position: int]
        """
        if [hotel, position] == TRANSFER_STATION:
            return hotel, position
        if [hotel, position] == LPX_UPDATING:
            return hotel, position

        try:
            hotel = str(int(hotel))
        except ValueError:
            raise ValueError(f"Hotel '{hotel}' is not numerical")
        try:
            position = int(position)
        except ValueError:
            raise ValueError(f"Position '{position}' is not numerical")
        try:
            hotel_size = self.configuration[f'{hotel}']['size']
        except KeyError:
            raise ValueError(f"LPX asked to check position which does not exist (bad hotel: {hotel}")

        if position > hotel_size or position < 0:
            raise ValueError(f"LPX asked to check position which does not exist (bad position: {position}")

        return hotel, position

    def check_location(self, hotel, position):
        """ Helper method, looks up a given location in the database

        :param hotel: The Hotel (a string of a number)
        :param position: The position (a number)
        :return: document['full_matches']
        :raises DatabaseRequestError: If DB request fails
        """
        hotel, position = self.validate_location(hotel, position)

        try:
            document, err_details, resp_code = dbi.query_location('SP', 'lpx', [hotel, position])
        except DatabaseRequestError as dre:
            exc_msg = "Ss DB request failed - plate location lookup"
            system_log.exception(exc_msg)
            raise dre

        return document['full_matches']

    def move_plate_on_db_by_location(self, source, destination):
        """ Updates a plate location on the database

        Notably, this moves a place between locations on the LPX by source and destination (not plate name)
        Compare: `update_db_location()`

        :param source: The previous location of a plate
        :param destination: The new/current location of a plate
        :return: None
        :raises DatabaseRequestError: If any DB request fails
        """
        source = self.validate_location(*source)
        destination = self.validate_location(*destination)

        try:
            document, _, _ = dbi.query_location("SP", 'lpx', ['lpx', source])
        except DatabaseRequestError as dre:
            exc_msg = "Ss DB request failed - plate location lookup"
            system_log.exception(exc_msg)
            raise dre

        full_matches = list(document['full_matches'].items())
        if len(full_matches) != 1:
            raise ValueError(f"Database gave more than one exact match for a specified location: '['lpx', {source}]'")
        match = full_matches[0][1]['container_name']

        try:
            dbi.update_field("SP", dbi.DBG_LOCATION, match, ['lpx', source], ['lpx', destination])
        except DatabaseRequestError as dre:
            system_log.exception("Ss DB request failed - plate location update")
            raise dre

    def update_db_location(self, plate_name, location, *, source=None):
        """ Updates a plate location on the database

        Notably, this moves a plate by name and can move it to any location.
        Compare: `move_plate_on_db_by_location()`

        :param plate_name: Plate name in DB
        :param location: [location_tag, [sublocation]]
        :param source: [location_tag, [sublocation]] e.g. ['lpx', [4, 7]] or ['lpx', ['updating', -1]].  The default
          value (None) is the Updating position.
        :return:
        """
        if source is None:
            source = ['lpx', LPX_UPDATING]

        if location == ['lpx', LPX_UPDATING]:
            pass
        elif location[0] == 'off_platform':
            pass
        elif location[0] == 'lpx':
            sub_location = location[1]
            sub_location = self.validate_location(*sub_location)
            location = ['lpx', sub_location]
        else:
            pass

        try:
            dbi.update_field("SP", dbi.DBG_LOCATION, plate_name, old_value=source, new_value=location)
        except DatabaseRequestError as dre:
            system_log.exception("Ss DB request failed - plate location update")
            raise dre

    def add_plate_to_db(self, container_name, labware_type, location):
        """ Creates a new plate in the Database

        :param container_name: The name of the labware
        :param labware_type: The type of the labware
        :param location: Its location
        :return: None
        :raises DatabaseRequestError: If plate addition fails
        """
        # add validation if ever called with something other than LPX_UPDATING
        new_entry = {dbi.DBG_CONTAINER_NAME: container_name,
                     dbi.DBG_LABWARE_TYPE: labware_type,
                     dbi.DBG_LOCATION: location,
                     dbi.DBG_DATE_CREATED: datetime.datetime.now().strftime(DATE_FORMAT),
                     }
        if 'DiTi SBS Waste' in labware_type:
            tip_stack = self.configuration['tip_stack']
            tip_detail = list(range(1, tip_stack + 1))
            new_entry.update({dbi.CON_TIP_LOCATIONS: tip_detail,
                              })
            collection = dbi.COL_CONSUMABLES
        elif '96 Well Filtration Plate' in labware_type:
            new_entry.update({dbi.DBG_CONTENTS: "Empty",
                              'consumable_status': "ready",
                              'consumable_number': 1,
                              })
            collection = dbi.COL_CONSUMABLES
        else:
            new_entry.update({dbi.DBG_CONTENTS: "Empty",
                              dbi.W_BARCODE: 'None',
                              })
            collection = dbi.COL_WELLPLATES

        try:
            _, err_details, resp_code = dbi.add_document('SP', collection, new_entry)
        except DatabaseRequestError as dre:
            system_log.exception("Ss request failed - LPX could not add plate to database")
            raise dre

    def is_plate_in_db(self, container_name, labware_type):
        """ Attempts to locate a container in the database

        :param container_name: The name of the labware
        :param labware_type: The type of the labware (database collection looked up from config file)
        :return: The document if present, otherwise None.  If a DB request fails, it will return None.
        :raises ValueError: If the labware is found in multiple collections
        """
        if labware_type not in self.configuration['wellplates'] + self.configuration['consumables']:
            raise ValueError(f"Labware type '{labware_type}' not recognized by LPX")

        in_wellplates = labware_type in self.configuration['wellplates']
        in_consumables = labware_type in self.configuration['consumables']

        wellplate_document = None
        consumables_document = None

        if in_wellplates:
            try:
                wellplate_document, _, _ = dbi.query_document("SP",
                                                              dbi.COL_WELLPLATES,
                                                              dbi.DBG_CONTAINER_NAME,
                                                              container_name)
            except DatabaseRequestError:
                wellplate_document = None
        if in_consumables:
            try:
                consumables_document, _, _ = dbi.query_document("SP",
                                                                dbi.COL_CONSUMABLES,
                                                                dbi.CON_NAME,
                                                                container_name)
            except DatabaseRequestError:
                consumables_document = None

        if wellplate_document and consumables_document:
            raise ValueError("Plate found in two separate collections!")
        elif wellplate_document:
            return wellplate_document
        elif consumables_document:
            return consumables_document
        else:
            return None

    def get_first_available_slot_from_db(self, plate_type):
        """ Attempts to find the first available slot on the LPX for a given plate type

        (This is the slot to use if storing labware of this type)

        :param plate_type: The type of labware
        :return: A validated location [str, int] or None if none found
        :raises DatabaseRequestError: If DB initialization fails
        :raises ValueError: If the labware type is invalid
        """
        try:
            candidate_hotels = [str(x) for x in self.configuration[plate_type]]
        except KeyError:
            raise ValueError(f"LPX provided invalid plate type: '{plate_type}'")

        try:
            document, err_details, resp_code = dbi.query_location("SP", 'lpx')
        except DatabaseRequestError as dre:
            system_log.exception("Ss DB request failed - plate location lookup")
            raise dre

        occupied_positions = list()
        for k in document['partial_matches'].keys():
            occupied_positions.append(document['partial_matches'][k][dbi.DBG_LOCATION][1])

        for hotel in candidate_hotels:
            if hotel in LPX_INACCESSIBLE:
                continue
            candidate_positions = range(1, self.configuration[hotel]['size']+1)  # (edit_tag: ABC123)
            for position in candidate_positions:
                if [hotel, position] not in occupied_positions:
                    return self.validate_location(hotel, position)
        return None

    def get_last_slot_of_type_from_db(self, plate_type):
        """ Attempts to find the last slot on the LPX for a given plate type that is occupied by empty labware

        (This is the slot to use if request clean labware of this type)

        :param plate_type: The labware type
        :return: A location [str, int] or None if none found
        :raises DatabaseRequestError: If DB initialization fails
        :raises ValueError: If the labware type is not in the configuration file
        """
        try:
            candidate_hotels = [str(x) for x in self.configuration[plate_type]]
        except KeyError:
            raise ValueError(f"LPX provided invalid plate type: '{plate_type}'")
        candidate_hotels.reverse()

        try:
            document, err_details, resp_code = dbi.query_location("SP", "lpx")
        except DatabaseRequestError as dre:
            exc_msg = "Ss DB request failed - plate location lookup"
            system_log.exception(exc_msg)
            raise dre

        occupied_positions = list()
        for k in document['partial_matches'].keys():
            # Is the plate of the correct type
            if document['partial_matches'][k][dbi.DBG_LABWARE_TYPE] == plate_type:
                plate_contents = document['partial_matches'][k].get(dbi.DBG_CONTENTS, None)
                # Is the plate empty (contents ="Empty", ={}, or =None; or 'contents' isn't even a key--like for tips)
                if plate_contents == 'Empty' or not plate_contents:
                    occupied_positions.append(document['partial_matches'][k][dbi.DBG_LOCATION][1])
        # print("DEBUG(get_last_slot_of_type_from_db):", occupied_positions)

        for hotel in candidate_hotels:
            candidate_positions = list(range(1, self.configuration[hotel]['size']+1))
            # ^^^ Think this fixes the access problem (edit_tag: ABC123)
            candidate_positions.reverse()
            for position in candidate_positions:
                if [hotel, position] in occupied_positions:
                    return [hotel, position]
        return None

    def disconnect(self, loud: bool = True):
        """ Sends a disconnect signal to the LPX and closes serial connection

        :param loud: Boolean if messages should be printed to the logger (True) or to only the console (False)
        :return: None
        """
        if loud:
            emitter = system_log.info
        else:
            emitter = print
        emitter("Disconnecting LPX...")
        self.lpx_port.write(bytes("CQ\r", ENCODING))
        self.lpx_port.close()
        emitter("...LPX disconnected")

    def execute(self, cmd, *, force_run=False, wait_time=0.250, poll_time=0.150, n_attempts=200):
        """ Attempts to execute a command on the LPX.  Checks the Ready Flag before submitting a command (unless forced)
        and will check a specified number of times (n_attempts).  Separate commands are spaced by a wait time (wait).
        This function also manages receiving the response codes from the LPX after a command is submitted.

        :param cmd: a string or bytes representing the command sent to port 'port'
        :param force_run: submits a command without checking if the port is ready
        :param wait_time: specifies the amount of time to wait before sending the command (>200 ms)
        :param poll_time: specifies the amount of time to wait between ready requests (~100 ms)
        :param n_attempts: specifies how many ready requests can be sent before timeout
        :return: (-1: state unknown), (0: command successfully sent),
        (1: command timed out), (2: error)
        """
        if not isinstance(cmd, bytes):
            cmd = cmd.encode(ENCODING)
        # Wait >200 ms from last command before requesting the ready bit
        time.sleep(wait_time)
        # Issue command
        if force_run:
            # If forced, send and return
            system_log.info("LPX: Writing command: " + cmd.decode(ENCODING).strip() + " (forced)")
            self.lpx_port.write(cmd)
            time.sleep(wait_time)
            if self.lpx_port.in_waiting > 0:
                print("Forced-run response:", "".join(self.lpx_port.readline().decode(ENCODING).split()))
            return -1
        else:
            # If not forced
            # Wait until system says it's ready unless n_attempts <= 0 (reserved for error handling)
            print(f"DEBUG:: Polling ready for cmd '{cmd.decode(ENCODING)}'")
            if self.poll_ready(n_attempts, poll_time):
                system_log.info("LPX: Writing command: " + cmd.decode(ENCODING).strip())
                self.lpx_port.write(cmd)
            else:
                return 1
            # Get response
            time.sleep(wait_time)
            while self.lpx_port.in_waiting > 0:
                self.last_response = "".join(self.lpx_port.readline().decode(ENCODING).split())
                system_log.info(f"LPX response: {self.last_response}")
                if cmd == "RD 1915\r".encode(ENCODING):
                    # If this was just a "RD 1915" ping, by virtue of self.poll_ready -> True, we are good
                    return 0
                elif self.last_response in ['0', 'OK']:
                    return 0
                elif (self.last_response[0] == 'E' or self.last_response.isnumeric()) \
                        and (cmd.decode(ENCODING)[:2] != "RD"):
                    rcp_key = str(self.last_response)
                    delim = "\n\t"
                    system_log.error(f"LPX ERROR: '{self.configuration[RCP].get(rcp_key, rcp_key)}'\n"
                                     f"Details:\n\t{delim.join([str(deets) for deets in self.get_error_details()])}")
                    return 2
                else:
                    return -1

    def poll_ready(self, n_attempts, poll_time):
        """ Periodically check the ready signal for a positive result

        :param n_attempts: How many times to check
        :param poll_time: How long to wait between checks (in seconds)
        :return: True if ready, False otherwise
        """
        if n_attempts == 0:
            return True
        for _ in range(n_attempts):
            self.lpx_port.write(CHECK_READY_SIGNAL)
            time.sleep(poll_time)
            if self.lpx_port.in_waiting > 0:
                poll_response = "".join(self.lpx_port.readline().decode(ENCODING).split())
                # If the LPX says it's ready, send it the command
                if poll_response == '1':
                    return True
        else:
            # 'else' only runs if the 'break' isn't invoked, which means a time-out
            system_log.warning(f"LPX: LPX not responding; last known response:\n"
                               f"\t{self.last_response}")
            return False

    def get_error_details(self):
        """ Sends the check error signal then attempts to report the error info in human-readable form

        :return: (Uses YIELD) The literal error response, a human-readable form (or the literal error response again
          if the lookup fails)
        """
        time.sleep(0.250)
        self.lpx_port.write(CHECK_ERROR_SIGNAL)
        time.sleep(0.150)
        counter = 0
        while self.lpx_port.in_waiting > 0:
            if counter > 5:
                yield '-1', 'LPX overflow'
                return
            err_response = "".join(self.lpx_port.readline().decode(ENCODING).split())
            err_response = str(err_response)
            counter += 1
            yield err_response, self.configuration[RCP].get(err_response, err_response)

    def _generate_transfer_sequence(self, hotel_number: int, position_number: int, direction: int):
        """ Generates a sequence of commands for the LPX to execute which accomplish the goal of transferring
        a plate from one location to another (between the transfer station and the storage unit).

        :param hotel_number: Hotel number 1-20
        :param position_number: Plate number 1-7, 1-10, or 1-22 (depending on hotel)
        :param direction: -1 for accessing a plate (plate leaves LPX) and +1 for stowing a plate (plate enters LPX)
        :return: The transfer sequence as a list(str)
        """
        hotel_size = self.configuration[str(hotel_number)]["size"]
        hotel_pitch = self.configuration[str(hotel_number)]["pitch"]
        slot = (hotel_number - 1) * hotel_size + position_number
        if direction > 0:
            action = "WR DM10 " + str(slot) + chr(13)
        elif direction < 0:
            action = "WR DM15 " + str(slot) + chr(13)
        else:
            # direction == 0
            system_log.error(f"Transfer sequence does not recognize direction {direction}")
            return []
        transfer_sequence = [
            "WR DM25 " + str(hotel_size) + chr(13),  # Set number of levels
            "WR DM23 " + str(hotel_pitch) + chr(13),  # Set pitch
            action,  # Import/export plate
            "RD 1915\r"  # Confirm success
        ]
        transfer_sequence = [bytes(item, ENCODING) for item in transfer_sequence]
        return transfer_sequence

    def generate_transfer_sequence(self, hotel_number, position_number, direction):
        """ Generates a sequence of commands for the LPX to execute which accomplish the goal of transferring
        a plate from one location to another (between the transfer station and the storage unit).

        :param hotel_number: Hotel number 1-20
        :param position_number: Plate number 1-7, 1-10, or 1-22 (depending on hotel)
        :param direction: -1 for accessing a plate (plate leaves LPX), +1 for stowing a plate (plate enters LPX), and 0
        to reposition the selected plate
        :return: The transfer sequence as a list(str)
        """
        hotel_number, position_number = self.validate_location(hotel_number, position_number)
        if hotel_number in LPX_INACCESSIBLE:
            system_log.warning(f"Transfer sequences not set up for {hotel_number}, must be a hotel #1-20")
            raise LPXPopupDialog
        hotel_number = int(hotel_number)
        position_number = int(position_number)
        direction = int(direction)

        if direction == 0:
            plate_type = self.configuration[str(hotel_number)]['occupancy'][position_number - 1]
            if plate_type == 0:
                system_log.warning("LPX asked to relocate plate which doesn't exist")
                return []
            new_hotel_number, new_position_number = self.get_first_available_slot_from_db(plate_type)
            # There is a better way to do this that is faster as it doesn't take up the transfer station but I don't
            # remember all the LPX commands
            part1 = self._generate_transfer_sequence(hotel_number, position_number, -1)
            part2 = self._generate_transfer_sequence(new_hotel_number, new_position_number, 1)
            return part1 + part2
        else:
            return self._generate_transfer_sequence(hotel_number, position_number, direction)

    def is_transfer_station_occupied(self):
        """ Used to check if the transfer station is clear or occupied

        :return: (-1 if error -> True), (0 if clear -> False), (1 if occupied -> True)
        """
        time.sleep(0.250)
        self.lpx_port.write(bytes("RD 1813\r", ENCODING))
        time.sleep(0.150)
        while self.lpx_port.in_waiting > 0:
            self.last_response = "".join(self.lpx_port.readline().decode(ENCODING).split())
            system_log.info(f"LPX: {self.last_response}")
            if re.search('[0-9]{5}', self.last_response) or (self.last_response[0] == 'E'):
                system_log.error(f"LPX Error:\n"
                                 f"\t{self.configuration[RCP][str(self.last_response)]}")
                return -1
            try:
                return int(self.last_response)  # It's 0 for no plate and 1 for plate
            except ValueError:
                system_log.exception(f"LPX Exception: When checking transfer station received an unexpected response: "
                                     f"'{self.last_response}' which caused an exception.")
                return -1

    def load_onto_carousel(self, ):
        """ Opens the storage carousel GUI

        :return: True if a change has occured, False otherwise
        """
        root = tkinter.Tk()
        ui = ui_storage_shell.StorageShellUserInterface(root, self)
        code = ui.run()
        return code

    def carousel_update_from_user(self, manifest, selected_hotel):
        """ Takes a manifest from the GUI and executes the required changes

        :param manifest: [[name, type, [hotel, slot]], ...]
        :param selected_hotel: The hotel being modified
        :return: True/False - if there was a change
        """

        there_was_a_change = False
        big_moves = list()

        # Run through DB and check if name exists
        #     If YES, set its position to ["lpx", ["updating", -1]]
        #     If NO, create a new plate in position ["lpx", ["updating", -1]]
        for plate in manifest:
            if plate[0] == '' or plate[1] == '':
                continue
            extant = self.is_plate_in_db(plate[0], plate[1])
            if not extant:
                self.add_plate_to_db(plate[0], plate[1], ['lpx', LPX_UPDATING])
            else:
                if 'lpx' not in extant['location']:
                    big_moves.append([extant['container_name'], extant['location']])

        relocated_in_place = [plate[0] for plate in manifest]

        # Garbage collect
        partial_matches = self.get_plates_from_db('partial_matches')
        for k in partial_matches.keys():
            hotel, slot = partial_matches[k]['location'][1]
            if hotel == selected_hotel:
                if partial_matches[k]['container_name'] in relocated_in_place:
                    self.update_db_location(partial_matches[k]['container_name'],
                                            ['lpx', LPX_UPDATING],
                                            source=['lpx', [hotel, slot]])
                else:
                    self.update_db_location(partial_matches[k]['container_name'],
                                            ['off_platform',
                                             ['removed by user during LPX reloading and not replaced', ]],
                                            source=['lpx', [hotel, slot]])
                there_was_a_change = True

        # Handle Big moves
        # (where the plate isn't moving within/on/off the LPX but is jumping from another location entirely)
        for name, loc in big_moves:
            self.update_db_location(name, ['lpx', LPX_UPDATING], source=loc)

        # Run through manifest and update the plate locations from ["lpx", ["updating", -1]] to ["lpx", [hotel, slot]]
        for plate in manifest:
            self.update_db_location(plate[0], ['lpx', plate[2]], source=['lpx', LPX_UPDATING])

        return there_was_a_change

    def emergency_stop(self):
        manifest = list()
        for step in SHUTDOWN_SEQUENCE:
            try:
                a = self.execute(step, force_run=True)
            except Exception as e:
                a = repr(e)
            manifest.append([step, a])
        return manifest

    @staticmethod
    def generate_barcode_sequence():
        sequence = ["ST 1714", "ST 1712", "ST 1713", "RD 1915"]
        return [bytes(item + chr(13), ENCODING) for item in sequence]


if __name__ == '__main__':
    from pprint import pprint

    cfg_path = r".\LPX.cfg"
    my_controller = LpxController("COM10", cfg_path, real_mode=True)

    extra = ['__del__', '__dict__', '__enter__', '__exit__', '__module__', '__weakref__',
             'state', 'config_file', "_lock", "lpx_port", "configuration", "lpx_serial_settings",
             "last_response", "last_error_code", "_generate_transfer_sequence", ]

    funcs = [x for x in dir(my_controller) if x not in dir(object) + extra]
    pprint(funcs)
    macros = {"$E": 'execute -cmd="RD 1813\r" -force_run=true',
              "$0": 'execute -cmd="ST 1903\r" -force_run=true',
              "$1": '_exc_seq -seq=["WR DM25 7\r", "WR DM23 2460\r", "WR DM15 126\r", "RD 1915\r"]',  # H:18, P:7, D-1
              "$2": '_exc_seq -seq=["WR DM25 7\r", "WR DM23 2460\r", "WR DM10 126\r", "RD 1915\r"]',  # H:18, P:7, D+1
              }

    # TODO: Just make a response parse/interpreter method

    while True:

        command = input("Enter Command (-exit to exit, -help for command list)\n>> ")

        if command == "-exit":
            break
        elif command == "-help":
            pprint(funcs)
            time.sleep(1)
            continue
        elif command in macros:
            command = macros[command].replace(":", "")

        kwargs = dict()

        if " -" in command:
            try:
                command, *args = command.split(" -")
            except ValueError:
                time.sleep(1)
                continue

            for arg in args:
                try:
                    try:
                        main_k, main_v = arg.split("=")
                    except ValueError:
                        main_k, main_v = arg.split(":")
                except:  # noqa
                    break

                kwargs[main_k] = main_v
                if main_v != main_v.strip("'"):
                    kwargs[main_k] = main_v.strip("'")
                elif main_v != main_v.strip("\""):
                    kwargs[main_k] = main_v.strip("\"")

                kwargs[main_k] = kwargs[main_k].replace("\\r", "\r")

        try:
            func = getattr(my_controller, command)

            if kwargs:
                pass
            else:
                print(getattr(func, "__doc__"))
                try:
                    kwargs = json.loads(input("Enter Kwargs as a json dictionary\n>> "))
                    print(kwargs)
                except json.decoder.JSONDecodeError as jde:
                    traceback.print_tb(jde.__traceback__)
                    continue
            ret_val = func(**kwargs)
            print(ret_val)
        except Exception as main_e:
            print(f"Exception\n{repr(main_e)}\n")
            traceback.print_tb(main_e.__traceback__)

        time.sleep(1)

    del my_controller

# WR DM123 63450 ; Activate service mode
# ST 1701 ; shovel out
# RS 1701 ; shovel in
# WR DM123 0 ; Leave service mode
#
# RD DM22; reads transfer height
# WR DM123 63450; puts it in the service mode
# WR DM6 [value of DM22]
# ST 1702; turns out
# ST 1701; shovel out
# RS 1701; shovel in
# RS 1702; turn back in
# WR DM123 0; exit service mode
#
# 'WR DM25 10' = Hotel 12
# 'WR DM23 1705' = Pitch
# 'WR DM15 113' = (size)*(12-1) + Position
# 'RD 1915'
#

# CR - Comm Req.  -> CC \cr \lf
# CQ - Comm Quit  -> CF \cr \lf
#
# ST - Set        -> "OK"
# RS - Reset      -> "OK"
# RD - Read       -> "&&&&&"
# WR - Write      -> "OK"
#                 -> or "E#"
# - - - -
# RD 1915         -> 1 (if ready) else 0

#
# ST and RS always return "OK"
#
