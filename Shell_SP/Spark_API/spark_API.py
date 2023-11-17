# Spark_Controller.py
# Created 12 Feb2022
# Creator: Ben Canty
# Python Version 3.7.2
# Needs to be Running As Admin in order to access the Spark's API
#
# NOTE: Elements of this code had to be redacted before public distribution, items noted with "Redacted:"

import subprocess
import sys
import time
from typing import Union

import custom_exceptions as cexc
import spectrometer_project as sp
from custom_classes import SliceDict, Narg
# from disposable_wrapper import disposable
from mcn_logging_manager import system_log
from photoreactor_API import PhotoreactorControl

try:
    import clr  # noqa # the module is 'pythonnet' but is imported as 'clr' for some reason
except ModuleNotFoundError:
    print("WARNING! This computer's pythonnet module is not good, could not import 'clr'")

# Paths and Directories
sys.path.append(r"C:\Program Files\Tecan\SparkControl\Clients")
_SPARK_AI_TARGET_DIR = r"/"
_SPARK_AI_PS1_PATH = r". 'C:\Program Files\Tecan\SparkControl\Clients\CopyAdditionalAutomationInterfaceFiles.ps1'"


# System codes
INSTRUMENT_STATE_ENUM = SliceDict()  # Redacted: Enum mapping for Spark instrument states
METHOD_STATUS_ENUM = SliceDict()  # Redacted: Enum mapping for Spark method states

# Spark
INITIAL_STATE = {
    "tray_position": "In",
    "thermostat": None,
    "encumbered": False,
    "expose": False,
}


# Redacted: Load the Spark interface

# As long as we Run As Administrator, we don't need this line
# ["powershell.exe", "Set-ExecutionPolicy -Scope CurrentUser -ExecutionPolicy Unrestricted -Force",...]
subprocess.Popen(["powershell.exe", _SPARK_AI_PS1_PATH, f"{_SPARK_AI_TARGET_DIR}"], stdout=sys.stdout).communicate()

class SparkAutomation:
    """ SparkAutomation (class: Hardware controller)

    Manager for the Spark instrument, hardware, and measurement functions
    """

    def __init__(self, event_logger=None, **kwargs):
        self.shutdown_success = None
        self.event_logger: sp.EventLogger = event_logger

        # Redacted: We begin communications with the Spark's interface...
        self.instrument = None  # ...we will map this to "self.instrument"
        self.photoreactor = PhotoreactorControl(**kwargs)

        self.state: dict = INITIAL_STATE
        """
        Dictionary for remembering the current state of the instrument:
          * tray_position - What position is the tray (In/Left/Right)
          * thermostat - Is T-control being used (None = No, # = yes)
          * encumbered - Is a plate on the tray (true/false)
          * expose - Is a plate in the photoreactor (true/false)
        """

        self.photoreactor.initialize(self.event_logger)
        self.move_plate_in()

    def _shutdown(self):
        try:
            # Redacted: We shut down the Spark's interface
            pass
        except:  # noqa (keep trucking)
            system_log.exception('Spark controller encountered an exception during shutdown (AIFactory)')
            self.shutdown_success = False

        self.shutdown_success = self.shutdown_success and self.photoreactor.shutdown()

    def shutdown(self):
        """ Shutdown protocol for the Spark.  Shuts down the Spark's interface and the forklift

        Exceptions are caught and set the self.shutdown_success flag to False
        """
        self.shutdown_success = True
        self._shutdown()

    def reset_instrument(self):
        """ Attempts to return the spark to an initial state """
        try:
            self.reserve_spark()
            self.move_plate_in()
            self.photoreactor.initialize()
            self.state.update(INITIAL_STATE)
        finally:
            self.release_spark()

    def __del__(self):
        """
        Wraps self._shutdown() with some additional logging

        Exceptions are caught and logged

        :return: None
        """
        system_log.info('Spark Controller Shutdown called')
        if not self.shutdown_success:
            self._shutdown()

    def reserve_spark(self):
        """ Reserves the Spark """
        # Redacted
        pass

    def release_spark(self):
        """ Release the Spark """
        # Redacted
        pass

    def get_instrument_status(self):
        """
        Used to populate instrument diagnostic data, returns the instrument value or repr(Exception) for each.

        :return: state, serial number, alias, type
        """
        # Redacted: Looks up the state, serial number, alias, type of the instrument to...
        return None, None, None, None  # ...return as a tuple

    def set_instrument_state(self, *,
                             tray_position: str = Narg(),
                             thermostat: Union[None, int, float] = Narg(),
                             encumbered: bool = Narg(),
                             expose: bool = Narg()):
        """
        Generalized Setter for the self.state property.

        :param tray_position: <str> Location of the tray: "In", "Right", or "Left"
        :param thermostat: <None or Number> If the thermostat is in use
        :param encumbered: <bool> If the plate is on the tray
        :param expose: <bool> If the plate is on the forklift
        :return: self.state
        """
        # Verify Arguments
        if not isinstance(tray_position, Narg):
            if not isinstance(tray_position, str):
                raise TypeError("The value of 'tray_position' must be a string")
            if tray_position not in ["In", "Right", "Left"]:
                raise ValueError(f"The value of 'tray_position' must be 'In', 'Left', or 'Right' not '{tray_position}'")
        if not isinstance(thermostat, Narg):
            if not isinstance(thermostat, (type(None), int, float)):
                raise TypeError("The value of 'thermostat' must be None or a number")
        if not isinstance(encumbered, Narg):
            if not isinstance(encumbered, bool):
                raise TypeError("The value of 'encumbered' must be a bool")
        if not isinstance(expose, Narg):
            if not isinstance(expose, bool):
                raise TypeError("The value of 'expose' must be a bool")

        # Set arguments
        if not isinstance(tray_position, Narg):
            self.state['tray_position'] = tray_position
        if not isinstance(thermostat, Narg):
            self.state['thermostat'] = thermostat
        if not isinstance(encumbered, Narg):
            self.state['encumbered'] = encumbered
        if not isinstance(expose, Narg):
            self.state['expose'] = expose

        return self.state

    def s_load_tray(self):
        """ Records a plate as being on the tray.

        :raises TrayOverload: If the tray is already occupied
        :return: True
        """
        if self.state['encumbered']:
            raise cexc.TrayOverload("Spark Tray")
        self.state['encumbered'] = True
        return True

    def s_unload_tray(self, force=False):
        """ Records a plate as being removed from the tray.

        :raises MissingPlate: If no plate was on the tray and force=False
        :param force: True-Ignore unloading an empty tray, False-Raise exception
        :return: True - success, False - If an exception would have been thrown
        """
        if not self.state['encumbered']:
            if force:
                return False
            raise cexc.MissingPlate("Spark Tray")
        self.state['encumbered'] = False
        return True

    def s_load_photoreactor(self, force=True):
        """ Records a plate as being transferred between the tray and photoreactor.

        :raises MissingPlate: If the tray had no plate
        :raises TrayOverload: If the photoreactor was already occupied
        :return: True
        """
        if not self.state['encumbered']:
            try:
                raise cexc.MissingPlate("Spark Tray")
            except cexc.MissingPlate as mp:
                if force:
                    system_log.exception("Continuing regardless")
                else:
                    raise mp
        if self.state['expose']:
            raise cexc.TrayOverload("Photoreactor")
        self.set_instrument_state(encumbered=False, expose=True)
        return True

    def s_unload_photoreactor(self, force=False):
        """ Records a plate as being transferred between the photoreactor and tray.

        :raises MissingPlate: If the photoreactor had no plate and force=False
        :raises TrayOverload: If the tray was already occupied
        :param force: True - ignore MissingPlate
        :return: True - success, False - If an exception would have been thrown
        """
        if self.state['encumbered']:
            raise cexc.TrayOverload("Spark Tray")
        self.set_instrument_state(encumbered=True)
        if not self.state['expose']:
            if force:
                return False
            raise cexc.MissingPlate("Photoreactor")
        self.set_instrument_state(expose=False)
        return True

    def move_plate_in(self, job_id=None):
        """ Moves the plate in

        :param job_id: Stamps the job's manifest
        :return: None
        """
        self.event_logger.tray_in(self, job_id)
        system_log.info("Spark moving plate in")
        # Redacted: Move spark's plate in
        self.set_instrument_state(tray_position="In")

    def move_plate_out(self, side="left", job_id=None):
        """
        Moves the plate out
        :param side: Allows the user to specify which side the plate exits from
        :param job_id: Stamps the job's manifest
        :return: True/False
        """
        system_log.info(f"Spark moving plate out on {side.lower()}")
        if side.lower() == "right":
            self.event_logger.tray_out_right(self, job_id)
            # Redacted: Move spark's plate out on the right side
            self.set_instrument_state(tray_position="Right")
            return True
        elif side.lower() == "left":
            self.event_logger.tray_out_left(self, job_id)
            # Redacted: Move spark's plate out on the left side
            self.set_instrument_state(tray_position="Left")
            return True
        else:
            return False

    def cancel(self, method_mngr):
        # Redacted: Asks Spark to cancel a measurement
        pass

    def transfer_to_photoreactor(self, job_id=None, wait_time=2.5):
        """ Moves a plate into the photoreactor

        Moves plate into spark, then moves out on right side, raises the lift, moves the tray in, then lowers the lift.
        Simultaneously with each, updates the internal status to track the plate.

        :param job_id: A key for the event logger
        :param wait_time: How long to wait between operations
        """
        self.event_logger.photoreactor_load(self, job_id)
        self.move_plate_in()  # Ensure tray is in (catches it being out-left, for example)
        time.sleep(wait_time)
        self.move_plate_out(side="right")  # tray out
        time.sleep(wait_time)
        self.photoreactor.lift_raise()  # fingers up
        time.sleep(wait_time)
        self.move_plate_in()  # tray in
        time.sleep(wait_time)
        self.photoreactor.lift_lower()  # fingers down
        self.s_load_photoreactor()

    def receive_from_photoreactor(self, job_id=None, wait_time=2.5):
        """ Moves a plate out from the photoreactor

        Moves plate into spark, raises the lift, moves tray out on right side, lowers the lift, then moves the tray in.
        Simultaneously with each, updates the internal status to track the plate.

        :param job_id: A key for the event logger
        :param wait_time: How long to wait between operations
        """
        self.event_logger.photoreactor_unload(self, job_id)
        self.s_unload_photoreactor()
        self.move_plate_in()  # Ensure tray is in (catches it being out-left, for example)
        time.sleep(wait_time)
        self.photoreactor.lift_raise()  # fingers up
        time.sleep(wait_time)
        self.move_plate_out(side="right")  # tray out
        time.sleep(wait_time)
        self.photoreactor.lift_lower()  # fingers down
        time.sleep(wait_time)
        self.move_plate_in()  # tray in

    def prepare_send(self):
        """ Prepares the Spark to send a plate.

        Moves the tray out and updates the status removing a plate from the tray.
        """
        self.move_plate_out()
        self.s_unload_tray(force=True)

    def prepare_receive(self):
        """ Prepares the Spark to receive a plate.

        Moves the tray out and updates the status adding a plate from the tray.
        """
        self.move_plate_out()
        self.s_load_tray()

    def run_plate(self, job, live_method: str) -> dict:
        """
        Runs individual scans

        :param job: dict
        :param live_method:
        :return: {run_results, error_codes, post: [instrument_state, method_return, partial_op_flag]}
        :raises RuntimeError:
        """
        system_log.info("Spark initiating a run")
        method_as_xml = live_method
        method_name = job['appello'].generate_method_name()

        i_status, *_ = self.get_instrument_status()
        if i_status in INSTRUMENT_STATE_ENUM[2:]:
            if i_status == INSTRUMENT_STATE_ENUM[2]:
                system_log.warning("Spark only partially operational")
                partial_flag = True
            else:
                partial_flag = False
        else:
            system_log.error(f"Spark not ready for operation: {i_status}")
            raise RuntimeError(i_status)

        # Redacted: Bind the method to the instrument...
        job['method_mngr'] = None # ... Save the manager ...
        is_good_method, error_codes = None, ["None", ] # ... Validate the loaded method

        if not is_good_method:
            system_log.warning("Specified method was bad")
            job['error_codes'].append(error_codes)
            return {'run_results': None,
                    'error_codes': " ".join(error_codes),
                    'post': None}  # Redacted: 'post', look up system state afterward

        system_log.info("Specified method is good")
        # Redacted: Run the method...
        run_results = None  # ...and save the result...
        state_code = None  # ...as well as the resulting system state
        self.release_spark()

        time.sleep(3)
        return {'run_results': run_results,
                'error_codes': " ".join(error_codes),
                'post': None} # Redacted: 'post', look up system state afterward
