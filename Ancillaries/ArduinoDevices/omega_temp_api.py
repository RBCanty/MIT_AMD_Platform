# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 18:16:59 2021
​
@author: Brent
@editor: Ben C (mostly for logging, locking, and comments)
"""

import minimalmodbus
import serial
import time
import datetime
from threading import Lock
from contextlib import contextmanager
from custom_classes import ParcelEvent
from mcn_logging_manager import system_log

# -----------------------------------------------------------------------------#
# ----------------- MEMORY REGISTERS BASED ON OMEGA REGISTERS -----------------#
# -------------- DO NOT CHANGE UNLESS YOU KNOW WHAT YOU ARE DOING -------------#
# -----------------------------------------------------------------------------#
memory_registers = {
    'process_value': int("4700", 16),  # Read-only
    'set_point': int("4701", 16),
    'upper_limit_alarm1': int("4702", 16),
    'lower_limit_alarm1': int("4703", 16),
    'upper_limit_alarm2': int("4704", 16),
    'lower_limit_alarm2': int("4705", 16),
    'upper_t_range': int("4706", 16),
    'lower_t_range': int("4707", 16),
    'pb_proportional_band': int("4708", 16),
    'ti_integral_time': int("4709", 16),
    'td_derivative_time': int("470A", 16),
    'heat_cool_hysteresis': int("470B", 16),
    'input_temp_sensor_type': int("4710", 16),
    'control_method': int("4711", 16),
    'heat_cool_control_cycle': int("4712", 16),
    'prop_control_offset_error_value': int("4713", 16),
    'temperature_regulation_value': int("4714", 16),
    'alarm_1_type': int("4715", 16),
    'alarm_2_type': int("4716", 16),
    'temp_unit_display_selection': int("4717", 16),
    'heat_cool_control_selection': int("4718", 16),
    'control_run_stop_setting': int("4719", 16),
    'comm_writein_selection': int("471A", 16),
    'software_version': int("471B", 16),
    'at_setting': int("4729", 16),
    'code_return_register': int("472B", 16),  # Read-only
    'ct_monitor_value': int("4733", 16)
}

mode_registers = {
    'output_status': 1,
    'input_status': 2,
    'read_register': 3,
    'force_coil': 5,
    'write_register': 6,
    'device_info': 16
}

code_returns = {
    0: 'Normal operation (No error)',
    1: 'Initial process',
    2: 'Initial status (Temperature is not stable)',
    3: 'Temperature sensor is not connected',
    4: 'Temperature sensor input error',
    5: 'Measured temperature value exceeds the temperature range',
    6: 'No Int. error',
    7: 'EEPROM Error'
}

# -----------------------------------------------------------------------------#


@contextmanager
def acquire(lock: Lock, blocking=True, timeout=-1, reentering=False):
    """
    Deploys a reentrant lock for the context manager

    :param lock: A Lock object being used by the context manager
    :param blocking: (Default: True) if the lock.acquire() call should be blocking
    :param timeout: (Default: -1) given as timeout param to lock.acquire()
    :param reentering: (Default: False) context management for if the call is reentrant
    :return: A boolean if the lock acquisition was successful (via yield for @contextmanager)
    """
    if reentering:
        success = True
    else:
        success = lock.acquire(blocking, timeout)
    if not success:
        raise RuntimeError(f"Lock failed to acquire after {timeout} seconds ({blocking})")
    try:
        yield success
    finally:
        if success and lock.locked():
            lock.release()


class OmegaTempController:
    """Class to Control Omega CN710/730/740 Controllers"""

    def __init__(self, port, safety_bounds=None):
        """
        Creates a minimalmodbus controller for an Omega CN710/730/740.

        :param port: The COM port (e.g. "COM7")
        :param safety_bounds: A tuple of (lower bound, upper bound) in degrees Celsius
        """
        self._lock = Lock()
        slave_address = 5
        self.controller = minimalmodbus.Instrument(port, slave_address, minimalmodbus.MODE_ASCII)
        self.controller.serial.baudrate = 9600
        self.controller.serial.parity = serial.PARITY_EVEN
        self.controller.serial.stopbits = 2
        self.controller.serial.bytesize = 7
        self.controller.serial.timeout = 1
        self.lock_timeout = 15  # seconds

        # load in safety bounds
        if isinstance(safety_bounds, (tuple, list)) and len(safety_bounds) == 2:
            try:
                self.safety_lower_bound_celsius = float(safety_bounds[0])
                self.safety_upper_bound_celsius = float(safety_bounds[1])
            except (TypeError, ValueError):
                system_log.exception(f"Failed to load in safety bounds for Omega specified by {safety_bounds}, "
                                     f"setting to restricted values: (5, 100)")
                self.safety_lower_bound_celsius = 5
                self.safety_upper_bound_celsius = 100
        else:
            self.safety_lower_bound_celsius = 5
            self.safety_upper_bound_celsius = 200

    def _unlock(self):
        if self._lock.locked():
            self._lock.release()

    def disconnect_instrument(self):
        """ Forces a release on the Lock, then attempts to set to initial state and close the port

        Deployed by setting the minimalmodbus's close_port_after_each_call to True then
        calling the self.initial_state() method

        :return: None
        """
        system_log.info("Disconnecting omega temperature controller")
        try:
            self._unlock()
        except RuntimeError:
            pass
        self.controller.close_port_after_each_call = True
        self.initial_state()

    def get_current_temperature(self):
        """ Requests the setpoint and current temperature

        :return: {'setpoint': float, 'temperature': float}
        """
        setpoint = self.controller.read_register(memory_registers['set_point'], 1)
        temperature = self.controller.read_register(memory_registers['process_value'], 1)
        return {'setpoint': setpoint, 'temperature': temperature}

    def set_temperature(self, temp_set_point, reentrant=False):
        """ Changes the setpoint of the controller

        :param temp_set_point: The new set point (int | float)
        :param reentrant: (Default: False) used if set_temperature() is called from a context-managed method
        :return: A dict with two keys 'old' and 'new', each maps to a dictionary of the form {'setpoint': float,
          'temperature': float} from self.get_current_temperature()
        """
        temp_set_point = float(temp_set_point)
        system_log.info(f"Setting temperature to {temp_set_point} ({reentrant})")

        with acquire(self._lock, timeout=self.lock_timeout, reentering=reentrant):
            old_values = self.get_current_temperature()

            # Make sure setpoint does not exceed safety bound
            if temp_set_point > self.safety_upper_bound_celsius:
                system_log.info(f"Set-point of {temp_set_point} exceeds maximum safe value, "
                                f"reducing set-point to upper bound ({self.safety_upper_bound_celsius})")
                temp_set_point = self.safety_upper_bound_celsius

            # A setpoint less than 5 *C is considered "off" (The Omega may enter an error state if it is given a
            #   setpoint less than or equal to zero)
            if temp_set_point <= self.safety_lower_bound_celsius:
                system_log.info(f"Turning heater off (set-point --> {self.safety_lower_bound_celsius})")
                temp_set_point = self.safety_lower_bound_celsius
            else:
                system_log.info(f"Turning heater on (set-point --> {temp_set_point})")

            # Update the setpoint
            try:
                self.controller.write_register(memory_registers['set_point'], temp_set_point, 1,
                                               functioncode=mode_registers['write_register'])
                time.sleep(1)
            except minimalmodbus.ModbusException:
                system_log.exception(f"Omega Temp Controller encountered Exception while adjusting "
                                     f"the setpoint to {temp_set_point}")

            return {'old': old_values, 'new': self.get_current_temperature()}

    def check_controller_status(self):
        """ Provides the status of the Omega controller

        References 'code_returns', a codex specified in this module.

        :return: A list(the value of the code return register, a string representation thereof)
        """
        with acquire(self._lock, timeout=self.lock_timeout):
            code_return = self.controller.read_register(memory_registers['code_return_register'])
            code_statement = code_returns.get(code_return, "Value of code return register not recognized")
            system_log.info(f"Temperature controller status: {code_return} & {code_statement}")
            return [code_return, code_statement]

    def initial_state(self):
        """ Sets the omega to a safe state (off/set-point=5)

        :return: self.get_current_temperature()
        """
        system_log.info("Resetting to initial state")
        with acquire(self._lock, timeout=self.lock_timeout):
            return self.set_temperature(5, reentrant=True)['new']

    def wait_for_temperature(self, target_temperature, threshold=5, tolerance=0.05, timeout=3600):
        """ Waits for the current temperature to match the target within a given threshold

        :param target_temperature: The temperature being targeted (note: does not have to be the setpoint)
        :param threshold: The allowable threshold for the difference between the current and target (in Celsius)
        :param tolerance: Not Implemented
        :param timeout: Maximum number of iterations to wait for the target to be reached (currently roughly seconds due
          to the time.sleep() call being 1 second)
        :return: [either "Error" or "Complete", [the time as a datetime.now(), the setpoint, the current temperature]]
        """

        _delta = timeout / 20.0
        rough_progress_bar = [int(i * _delta + 0.5) for i in range(0, 21)]

        with acquire(self._lock, timeout=self.lock_timeout):
            # A record for the progression of the temperature
            func_status_list = []
            # Counts the number of iterations, used for timeouts
            counter = 1
            # Check the temperature
            while True:
                current_temp = self.get_current_temperature()
                func_status_list.append([datetime.datetime.now(),
                                         current_temp['setpoint'],
                                         current_temp['temperature']])
                if counter in rough_progress_bar:
                    system_log.info(f"Temperature progress at "
                                    f"{current_temp['temperature']}/({current_temp['setpoint']} +/- {threshold}) "
                                    f"({int(100*round(counter/timeout, 2))}% of timeout)")
                if abs(current_temp['temperature'] - target_temperature) < threshold:
                    break
                elif counter > timeout:
                    return ['Error', func_status_list]
                else:
                    time.sleep(1)
                    counter += 1
            return ['Complete', func_status_list]

    def hold_for_time(self, target_time):
        """ Waits for a given period of time, but collects current data periodically to build a history

        Prints progress updates every ~5% of duration

        :param target_time: The duration of the hold in seconds
        :return: ["Complete", [the time as a datetime.now(), the setpoint, the current temperature]]
        """
        target_time = int(target_time)
        _delta = target_time / 20.0
        rough_progress_bar = [i * _delta for i in range(0, 21)]

        func_status_list = []
        starting_time = datetime.datetime.now()
        for i in range(0, target_time):
            current_temp = self.get_current_temperature()
            current_time = datetime.datetime.now()
            func_status_list.append([current_time,
                                     current_temp['setpoint'],
                                     current_temp['temperature']])
            elapsed_time = (current_time - starting_time).total_seconds()

            # Progress "bar"
            if any(t - elapsed_time <= 0 for t in rough_progress_bar):
                rough_progress_bar = [t for t in rough_progress_bar if t - elapsed_time > 0]  # was .pop(0) but just
                # incase multiple checkpoints are passed in one cycle of the loop, pop all non-positives
                progress = int(100*round(elapsed_time/target_time, 2))
                remaining = int((target_time - elapsed_time)/60.0 + 0.5)
                system_log.info(f"Current hold progress at {progress}% ({remaining} minutes remaining)")

            # Check criterion
            if elapsed_time >= target_time:
                system_log.info('Completed time hold!')
                break
            else:
                time.sleep(1)
        else:
            system_log.warning("FOR-Range loop exited without break, actual wait time may not have been reached "
                               "(this should be impossible as it would require a negative latency)")
        return ['Complete', func_status_list]

    def hold_for_end(self, event: ParcelEvent):
        """ Waits for a period of time specified by an Event, but collects current data periodically to build a history

        :param event: The event which governs the wait and stores the temperature profile history
        :return: ["Complete", [the time as a datetime.now(), the setpoint, the current temperature]]
        """
        func_status_list = []
        while True:
            current_temp = self.get_current_temperature()
            func_status_list.append([datetime.datetime.now(),
                                     current_temp['setpoint'],
                                     current_temp['temperature']])
            # When the event is set, save the history and exit
            if event.is_set():
                event.set_parcel(func_status_list)
                break
            else:
                time.sleep(1)
        return ['Complete', func_status_list]


if __name__ == "__main__":
    my_control = OmegaTempController("COM15")
    print(my_control.initial_state())
    my_control.disconnect_instrument()
