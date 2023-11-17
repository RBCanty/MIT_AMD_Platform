"""
@authors: Brent K, Ben C
"""

import datetime
import queue
import threading
import time

import pandas as pd

import omega_temp_api
import thermoreactor_arduino
from custom_classes import Nop, directory_check
from mcn_logging_manager import system_log


class ThermoController:
    """ Main controller for the thermoreactor """
    def __init__(self, door_com, pcb_com, valve_com, omega_com, real_mode=True):
        """ Creates a thermoreactor controller

        :param door_com: The com port for the door (if None, simulate)
        :param pcb_com: The com port for the PCB that controls the fans, piston, and heater (if None, simulate)
        :param valve_com: The com port for the valves (if None, simulate)
        :param omega_com: The com port for the omega heater
        :param real_mode: Flag to bypass and simulate all action
        """
        self._lock = threading.Lock()
        self.real_mode = real_mode
        if real_mode:
            self.temp_controller = omega_temp_api.OmegaTempController(omega_com, safety_bounds=(5, 200))
            self.pcb_controller = thermoreactor_arduino.ThermoArduinoController(door=door_com,
                                                                                pcb=pcb_com,
                                                                                valve=valve_com)
            # self.pcb_controller.open_valve(2)
        else:
            self.temp_controller = Nop()
            self.pcb_controller = Nop()
        self.my_status = queue.Queue()
        self.my_status.put({'status': 'idle', 'door': 'closed',
                            'piston': 'open', 'tcontrol': 'idle',
                            'valve1': 'closed', 'valve2': 'closed'})

    def __enter__(self):
        self._lock.acquire()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type is not None:
            system_log.warning(f"LpxController:: \n"
                               f"\t- type = {exc_type}\n"
                               f"\t- value = {exc_val}\n"
                               f"\t- traceback = {exc_tb}",
                               exc_info=True)
        self._lock.release()

    def __del__(self):
        if not self.real_mode:
            print("Simulated ThermoController:: __del__ called")
            return
        system_log.debug("ThermoController Shutdown called")
        try:
            self.temp_controller.set_temperature(-1)
            self.pcb_controller.all_off()  # door open, valves closed, heat off, fan off, piston up
            # TODO What is the best state for the gas valves to be in?
            self.pcb_controller.disconnect_controller()
        finally:
            pass

    def check_for_status(self, use_instead=None):
        """ Requests a copy of the instrument status

        :param use_instead: A queue object to use in lieu of self.my_status
        :return: [status: str, my_status: dict]
        """
        if use_instead:
            thermo_reactor_status = use_instead.queue[0]
        else:
            thermo_reactor_status = self.my_status.queue[0]
        if thermo_reactor_status['status'] == 'idle':
            return ['idle', thermo_reactor_status]
        elif thermo_reactor_status['status'] == 'busy':
            return ['busy', thermo_reactor_status]
        else:
            return ['error', thermo_reactor_status]

    def prepare_for_transfer(self):
        """ Prepares thermoreactor for a transfer

        :return: [status: str, my_status: dict]
        """
        if not self.real_mode:
            print("Simulated ThermoController:: prepare_for_transfer called")
            return ['success', {'status': 'idle', 'door': 'closed', 'piston': 'open', 'tcontrol': 'idle'}]
        thermo_reactor_status = self.my_status.queue[0]
        if thermo_reactor_status['status'] == 'busy':
            return ['busy', thermo_reactor_status]
        elif thermo_reactor_status['status'] == 'idle':
            if thermo_reactor_status['piston'] == 'closed':
                self.pcb_controller.release()
                time.sleep(10)
                thermo_reactor_status['piston'] = 'open'
            if thermo_reactor_status['piston'] == 'intermediate':
                self.pcb_controller.release()
                time.sleep(10)
                thermo_reactor_status['piston'] = 'open'
            if thermo_reactor_status['door'] == 'closed':
                self.pcb_controller.open_door()
                time.sleep(10)
                thermo_reactor_status['door'] = 'open'
            if thermo_reactor_status['valve1'] == 'closed':
                self.pcb_controller.open_valve(1)
                time.sleep(10)
                thermo_reactor_status['valve1'] = 'open'
            if thermo_reactor_status['valve2'] == 'closed':
                self.pcb_controller.open_valve(2)
                time.sleep(10)
                thermo_reactor_status['valve2'] = 'open'
            self.my_status.get()
            self.my_status.put(thermo_reactor_status)
            return ['success', thermo_reactor_status]
        else:
            return ['error', thermo_reactor_status]

    def prepare_for_film_transfer(self):
        """ Prepares thermoreactor for a transfer of the PFA film

        :return: [status: str, my_status: dict]
        """
        if not self.real_mode:
            print("Simulated ThermoController:: prepare_for_film_transfer called")
            return ['success', {'status': 'idle', 'door': 'closed', 'piston': 'open', 'tcontrol': 'idle'}]
        thermo_reactor_status = self.my_status.queue[0]
        if thermo_reactor_status['status'] == 'busy':
            return ['busy', thermo_reactor_status]
        elif thermo_reactor_status['status'] == 'idle':
            if thermo_reactor_status['piston'] == 'closed':
                self.pcb_controller.intermediate_position()
                time.sleep(10)
                thermo_reactor_status['piston'] = 'intermediate'
            if thermo_reactor_status['piston'] == 'open':
                self.pcb_controller.intermediate_position()
                time.sleep(10)
                thermo_reactor_status['piston'] = 'intermediate'
            self.my_status.get()
            self.my_status.put(thermo_reactor_status)
            return ['success', thermo_reactor_status]
        else:
            return ['error', thermo_reactor_status]

    def thermoreactor_thread_func(self, heating_operations, db_queue_name):
        """ Executes the experiment specified by heating_operations for the given queue

        heating_operations is defined by:
          - venting_time: How long to purge with N2 before heating
          - heat_temperature: The setpoint for heating
          - heat_hold_time: How long to hold at temperature before evaporating
          - evap_temperature: The setpoint for evaporating off solvents
          - evap_hold_time: How long to hold at temperature before allowing to cool
          - safe_temperature: At what temperature is it safe to end the method

        :param heating_operations: Operational details
        :param db_queue_name: Associated queue (used to define output data file name only)
        :return: None
        """
        if not self.real_mode:
            print("Simulated ThermoController:: thermoreactor_thread_func called")
            return
        start_time = datetime.datetime.now()
        thermo_reactor_status = self.my_status.queue[0]
        thermo_reactor_status['status'] = 'busy'
        self.my_status.get()
        self.my_status.put(thermo_reactor_status)
        temperature_profile = []
        if thermo_reactor_status['door'] == 'open':
            self.pcb_controller.close_door()
            thermo_reactor_status['door'] = 'closed'
            self.my_status.get()
            self.my_status.put(thermo_reactor_status)
            time.sleep(10)
            system_log.debug(str(thermo_reactor_status))
        if thermo_reactor_status['valve1'] == 'closed':
            self.pcb_controller.open_valve(1)
            time.sleep(10)
            thermo_reactor_status['valve1'] = 'open'
        if 'venting_time' in heating_operations.keys():
            if heating_operations['venting_time'] != -1:
                system_log.debug('Starting venting hold!')
                holding_return = self.temp_controller.hold_for_time(heating_operations['venting_time'])
                if holding_return[0] == 'Error':
                    # TODO: Need to go to a safe state...
                    pass
                else:
                    temperature_profile.extend(holding_return[1])
                    system_log.debug('Finished venting hold!')

        if thermo_reactor_status['piston'] == 'open':
            self.pcb_controller.press()
            time.sleep(10)
            thermo_reactor_status['piston'] = 'closed'
            self.my_status.get()
            self.my_status.put(thermo_reactor_status)
            time.sleep(10)
            system_log.debug('Piston closed!')

        if thermo_reactor_status['piston'] == 'intermediate':
            self.pcb_controller.press()
            time.sleep(10)
            thermo_reactor_status['piston'] = 'closed'
            self.my_status.get()
            self.my_status.put(thermo_reactor_status)
            time.sleep(10)
            system_log.debug('Piston closed!')

        if 'heat_hold_time' in heating_operations.keys():
            if heating_operations['heat_temperature'] != -1 and heating_operations['heat_hold_time'] != -1:
                self.temp_controller.set_temperature(heating_operations['heat_temperature'])
                waiting_return = self.temp_controller.wait_for_temperature(heating_operations['heat_temperature'])
                temperature_profile.extend(waiting_return[1])
                system_log.debug('Reached set point!')
                holding_return = self.temp_controller.hold_for_time(heating_operations['heat_hold_time'])
                temperature_profile.extend(holding_return[1])
                system_log.debug('Finished heating hold!')

        if 'evap_hold_time' in heating_operations.keys():
            if heating_operations['evap_temperature'] != -1 and heating_operations['evap_hold_time'] != -1:
                system_log.debug('Starting evaporation hold!')
                self.temp_controller.set_temperature(heating_operations['evap_temperature'])
                waiting_return = self.temp_controller.wait_for_temperature(heating_operations['evap_temperature'])
                temperature_profile.extend(waiting_return[1])
                system_log.debug('Reached set point!')
                if thermo_reactor_status['piston'] == 'closed':
                    self.pcb_controller.intermediate_position()
                    time.sleep(10)
                    thermo_reactor_status['piston'] = 'intermediate'
                    self.my_status.get()
                    self.my_status.put(thermo_reactor_status)
                    time.sleep(10)
                    system_log.debug('Piston at intermediate position!')
                holding_return = self.temp_controller.hold_for_time(heating_operations['evap_hold_time'])
                temperature_profile.extend(holding_return[1])
                system_log.debug('Finished evaporation hold!')
            else:
                if thermo_reactor_status['piston'] == 'closed':
                    self.pcb_controller.intermediate_position()
                    time.sleep(10)
                    thermo_reactor_status['piston'] = 'intermediate'
                    self.my_status.get()
                    self.my_status.put(thermo_reactor_status)
                    time.sleep(10)
                    system_log.debug('Piston at the intermediate position!')
        else:
            if thermo_reactor_status['piston'] == 'closed':
                self.pcb_controller.intermediate_position()
                time.sleep(10)
                thermo_reactor_status['piston'] = 'intermediate'
                self.my_status.get()
                self.my_status.put(thermo_reactor_status)
                time.sleep(10)
                system_log.debug('Piston at the intermediate position!')

        # Now we are preparing for the safe transfer temperature
        self.temp_controller.set_temperature(10)
        system_log.debug('Waiting for safe handling temperature!')
        waiting_return = self.temp_controller.wait_for_temperature(heating_operations['safe_temperature'])
        temperature_profile.extend(waiting_return[1])
        system_log.debug('Reached safe handling temperature')

        dataset_output_path = r'C:\Users\kfj_AMD_FTIR\Desktop\Thermal Reactor Datasets\%s' % (
            start_time.strftime('%Y%m%d'))
        directory_check(dataset_output_path)
        output_filename = dataset_output_path + r'\%s_%s_temp_profile' % (db_queue_name, start_time.strftime('%Y%m%d%H%M%S'))

        time_data = [(time_point[0] - start_time).total_seconds() for time_point in temperature_profile]
        set_point_data = [time_point[1] for time_point in temperature_profile]
        pv_temperature_data = [time_point[2] for time_point in temperature_profile]

        thermo_reactor_status['status'] = 'idle'
        self.my_status.get()
        self.my_status.put(thermo_reactor_status)

        data_frame = pd.DataFrame(list(zip(time_data, set_point_data, pv_temperature_data)),
                                  columns=['time', 'setpoint', 'temperature'])
        data_frame.to_csv(output_filename, sep=',', encoding='utf-8')

        system_log.debug('Finished!')

    def thermo_thread_spawn(self, heating_operations):
        """ Creates a thread to run the operation inside

        :param heating_operations: Specifications for the operation
        :return: ['Thermoreactor thread spawned']
        """
        thermoreactor_thread = threading.Thread(
            target=self.thermoreactor_thread_func,
            args=(heating_operations['heating_profile'], heating_operations['queue_name']),
            daemon=True
        )
        thermoreactor_thread.start()
        return ['Thermoreactor thread spawned']

    def thermo_reactor_profile_selector(self, heating_operations: dict):
        """ Container for thermo_thread_spawn

        :param heating_operations: Specifications for the operation
        :return: [key message(, further details if key message was "busy")]
        """
        # First check the queue to make sure that things are idle
        thermo_reactor_queue = self.my_status.queue[0]
        if thermo_reactor_queue['status'] == 'busy':
            return ['busy', thermo_reactor_queue]
        start_return = self.thermo_thread_spawn(heating_operations)
        if 'Thermoreactor thread spawned' not in start_return:
            return start_return

        while True:
            thermo_reactor_queue = self.my_status.queue[0]
            if thermo_reactor_queue['status'] == 'idle':
                return ['Finished operation']
            elif thermo_reactor_queue['status'] == 'error':
                return ['Thermoreactor encountered an error']
            elif thermo_reactor_queue['status'] == 'busy':
                time.sleep(5)
            else:
                time.sleep(5)

    def go_to_initial_state(self):
        """ Returns the thermoreactor to a default state

        :return: [status: str, my_status: dict]
        """
        if not self.real_mode:
            print("Simulated ThermoController:: go_to_initial_state called")
            return {'status': 'idle', 'door': 'closed', 'piston': 'open', 'tcontrol': 'idle'}
        thermo_reactor_status = self.my_status.queue[0]
        if thermo_reactor_status['status'] == 'busy':
            return ['busy', thermo_reactor_status]
        if thermo_reactor_status['door'] == 'open':
            self.pcb_controller.close_door()
            thermo_reactor_status['door'] = 'closed'
            time.sleep(10)
        if thermo_reactor_status['piston'] == 'closed':
            self.pcb_controller.release()
            thermo_reactor_status['piston'] = 'open'
            time.sleep(10)
        if thermo_reactor_status['piston'] == 'intermediate':
            self.pcb_controller.release()
            thermo_reactor_status['piston'] = 'open'
            time.sleep(10)
        if thermo_reactor_status['valve1'] == 'open':
            self.pcb_controller.close_valve(1)
            time.sleep(10)
            thermo_reactor_status['valve1'] = 'closed'
        self.my_status.get()
        self.my_status.put(thermo_reactor_status)
        if thermo_reactor_status['door'] == 'closed' and thermo_reactor_status['piston'] == 'open':
            return ["success", thermo_reactor_status]
        return ["error", thermo_reactor_status]
