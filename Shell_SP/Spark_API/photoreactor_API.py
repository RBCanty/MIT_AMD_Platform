""" Control for the photoreactor (attached to plate reader)
@authors: Ben C, Brent K
"""
import time
from threading import Thread

from arduino_device import AnemometerController, ValveController, SwitchController, ForkliftController
from custom_classes import ParcelEvent, Nop
from mcn_logging_manager import system_log
from omega_temp_api import OmegaTempController
from spectrometer_project import EventLogger


class PhotoreactorControl:
    """ Controller class for the photoreactor """
    def __init__(self, lift_com=None, solar_button=None, solar_omega=None, solar_valve=None, anem_valve=None, **__):
        """ Creates a controller for the plate-reader--adjacent photoreactor

        :param lift_com: Com port for the list (simulates if None)
        :param solar_button: Com port for the button pusher (simulates if None)
        :param solar_omega: Com port for the omega heater (heater in gas line)
        :param solar_valve: Com port for the gas flow valve (simulates if None)
        :param anem_valve: Com port for gas flow measurement (simulates if None)
        """
        self.forklift = ForkliftController(com=lift_com)
        self.lamp = SwitchController(com=solar_button)
        if solar_omega:
            self.heater: OmegaTempController = OmegaTempController(solar_omega, safety_bounds=(5, 70))
        else:
            print("Simulated OmegaTempController:: serial opened on None")
            self.heater: OmegaTempController = Nop()
        self.valves = ValveController(com=solar_valve, n_valves=1)
        self.anemometer = AnemometerController(com=anem_valve)

        self.event_logger: EventLogger = None
        self.set_temp = 0

    def initialize(self, event_logger: EventLogger = None):
        """ Sets photoreactor into default state (lamp off, heaters off, anenometer zero point of 0)

        :param event_logger: used for collecting logs for the spectrometer projects
        :return: None
        """
        # self.forklift
        self.lamp.off()
        self.heater.initial_state()
        # self.valves.close_valve()
        time.sleep(1)
        self.anemometer.calibrate_zero_point(override=0)
        if event_logger:
            self.event_logger = event_logger

    def shutdown(self):
        """ Attempts a safe shutdown

        :return: If the shutdown was a success (bool)
        """
        shutdown_success = True

        # Forklift
        try:
            self.forklift.disconnect_controller()
        except:
            system_log.exception('Spark controller encountered an exception during shutdown (Forklift)')
            shutdown_success = False

        # Lamp
        self.set_lamp(False)
        try:
            self.lamp.disconnect_controller()
        except:
            system_log.exception('Spark controller encountered an exception during shutdown (Lamp)')
            shutdown_success = False

        # Heater
        try:
            self.heater.disconnect_instrument()
        except:
            system_log.exception('Spark controller encountered an exception during shutdown (Heater)')
            shutdown_success = False

        # Valve
        try:
            self.valves.disconnect_controller()
        except:
            system_log.exception('Spark controller encountered an exception during shutdown (Valve)')
            shutdown_success = False

        return shutdown_success

    def set_lamp(self, lamp_on=False):
        """
        Changes the state of the lamp (auto-idempotent)

        :param lamp_on: The state that the lamp should be in
        :return: None or a Serial Exception object
        """
        if lamp_on:
            return self.lamp.on()
        else:
            return self.lamp.off()

    def lift_raise(self):
        self.forklift.lift()

    def lift_lower(self):
        self.forklift.lower()

    def set_temperature(self, setpoint):
        self.set_temp = setpoint
        self.heater.set_temperature(setpoint)

    def wait_for_temperature(self, **kwargs):
        """
        References: :method:`PhotoreactorControl.set_temperature`

        :keyword target_temperature: (Optional, defaults to last call to :method:`set_temperature`)
        :keyword threshold: The acceptable window of temperatures, centered on target_temperature
        :keyword timeout: Number of seconds to allow for the target_temperature +/- threshold to be reached
        :return: ['Complete' or 'Error', [[datetime.now(), setpoint, current_temp], ...]]
        """
        if 'target_temperature' not in kwargs:
            kwargs['target_temperature'] = self.set_temp
        return self.heater.wait_for_temperature(**kwargs)

    def start_temperature_tracker(self, job_id, event: ParcelEvent):

        def track_temperature(_job_id: str, _event: ParcelEvent, logger: EventLogger):
            if logger:
                logger.stamp_manifest(
                    caller="temperature tracker",
                    message=self.heater.hold_for_end(_event),
                    job_id=_job_id,
                )
            else:
                self.heater.hold_for_end(_event)

        tracker = Thread(target=track_temperature, args=(job_id, event, self.event_logger), daemon=True)
        tracker.start()

    def is_gas_flowing(self):
        try:
            return self.anemometer.check_gas_flow()
        except (TypeError, ValueError) as bad_anemometer_exception:
            raise RuntimeError from bad_anemometer_exception

    def unlock_omega(self):
        self.heater._unlock()


if __name__ == "__main__":
    test = PhotoreactorControl()
    print("Hello to a world")
    print(test.is_gas_flowing())
    print("Goodbye to a world")
