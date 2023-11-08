# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 14:00:08 2022

Author of original arduino-interfacing serial commands that defined each controller: Brooke
Ben Canty edited the methods to employ them as classes with additional safety and compatibility features for use on
the AMD platform
"""

import time
from typing import Union

from serial_intrument_interface import SerialInstrument, serial_protector
from ui_quick_dialog import QuickUI, tk
from mcn_logging_manager import system_log


def sim_override(self, caller, *_, **__):  # noqa: self as args[0] required by signature
    """
    Placeholder until the valve controller is installed

    :param self: Allows overrides to be nonstatic (not used in this case)
    :param caller: the name of the calling function (func.__name__), used to determine context
    :return: Raises RuntimeError or returns None
    """
    if caller == "open_valve":
        prompt = "Please confirm that gas is flowing"
    elif caller == "close_valve":
        prompt = "Please close valve"
    else:
        raise ValueError(f"Valve controller does not recognize caller '{caller}'")
    root = tk.Tk()
    prompt = QuickUI(root,
                     title="Gas Flow",
                     dialog=prompt,
                     buttons={"Abort": lambda *_, **__: None},
                     ret_if_ok="ok")
    ret_val = prompt.run()
    if ret_val != "ok":
        raise RuntimeError("User unable to confirm valve state")


class ExampleController(SerialInstrument):
    """
    Example for how SerialInstruments are built
    """
    def __init__(self, com, is_sim=False, additional_args="foo", **kwargs):
        super(ExampleController, self).__init__(com, is_sim, **kwargs)  # Call super to build the SerialInstrument
        self.private_variable = None  # any property should be declared here...
        self.sub_init(additional_args)  # ...and initialized here

    def sub_init(self, arg):
        # initialize properties in this method
        # if there are no properties to initialize, then this method can be omitted
        self.private_variable = arg + "-bar"

    # Create a serial protector for the class
    # (See documentation on serial_protector for more details)
    serial_protect = serial_protector((), timeout=5, post_sleep=2)

    # example of a protected method
    @serial_protect
    def say_hello(self):
        # Because the method is protected, it can just be a write() command
        # serial_protect will take care of Exceptions and waiting times as well as handling if the SerialInstrument
        #   is in simulated/real mode
        self.ser.write(b"Hello world")
        # Be careful to not call another protected method from any other protected method, or it will freeze

    # example of an unprotected method
    def say_hello_to_everyone(self, punctuation="!"):
        print(f"Hello world{punctuation}")
        self.say_hello()  # unprotected methods can call protected methods
        # (But a protected method cannot--it will freeze)


class ValveController(SerialInstrument):
    """  ValveController class (descendant of :class:`SerialInstrument`)
    Issues serial commands of the form "<&,#>" where & is a character (starting at 'L') to identify which valve and
    # is either 1 (open) or 2(close).

    Implements: open_valve(), close_valve(), and close_all_valves()

    Note: if changing the names of open_valve or close_valve, make sure to update in the sim_override method as well.
    """

    def __init__(self, n_valves, basis='L', **kwargs):
        """ Creates a multi-valve controller on a serial device.

        :param n_valves: Number of valves
        :param basis: The command basis
        :keyword com: The com port (e.g. COM7)
        :keyword is_sim: If the device is to be simulated (is_sim-->True if com is None)
        :param kwargs: Additional parameters to SerialInstrument class
        """
        super(ValveController, self).__init__(**kwargs)
        self.n_valves = n_valves
        if n_valves > 3 or n_valves < 1:
            raise ValueError("Valve controller must have between 1 and 3 (inclusive) valves connected")
        base_index = ord(basis)
        self.valve_codex = {
            **{(i + 1, 'open'): f"<{chr(base_index + i)},1>".encode() for i in range(n_valves)},
            **{(i + 1, 'close'): f"<{chr(base_index + i)},2>".encode() for i in range(n_valves)},
        }
        system_log.info(f"ValveController initialized with n_valves={n_valves} and basis={basis}")

    serial_protect = serial_protector((), timeout=5, post_sleep=2, alternate_method=sim_override)

    @serial_protect
    def open_valve(self, n=1):
        """
        Opens a valve

        :param n: The index of the valve (valves are 1-indexed)
        :return: None or Exception Object
        """
        system_log.debug(f"Opening valve {n}")
        self.ser.write(self.valve_codex[(n, 'open')])

    @serial_protect
    def close_valve(self, n=1):
        """
        Closes a valve.

        :param n: The index of the valve (valves are 1-indexed)
        :return: None or Exception Object
        """
        system_log.debug(f"Closing valve {n}")
        self.ser.write(self.valve_codex[(n, 'close')])

    def close_all_valves(self):
        for n in range(self.n_valves):
            self.close_valve(n + 1)


class AnemometerController(SerialInstrument):
    """  AnemometerController class (descendant of :class:`SerialInstrument`)
    Manages reading gas flow velocity from an anemometer.

    Implements: read_wind_speed(), calibrate_zero_point(), and check_gas_flow()
    Implements a method to check for gas flow to verify the state of the valve.

    ZERO POINT: -1.72 +/- 0
    STD.DEV: About 0.3 (at 25), 0.6 (at 41), and 1.3 (at 0.51)
    Brent-approved reading: ~25
    """

    def __init__(self, **kwargs):
        """ Creates an anemometer controller on a serial device.

        :keyword com: The com port (e.g. COM7)
        :keyword is_sim: If the device is to be simulated (is_sim-->True if com is None)
        :param kwargs: Additional parameters to SerialInstrument class
        """
        super(AnemometerController, self).__init__(**kwargs)
        self.zero = 0
        system_log.info(f"AnemometerController initialized")

    serial_protect = serial_protector((IOError, KeyError), timeout=5, post_sleep=0.9)

    @serial_protect
    def read_wind_speed(self):
        """ Protected by serial_protect

        Reports the wind speed as a float

        :return: measured wind speed
        """
        system_log.debug(f"Reading wind speed")
        self.ser.write(b'<L,2>')
        time.sleep(0.1)
        return float(self.ser.readline()
                     # Convert binary to string
                     .decode()
                     # Clean string for float()
                     .strip()
                     )

    def calibrate_zero_point(self, n=5, override=None):
        """
        Zeros the wind speed based on multiple measurements about 1 second apart

        :param n: Number of measurements
        :param override: Set's the zero point to this value manually
        :return: (mean, sample standard deviation) of calibration experiment.  If sim: returns (0,0)
        """
        # Use of override manually sets the zero point
        if override is not None:
            system_log.debug(f"Manually setting zero point to {override}")
            self.zero = float(override)
            return float(override), 0.0

        system_log.debug(f"Calibrating zero point using n={n}")

        # If simulated or n is not valid, return (0,0)
        if self.is_sim or n < 1:
            return 0.0, 0.0

        # Collect measurements and set zero point
        measurements = [self.read_wind_speed() for _ in range(n)]
        mean_gas_flow: float = sum(measurements, 0) / n
        self.zero = sum(measurements, 0) / n
        system_log.debug(f"Zero point set to {self.zero}")

        # If n permits, calculate the standard deviation
        if n < 2:
            return mean_gas_flow, None
        square_error = [(exp_gas_flow - mean_gas_flow) ** 2 for exp_gas_flow in measurements]
        std_dev: float = (sum(square_error, 0) / (n - 1)) ** 0.5

        return mean_gas_flow, std_dev

    def check_gas_flow(self, n=4, z_score=3) -> bool:
        """
        Calculates wind speed based on multiple measurements about 1 second apart

        :param n: Number of measurements
        :param z_score: Desired Z-score difference of measured from zero
        :return: bool( mean - Z*sigma > zero ).  If sim: returns True
        """
        system_log.debug(f"Checking gas flow with n={n} and a z-score threshold of {z_score}")
        if self.is_sim:
            root = tk.Tk()
            prompt = QuickUI(root,
                             title="Gas Flow",
                             dialog="Please confirm that gas is flowing",
                             buttons={"Abort": lambda *_, **__: None},
                             ret_if_ok="ok")
            ret_val = prompt.run()
            if ret_val == "ok":
                return True
            else:
                return False
        measurements = [self.read_wind_speed() for _ in range(n)]
        mean_gas_flow: float = sum(measurements, 0) / n
        noise = [exp_gas_flow - mean_gas_flow for exp_gas_flow in measurements]
        std_dev: float = (sum([n * n for n in noise], 0) / (n - 1)) ** 0.5
        return abs(mean_gas_flow - self.zero) > z_score * std_dev


class SwitchController(SerialInstrument):
    """ SwitchController class (descendant of :class:`SerialInstrument`)
    Controls the binary switch for the photoreactor

    Implements: on() and off()
    """

    def __init__(self, **kwargs):
        """ Creates a switch controller on a serial device.

        :keyword com: The com port (e.g. COM7)
        :keyword is_sim: If the device is to be simulated (is_sim-->True if com is None)
        :param kwargs: Additional parameters to SerialInstrument class
        """
        super(SwitchController, self).__init__(**kwargs)
        system_log.info(f"SwitchController initialized")

    serial_protect = serial_protector((IOError,), timeout=5, post_sleep=2)

    @serial_protect
    def on(self):
        """ Turns the lamp on

        :return: None (or IOError from serial protect)
        """
        system_log.debug("Switch turned ON")
        self.ser.write(b'<L,2>')

    @serial_protect
    def off(self):
        """ Turns the lamp off

        :return: None (or IOError from serial protect)
        """
        system_log.debug("Switch turned OFF")
        self.ser.write(b'<M,2>')


class ForkliftController(SerialInstrument):
    """ ForkliftController class (descendant of :class:`SerialInstrument`)
    Controls the Spark's photoreactor's forklift

    ImplementsL lift() and lower()
    """

    def __init__(self, is_lifted=False, **kwargs):
        """ Creates a forklift controller on a serial device.

        :param is_lifted: If the device initializes in a lifted state
        :keyword com: The com port (e.g. COM7)
        :keyword is_sim: If the device is to be simulated (is_sim-->True if com is None)
        :param kwargs: Additional parameters to SerialInstrument class
        """
        super(ForkliftController, self).__init__(**kwargs)
        self.is_lifted = is_lifted
        system_log.info(f"Initializing ForkliftController with is_lifted={is_lifted}")

    serial_protect = serial_protector((IOError,), timeout=5, post_sleep=2)

    @serial_protect
    def lift(self):
        """
        Raises the forklift if it is not already

        :return: None, A soft error message, (Or IOError from serial protect)
        """
        system_log.debug("Raising lift")
        if self.is_lifted:
            return "Forklift Already Raised"
        self.ser.write(b'<L,1>')
        self.is_lifted = True

    @serial_protect
    def lower(self):
        """
        Lowers the forklift if it is not already

        :return: None, A soft error message, (Or IOError from serial protect)
        """
        system_log.debug("Lowering lift")
        if not self.is_lifted:
            return "Forklift Already Lowered"
        self.ser.write(b'<L,2>')
        self.is_lifted = False


class DoorController(SerialInstrument):
    """ DoorController class (descendant of :class:`SerialInstrument`)
    Controls doors on custom reactors

    Implements: open_door() and close_door()
    """

    def __init__(self, basis=1, init_ajar=True, **kwargs):
        """ Creates a door controller on a serial device.

        :param basis: Defines if D1 is open or close (D2 is set opposite)
        :param init_ajar: If the door initializes open
        :keyword com: The com port (e.g. COM7)
        :keyword is_sim: If the device is to be simulated (is_sim-->True if com is None)
        :param kwargs: Additional parameters to SerialInstrument class
        """
        super(DoorController, self).__init__(**kwargs)
        self.is_open = init_ajar
        if basis == 1:
            self.door_codex = {'open': b'<D,1>', 'close': b'<D,2>'}
        elif basis == 2:
            self.door_codex = {'open': b'<D,2>', 'close': b'<D,1>'}
        else:
            raise ValueError(f"Door basis must be 1 or 2, not {basis}")
        system_log.info(f"Initializing DoorController with basis={basis}, init_ajar={init_ajar}")

    serial_protect = serial_protector((IOError,), timeout=5, post_sleep=3)

    @serial_protect
    def open_door(self):
        """
        Opens the door if not already open

        :return: None, A soft error message, (or IOError from serial protect)
        """
        system_log.debug("Opening door")
        if self.is_open:
            return "Door Already Open"
        self.ser.write(self.door_codex['open'])
        self.is_open = True

    @serial_protect
    def close_door(self):
        """
        Closes the door if not already closed

        :return: None, A soft error message, (or IOError from serial protect)
        """
        system_log.debug("Closing door")
        if not self.is_open:
            return "Door Already Closed"
        self.ser.write(self.door_codex['close'])
        self.is_open = False


class PistonController(SerialInstrument):
    """ PistonController class (descendant of :class:`SerialInstrument`)
    Controls the piston used in the thermoreactor

    Implements: press(), release(), and intermediate_position()
    """

    def __init__(self, **kwargs):
        """ Creates a piston controller on a serial device.

        :keyword com: The com port (e.g. COM7)
        :keyword is_sim: If the device is to be simulated (is_sim-->True if com is None)
        :param kwargs: Additional parameters to SerialInstrument class
        """
        super(PistonController, self).__init__(**kwargs)
        system_log.info(f"Initializing PistonController")

    serial_protect = serial_protector((IOError,), timeout=5, post_sleep=2)

    @serial_protect
    def press(self):
        """
        Moves the piston to a pressed/closed position

        :return: None (or IOError from serial protect)
        """
        system_log.debug("Pressing Piston")
        self.ser.write(b'<L,1>')

    @serial_protect
    def release(self):
        """
        Moves the piston to a raised/fully-open position

        :return: None (or IOError from serial protect)
        """
        system_log.debug("Releasing Piston")
        self.ser.write(b'<L,2>')

    @serial_protect
    def intermediate_position(self):
        """
        Moves the piston to a cracked/partially-open position

        :return: None (or IOError from serial protect)
        """
        system_log.debug("Middling Piston")
        self.ser.write(b'<L,3>')


class FanController(SerialInstrument):
    """ FanController class (descendant of :class:`SerialInstrument`)
    Controls fans on custom reactors

    Implements: set_fan_speed() and fans_off()
    """

    def __init__(self, n_fans, **kwargs):
        """ Creates a multi-fan controller on a serial device.

        :param n_fans: Number of fans
        :keyword com: The com port (e.g. COM7)
        :keyword is_sim: If the device is to be simulated (is_sim-->True if com is None)
        :param kwargs: Additional parameters to SerialInstrument class
        """
        super(FanController, self).__init__(**kwargs)
        self.n_fans = n_fans
        system_log.info(f"Initializing FanController with n_fans={n_fans}")

    serial_protect = serial_protector((IOError,), timeout=5, post_sleep=2)

    @serial_protect
    def set_fan_speed(self, fan_num, speed):
        """
        Sets a fan's speed.

        :param fan_num: Index (starts at 1) for the fan
        :param speed: The speed of the fan (units unknown)
        :return: None (or IOError from serial protect)
        """
        system_log.debuig(f"Setting fan {fan_num} to {speed}")
        message = f'<F,{fan_num},{speed}>'
        self.ser.write(message.encode())

    def fans_off(self):
        """
        Sets all fan speeds to zero

        :return:
        """
        for n in range(self.n_fans):
            self.set_fan_speed(n + 1, 0)


class ArduinoHeaterController(SerialInstrument):
    """ ArduinoHeaterController class (descendant of :class:`SerialInstrument`)
    Controls the heater on the thermoreactor

    Disambiguation: This is for the Peltiers on the thermoreactor, BinaryPeltierController is for the Peltiers on the
      photoreactor, and omega_temp_api is for the gas heaters.

    Implements: set_heat() and heat_off()
    """

    def __init__(self, **kwargs):
        """ Creates an arduino header controller on a serial device.

        :keyword com: The com port (e.g. COM7)
        :keyword is_sim: If the device is to be simulated (is_sim-->True if com is None)
        :param kwargs: Additional parameters to SerialInstrument class
        """
        super(ArduinoHeaterController, self).__init__(**kwargs)
        system_log.info(f"Initializing ArduinoHeaterController")

    serial_protect = serial_protector((IOError,), timeout=5, post_sleep=2)

    @serial_protect
    def set_heat(self, heat_val):
        """
        Sets the temperature setpoint

        :param heat_val: The temperature setpoint (units unknown)
        :return: None (or IOError from serial protect)
        """
        system_log.debug(f"Setting heater to {heat_val}")
        message = f'<H,{heat_val}>'
        self.ser.write(message.encode())

    @serial_protect
    def heat_off(self):
        """
        Sets the temperature setpoint to 0

        :return: None (or IOError from serial protect)
        """
        system_log.debug(f"Turning Heater OFF")
        self.ser.write(b'<H,0>')


class LEDController(SerialInstrument):
    """ LEDController class (descendant of :class:`SerialInstrument`)
    Controls LEDs for photoreactor and includes the ability to read a photodiode for validation

    Implements: led_off(), set_led_powe(), and read_photodiode()
    """

    def __init__(self, **kwargs):
        """ Creates an LED controller on a serial device.

        :keyword com: The com port (e.g. COM7)
        :keyword is_sim: If the device is to be simulated (is_sim-->True if com is None)
        :param kwargs: Additional parameters to SerialInstrument class
        """
        super(LEDController, self).__init__(**kwargs)
        system_log.info(f"Initializing LEDController")

    serial_protect = serial_protector((IOError,), timeout=5, post_sleep=2)

    @serial_protect
    def led_off(self):
        """ Turns LEDs off.
        :return: None (or IOError from serial protect)
        """
        system_log.debug(f"Turning LED off")
        self.ser.write(b'<B,1>')

    @serial_protect
    def set_led_power(self, power):
        """ Turns the LEDs on at a given power level.
        :param power: Power level for LEDs
        :return: None (or IOError from serial protect)
        """
        system_log.debug(f"Setting LED power to {power}")
        message = f'<F,{power}>'
        self.ser.write(message.encode())

    @serial_protect
    def read_photodiode(self, verbose=False):
        """ Measures the signal from the photodiode
        :return: A float value read from the photodiode
        :raises ValueError: If response cannot be decoded and parsed as a float
        """
        system_log.debug(f"Reading photodiode signal")
        self.ser.write(b'<H,2>')
        time.sleep(0.1)
        photo_val = float(self.ser.readline().decode().strip())
        if verbose:
            print(photo_val)
        return photo_val


class BinaryPeltierController(SerialInstrument):
    """ LEDController class (descendant of :class:`SerialInstrument`)
    Controls LEDs for photoreactor

    Disambiguation: This is for the Peltiers on the photoreactor, ArduinoHeaterController is for the Peltiers on the
      thermoreactor, and omega_temp_api is for the gas heaters.

    Implements: peltier_on() and peltier_off()
    """

    def __init__(self, **kwargs):
        """ Creates an on/off peltier controller on a serial device.

        :keyword com: The com port (e.g. COM7)
        :keyword is_sim: If the device is to be simulated (is_sim-->True if com is None)
        :param kwargs: Additional parameters to SerialInstrument class
        """
        super(BinaryPeltierController, self).__init__(**kwargs)
        system_log.info(f"Initializing BinaryPeltierController")

    serial_protect = serial_protector((IOError,), timeout=5, post_sleep=2)

    @serial_protect
    def peltier_on(self):
        """Turns peltiers on (only setpoint is 0*C).
        :return: None (or IOError from serial protect)
        """
        system_log.debug(f"Turning peltier ON")
        self.ser.write(b'<L,1>')

    @serial_protect
    def peltier_off(self):
        """ Turns peltiers off.
        :return: None (or IOError from serial protect)
        """
        system_log.debug(f"Turning peltier OFF")
        self.ser.write(b'<L,2>')


def spawn_multi_controller(com: str = None, is_sim: bool = False, serial_kwargs: dict = None,
                           classes: tuple = None, class_kwargs: dict = None,
                           serial_protect_kwargs: dict = None):
    """
    Used to merge various SerialInstruments into a single class (when one arduino/com-port controls multiple systems)

    Merging performed with type().  If multiple classes share a name (parameter or method), the instance called by
      the multi-controller will be the first in the order defined by the 'classes' parameter.

    :param com: The COM port assigned to the arduino
    :param is_sim: True/False if simulated
    :param serial_kwargs: kwargs passed to the SerialInstrument for creating the serial object
    :param classes: A tuple of tuples for each class being merged.
    :param class_kwargs: A dictionary of the kwargs passed to each class's __init__ method (can be empty)
    :param serial_protect_kwargs: A dictionary of kwarg pairs for serial_protect
    :return: a MultiController class instance with all the methods and properties of the merged classes
    """
    if serial_kwargs is None:
        serial_kwargs = {}
    serial_kwargs.update({"com": com, "is_sim": is_sim})
    if classes is None:
        classes = ()
    if class_kwargs is None:
        class_kwargs = {}
    if serial_protect_kwargs is None:
        serial_protect_kwargs = {'handled_exceptions': (), 'timeout': 5, 'post_sleep': 2}

    # create parent
    multi_controller = type("MultiController", classes, {})(**class_kwargs, **serial_kwargs)

    # create the serial_protector for the class
    multi_controller.serial_protect = serial_protector(**serial_protect_kwargs)

    return multi_controller


# Convenience definitions for the devices on the platform
Th_Classes = (PistonController, FanController, ArduinoHeaterController)
""" Classes deployed by the thermoreactor are a :class:`PistonController`, a :class:`FanController`, 
and an :class:`ArduinoHeaterController`.

The Valve and Door are on separate COM ports.

To import Th_Classes as a type, use Th_Type from this same module.
"""
Th_Type = Union[Th_Classes]

Ph_Classes = (DoorController, LEDController, BinaryPeltierController)
""" Classes deployed by the thermoreactor are a :class:`DoorController`, a :class:`LEDController`, 
and an :class:`BinaryPeltierController`.

The Valve is on a separate COM port.

To import Ph_Classes as a type, use Ph_Type from this same module.
"""
Ph_Type = Union[Ph_Classes]


if __name__ == '__main__':
    my_mc: Union[SerialInstrument, DoorController, LEDController, BinaryPeltierController, ValveController] = \
        spawn_multi_controller(
        com=None,
        classes=(DoorController, LEDController, BinaryPeltierController, ValveController),
        class_kwargs={'basis': 2, 'init_ajar': False, 'n_valves': 2}
    )

    my_mc.open_door()
    my_mc.close_door()
    my_mc.open_valve(2)
    my_mc.peltier_on()
    my_mc.set_led_power(100)
    time.sleep(5)
    my_mc.led_off()
    my_mc.peltier_off()
    my_mc.close_valve(2)
    my_mc.open_door()
    my_mc.disconnect_controller()

    print(my_mc.serial_protect, dir(my_mc.serial_protect))
