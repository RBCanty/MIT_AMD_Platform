# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 12:31:37 2020
​
@author: Brent K, Ben C
"""

import serial
import time
from threading import Lock
from custom_exceptions import ControllerError
from arduino_device import ValveController, DoorController, Th_Classes, spawn_multi_controller


class ThermoArduinoController:
    def __init__(self, *, door: str, pcb: str, valve: str, n_valves=2, inheritance: dict = None):
        """
        Creates a thermo-reactor like controller.

        :param door: COM name for the arduino which controls the door
        :param pcb: COM name for the arduino which controls the heater & fan
        :param valve: COM name for the arduino which controls the valves
        :param inheritance: A dictionary of Controller objects which should be used instead of creating new
        objects (keys: 'door', 'pcb', and 'valve' and optionally 'offset').  The 'offset' keyword is to allow the
        inheriting device to internally refer to valves starting from 1 instead of an arbitrary number (and so if valves
        are added/removed from a parent, the valve numbers do not need to be recalculated for every call).
        """
        self._lock = Lock()
        self.ports = dict()
        self.offset = 0
        self.n_valves = n_valves
        if inheritance and isinstance(inheritance, dict):
            self.offset = inheritance.pop("offset", 0)
            self.ports.update(inheritance)
        self.bind_door(door)
        self.bind_pcb(pcb)
        self.bind_valve(valve, n_valves=self.n_valves)
        time.sleep(5)

    # Construction # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    def bind_door(self, door_com):
        self.ports.setdefault('door', DoorController(com=door_com, init_ajar=False))

    def bind_pcb(self, pcb_com):
        self.ports.setdefault('pcb',
                              spawn_multi_controller(
                                  com=pcb_com,
                                  classes=Th_Classes,
                                  class_kwargs={'n_fans': 2},
                                  serial_protect_kwargs={'handled_exceptions': (IOError,),
                                                         'timeout': 5,
                                                         'post_sleep': 2}
                              )
                              )

    def bind_valve(self, valve_com, n_valves):
        self.ports.setdefault('valve', ValveController(com=valve_com, n_valves=n_valves))

    def get_inheritance_kwarg(self, door=False, pcb=False, valve=False):
        inheritance = dict()
        if door:
            inheritance.update({'door': self.ports['door']})
        if pcb:
            inheritance.update({'pcb': self.ports['pcb']})
        if valve:
            inheritance.update({'valve': self.ports['valve'], 'offset': self.n_valves})
        return inheritance

    # Deconstruction # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    def disconnect_controller(self):
        self._lock.acquire()
        try:
            for port in self.ports:
                self.ports[port].close()
        finally:
            self._lock.release()

    def all_off(self):
        self._lock.acquire()
        clean = True
        shutdown_funcs = [
            # function, args as tuple
            [self.set_heat, (0,)],
            [self.close_valve, (-1,)],
            [self.set_fan_speed, (-1,)],
            [self.release, ()],
            [self.open_door, ()],
        ]
        try:
            for _func, *_args in shutdown_funcs:
                try:
                    _func(*_args)
                    time.sleep(0.1)
                except ControllerError:
                    # The controller was never set up, ignore
                    pass
                except serial.SerialException:
                    # The controller was set up and had an error
                    clean = False
        finally:
            self._lock.release()
            return clean

    # Door # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    def open_door(self):
        if 'door' not in self.ports:
            raise ControllerError(f"Door accessed without door initialized")
        return self.ports['door'].open_door()

    def close_door(self):
        if 'door' not in self.ports:
            raise ControllerError(f"Door accessed without door initialized")
        return self.ports['door'].close_door()

    def get_door_state(self):
        if 'door' not in self.ports:
            raise ControllerError(f"Door accessed without door initialized")
        return self.ports['door'].is_open

    # Valve # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    def open_valve(self, n=1):
        if 'valve' not in self.ports:
            raise ControllerError(f"Valve {n+self.offset} accessed without valve initialized")
        return self.ports['valve'].open_valve(n+self.offset)

    def close_valve(self, n=1):
        if 'valve' not in self.ports:
            raise ControllerError(f"Valve {n+self.offset} accessed without valve initialized")
        if n < 0:
            self.ports['valve'].close_all_valves()
            return
        return self.ports['valve'].close_valve(n+self.offset)

    def get_valve_internals(self):
        if 'valve' not in self.ports:
            raise ControllerError(f"Valve accessed without valve initialized")
        return self.ports['valve'].valve_codex

    # PCB.Fan # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    def set_fan_speed(self, fan_num, speed=0):
        if 'pcb' not in self.ports:
            raise ControllerError(f"Fan accessed without pcb initialized")
        if fan_num < 0:
            self.ports['pcb'].fans_off()
            return
        return self.ports['pcb'].set_fan_speed(fan_num, speed)

    def get_fan_internals(self):
        if 'pcb' not in self.ports:
            raise ControllerError(f"Fan accessed without pcb initialized")
        return self.ports['pcb'].n_fans

    # PCB.Heat # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    def set_heat(self, heat_val):
        if 'pcb' not in self.ports:
            raise ControllerError(f"Heater accessed without pcb initialized")
        return self.ports['pcb'].set_heat(heat_val)

    # PCB.Piston # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    def press(self):
        if 'pcb' not in self.ports:
            raise ControllerError(f"Piston accessed without pcb initialized")
        return self.ports['pcb'].press()

    def release(self):
        if 'pcb' not in self.ports:
            raise ControllerError(f"Piston accessed without pcb initialized")
        return self.ports['pcb'].release()

    def intermediate_position(self):
        if 'pcb' not in self.ports:
            raise ControllerError(f"Piston accessed without pcb initialized")
        return self.ports['pcb'].intermediate_position()


if __name__ == '__main__':
    test = ThermoArduinoController(door=None, pcb=None, valve=None)
    test.set_heat(42)
    test.open_valve(n=2)
    test.set_fan_speed(1, 11)
    test.press()
    test.all_off()
