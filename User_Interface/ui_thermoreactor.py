""" A GUI element for the thermoreactor
@author: Ben C
"""

import tkfire as tkf
from tkfire import TkFire, Dispatcher, LAYOUT, TYPE, CHILDREN
from gui_constants import *
import json
from thermoreactor_control import ThermoController
from thermoreactor_arduino import ThermoArduinoController
from omega_temp_api import OmegaTempController
from mcn_logging_manager import system_log
from yaml_box import YamlBox

# Add custom GUI element to the tkfire interpreter
thermo_dispatcher = Dispatcher(custom_elements={'YamlBox': YamlBox})


# Define GUI structure
thermo_mother = {
    'state_panel': {
        TYPE: ["LabelFrame", {'text': "Status"}],
        LAYOUT: ['pack', LB33E],
        CHILDREN: {
            'IO_frame': {
                TYPE: ['Frame', {}],
                LAYOUT: ['pack', LB33E],
                CHILDREN: {
                    'Tabula': {
                        TYPE: ["YamlBox", {'w': 20, 'h': 8, 'stretchy': True}],
                        LAYOUT: ['pack', TB33E]
                    }
                },
            },
            'Buttons_frame': {
                TYPE: ['Frame', {}],
                LAYOUT: ['pack', LB33],
                CHILDREN: {
                    'B_Load': {
                        TYPE: ['Button', {'text': "Load"}],
                        LAYOUT: ['pack', TB33],
                    },
                    'B_Save': {
                        TYPE: ['Button', {'text': "Save"}],
                        LAYOUT: ['pack', TB33],
                    },
                },
            }
        },
    },
    'control_panel': {
        TYPE: ["LabelFrame", {'text': "Control"}],
        LAYOUT: ['pack', LB33],
        CHILDREN: {
            "kwargs_frame": {
                TYPE: ["Frame", {}],
                LAYOUT: ['pack', TB33],
                CHILDREN: {
                    "kwargs_input": {
                        TYPE: ["Entry", {}],
                        LAYOUT: ['pack', TB33],
                    }
                }
            },
            "buttons_frame": {
                TYPE: ["Frame", {}],
                LAYOUT: ['pack', TB33],
                CHILDREN: {
                    "B_InitialAll": {
                        TYPE: ["Button", {'text': "Go to Initial (All)", 'width': 16}],
                        LAYOUT: ['grid', tkf.grid_arg(row=0, col=0)],
                    },
                    "B_DisconnectT": {
                        TYPE: ["Button", {'text': "Disconnect (Omega)", 'width': 16}],
                        LAYOUT: ['grid', tkf.grid_arg(row=0, col=1)],
                    },
                    "B_GCT": {
                        TYPE: ["Button", {'text': "Get Omega Temp.", 'width': 16}],
                        LAYOUT: ['grid', tkf.grid_arg(row=1, col=1)],
                    },
                    "B_SetT": {
                        TYPE: ["Button", {'text': "Set Omega Temp.", 'width': 16}],
                        LAYOUT: ['grid', tkf.grid_arg(row=2, col=1)],
                    },
                    "B_Check": {
                        TYPE: ["Button", {'text': "Check Omega", 'width': 16}],
                        LAYOUT: ['grid', tkf.grid_arg(row=3, col=1)],
                    },
                    "B_InitialOmega": {
                        TYPE: ["Button", {'text': "Go to Initial (Omega)", 'width': 16}],
                        LAYOUT: ['grid', tkf.grid_arg(row=4, col=1)],
                    },
                    "B_Offset": {
                        TYPE: ["Button", {'text': "Get Offset", 'width': 16}],
                        LAYOUT: ['grid', tkf.grid_arg(row=0, col=2)],
                    },
                    "B_Valves": {
                        TYPE: ["Button", {'text': "Get Valve Codex", 'width': 16}],
                        LAYOUT: ['grid', tkf.grid_arg(row=1, col=2)],
                    },
                    "B_Off": {
                        TYPE: ["Button", {'text': "All Off (pcb)", 'width': 16}],
                        LAYOUT: ['grid', tkf.grid_arg(row=2, col=2)],
                    },
                    "B_DisconnectPCB": {
                        TYPE: ["Button", {'text': "Disconnect (pcb)", 'width': 16}],
                        LAYOUT: ['grid', tkf.grid_arg(row=3, col=2)],
                    },
                    "B_OpenD": {
                        TYPE: ["Button", {'text': "Open Door", 'width': 16}],
                        LAYOUT: ['grid', tkf.grid_arg(row=4, col=2)],
                    },
                    "B_CloseD": {
                        TYPE: ["Button", {'text': "Close Door", 'width': 16}],
                        LAYOUT: ['grid', tkf.grid_arg(row=5, col=2)],
                    },
                    "B_Door": {
                        TYPE: ["Button", {'text': "Get Door State", 'width': 16}],
                        LAYOUT: ['grid', tkf.grid_arg(row=6, col=2)],
                    },
                    "B_OpenV": {
                        TYPE: ["Button", {'text': "Open Valve", 'width': 16}],
                        LAYOUT: ['grid', tkf.grid_arg(row=7, col=2)],
                    },
                    "B_CloseV": {
                        TYPE: ["Button", {'text': "Close Valve", 'width': 16}],
                        LAYOUT: ['grid', tkf.grid_arg(row=8, col=2)],
                    },
                    "B_Valve": {
                        TYPE: ["Button", {'text': "Get # Valves", 'width': 16}],
                        LAYOUT: ['grid', tkf.grid_arg(row=9, col=2)],
                    },
                    "B_SetFan": {
                        TYPE: ["Button", {'text': "Set Fan", 'width': 16}],
                        LAYOUT: ['grid', tkf.grid_arg(row=0, col=3)],
                    },
                    "B_Fan": {
                        TYPE: ["Button", {'text': "Get # Fans", 'width': 16}],
                        LAYOUT: ['grid', tkf.grid_arg(row=1, col=3)],
                    },
                    "B_Heat": {
                        TYPE: ["Button", {'text': "Set Heat", 'width': 16}],
                        LAYOUT: ['grid', tkf.grid_arg(row=2, col=3)],
                    },
                    "B_Press": {
                        TYPE: ["Button", {'text': "Press Piston", 'width': 16}],
                        LAYOUT: ['grid', tkf.grid_arg(row=3, col=3)],
                    },
                    "B_Release": {
                        TYPE: ["Button", {'text': "Release Piston", 'width': 16}],
                        LAYOUT: ['grid', tkf.grid_arg(row=4, col=3)],
                    },
                    "B_Middle": {
                        TYPE: ["Button", {'text': "Middle Piston", 'width': 16}],
                        LAYOUT: ['grid', tkf.grid_arg(row=5, col=3)],
                    }
                }
            }
        }
    }
}


class ThermoreactorGUI:
    """ A user interface for the thermoreactor.
    Provides access to the instrument's state and basic commands.
    """
    def __init__(self, root: tkf.tk.Tk, thermo_controller: ThermoController):
        """ Generate a UI for the thermoreactor

        :param root: A base tkinter element
        :param thermo_controller: A ThermoController object of the device being controlled
        """
        self.core = TkFire(root, {}, thermo_mother, dispatcher=thermo_dispatcher)
        self.ctrl = thermo_controller
        self.kwarg_input: tkf.tk.Entry = self.core.gui["control_panel!kwargs_frame!kwargs_input"]
        self.tabula: YamlBox = self.core.gui["state_panel!IO_frame!Tabula"]
        self.bind_state_buttons()
        self.bind_ctrl_buttons()

    def _get_ctrl_button(self, name) -> tkf.tk.Button:
        """ Fetches a control-level button by name

        :param name: the name of the button as stated in mother
        :return: self.core.gui[f"control_panel!buttons_frame!{name}"]
        """
        return self.core.gui[f"control_panel!buttons_frame!{name}"]

    def _bind_ctrl_button(self, name, command):
        """ Binds a command to a control-level button by name

        :param name: the name of the button as stated in mother
        :param command: self._get_ctrl_button(name)['command'] = command
        :return: None
        """
        self._get_ctrl_button(name)['command'] = command

    def bind_state_buttons(self):
        load_button: tkf.tk.Button = self.core.gui[f"state_panel!Buttons_frame!B_Load"]
        save_button: tkf.tk.Button = self.core.gui[f"state_panel!Buttons_frame!B_Save"]
        load_button['command'] = self.load_state
        save_button['command'] = self.save_state

    def load_state(self):
        with self.ctrl as th:
            thermo_reactor_status = th.my_status.queue[0]
        self.tabula.delete(*tkf.ST_CONTENTS)
        self.tabula.insert("1.0", thermo_reactor_status)

    def save_state(self):
        key_items = ['status', 'door', 'piston', 'tcontrol', 'valve1', 'valve2']
        input_dict = self.tabula.get(*tkf.ST_CONTENTS)
        try:
            thermo_reactor_status = {k: input_dict[k] for k in key_items}
        except (KeyError, TypeError):
            return

        with self.ctrl as th:
            th.my_status.get()
            th.my_status.put(thermo_reactor_status)

    def bind_ctrl_buttons(self):
        button_map = [
            ("B_InitialAll", self.b_initial_all),
            ("B_DisconnectT", self.b_disconnect_t),
            ("B_GCT", self.b_get_current_temperature),
            ("B_SetT", self.b_set_temperature),
            ("B_Check", self.b_check),
            ("B_InitialOmega", self.b_initial_omega),
            ("B_Offset", self.b_offset),
            ("B_Valves", self.b_valves),
            ("B_Off", self.b_all_off),
            ("B_DisconnectPCB", self.b_disconnect_pcb),
            ("B_OpenD", self.b_open_door),
            ("B_CloseD", self.b_close_door),
            ("B_Door", self.b_door),
            ("B_OpenV", self.b_open_valve),
            ("B_CloseV", self.b_close_valve),
            ("B_Valve", self.b_valve),
            ("B_SetFan", self.b_set_fan),
            ("B_Fan", self.b_fan),
            ("B_Heat", self.b_heat),
            ("B_Press", self.b_press_piston),
            ("B_Release", self.b_release_piston),
            ("B_Middle", self.b_middle_piston)
        ]
        for n, b in button_map:
            self._bind_ctrl_button(n, b)

    # ### All the Buttons for the UI ### #

    def b_initial_all(self):
        with self.ctrl as th:
            result = th.go_to_initial_state()
        system_log.info(result)

    def b_disconnect_t(self):
        with self.ctrl as th:
            t_control: OmegaTempController = th.temp_controller
            t_control.disconnect_instrument()

    def b_get_current_temperature(self):
        with self.ctrl as th:
            t_control: OmegaTempController = th.temp_controller
            result = t_control.get_current_temperature()
        system_log.info(result)

    def b_set_temperature(self):
        input_value = self.kwarg_input.get()
        if not input_value:
            system_log.info("Missing input (set_temperature expects a float)")
            return
        try:
            input_value = float(input_value)
        except (ValueError, TypeError, AttributeError):
            system_log.info("Input improperly formatted (expected setpoint as a float)")
            return

        with self.ctrl as th:
            t_control: OmegaTempController = th.temp_controller
            result = t_control.set_temperature(temp_set_point=input_value, reentrant=True)
        system_log.info(result)

    def b_check(self):
        with self.ctrl as th:
            t_control: OmegaTempController = th.temp_controller
            result = t_control.check_controller_status()
        system_log.info(result)

    def b_initial_omega(self):
        with self.ctrl as th:
            t_control: OmegaTempController = th.temp_controller
            result = t_control.initial_state()
        system_log.info(result)

    def b_offset(self):
        with self.ctrl as th:
            pcb_control: ThermoArduinoController = th.pcb_controller
            result = pcb_control.offset
        system_log.info(result)

    def b_valves(self):
        with self.ctrl as th:
            pcb_control: ThermoArduinoController = th.pcb_controller
            result = pcb_control.n_valves
        system_log.info(result)

    def b_all_off(self):
        with self.ctrl as th:
            pcb_control: ThermoArduinoController = th.pcb_controller
            result = pcb_control.all_off()
        system_log.info(result)

    def b_disconnect_pcb(self):
        with self.ctrl as th:
            pcb_control: ThermoArduinoController = th.pcb_controller
            pcb_control.disconnect_controller()

    def b_open_door(self):
        with self.ctrl as th:
            pcb_control: ThermoArduinoController = th.pcb_controller
            result = pcb_control.open_door()
        system_log.info(result)

    def b_close_door(self):
        with self.ctrl as th:
            pcb_control: ThermoArduinoController = th.pcb_controller
            result = pcb_control.close_door()
        system_log.info(result)

    def b_door(self):
        with self.ctrl as th:
            pcb_control: ThermoArduinoController = th.pcb_controller
            result = pcb_control.get_door_state()
        system_log.info(result)

    def b_open_valve(self):
        input_value = self.kwarg_input.get()
        if not input_value:
            input_value = 1
        try:
            input_value = int(input_value)
        except (ValueError, TypeError):
            system_log.info("Input improperly formatted (expected valve number as an int)")
            return

        with self.ctrl as th:
            pcb_control: ThermoArduinoController = th.pcb_controller
            result = pcb_control.open_valve(input_value)
        system_log.info(result)

    def b_close_valve(self):
        input_value = self.kwarg_input.get()
        if not input_value:
            input_value = 1
        try:
            input_value = int(input_value)
        except (ValueError, TypeError):
            system_log.info("Input improperly formatted (expected valve numbder as an int)")
            return

        with self.ctrl as th:
            pcb_control: ThermoArduinoController = th.pcb_controller
            result = pcb_control.close_valve(input_value)
        system_log.info(result)

    def b_valve(self):
        with self.ctrl as th:
            pcb_control: ThermoArduinoController = th.pcb_controller
            result = pcb_control.get_valve_internals()
        system_log.info(result)

    def b_set_fan(self):
        input_value = self.kwarg_input.get()
        if not input_value:
            system_log.info("Missing input (set_fan_speed expects an int and a float)")
            return
        try:
            input_value = json.loads(input_value)
            if isinstance(input_value, dict):
                input_value = (input_value['fan_num'], input_value['speed'])
            elif isinstance(input_value, (list, tuple)):
                pass
            else:
                system_log.info("Input improperly formatted (expected json dict or list for fan_num and speed)")
                return
        except (KeyError, ValueError, TypeError, json.JSONDecodeError):
            system_log.info("Input improperly formatted (expected json dict or list for fan_num and speed)")
            return

        with self.ctrl as th:
            pcb_control: ThermoArduinoController = th.pcb_controller
            result = pcb_control.set_fan_speed(*input_value)
        system_log.info(result)

    def b_fan(self):
        with self.ctrl as th:
            pcb_control: ThermoArduinoController = th.pcb_controller
            result = pcb_control.get_fan_internals()
        system_log.info(result)

    def b_heat(self):
        input_value = self.kwarg_input.get()
        if not input_value:
            system_log.info("Missing setpoint for heater")
            return
        try:
            input_value = float(input_value)
        except (ValueError, TypeError):
            system_log.info("Input improperly formatted (expected a number)")
            return

        with self.ctrl as th:
            pcb_control: ThermoArduinoController = th.pcb_controller
            result = pcb_control.set_heat(input_value)
        system_log.info(result)

    def b_press_piston(self):
        with self.ctrl as th:
            pcb_control: ThermoArduinoController = th.pcb_controller
            result = pcb_control.press()
        system_log.info(result)

    def b_release_piston(self):
        with self.ctrl as th:
            pcb_control: ThermoArduinoController = th.pcb_controller
            result = pcb_control.release()
        system_log.info(result)

    def b_middle_piston(self):
        with self.ctrl as th:
            pcb_control: ThermoArduinoController = th.pcb_controller
            result = pcb_control.intermediate_position()
        system_log.info(result)


if __name__ == '__main__':
    my_controller = ThermoController(None, None, None, None, False)

    thermo_gui_core = tkf.tk.Tk()
    ThermoreactorGUI(thermo_gui_core, my_controller)
    thermo_gui_core.mainloop()
