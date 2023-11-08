""" Defines functionality for the Special Processes system
MAIN METHOD FOR SP

@author: Ben C
"""

from pprint import pformat
# imports the necessary modules for the FTIR, Storage Unit, Thermoreactor
import Pr_submodule as Pr
import Ss_submodule as Ss
import Th_submodule as Th
import functionals
import mcn_status as mcs
from AMD_LPX_Control import LpxController
from AMD_Spark_Control_v3 import SparkController
from mcn_logging_manager import system_log
from thermoreactor_control import ThermoController
import database_interface as dbi
import data_repository_interface as dri


def get_functions():
    """
    Called to populate the command_list of the System Object
    Think of this like a *.h method from the age of C
    :return: A dictionary of keywords : function handles
    """
    function_list = {
        'Initialize': initialize,

        'Ss.request': Ss.lpx_request,
        'Ss.stow': Ss.lpx_stow,
        'Ss.reorganize': Ss.reorganize,
        'Ss.stop': Ss.hard_stop,
        "Ss.reconnect": Ss.reconnect_lpx,
        "Ss.debug.get_locations": Ss.lpx_locations,

        'Th.prepare_to_send_receive': Th.thermo_prepare,
        'Th.run': Th.thermo_run,
        'Th.check_status': Th.thermo_check,
        'Th.go_to_initial_state': Th.thermo_initial_state,

        'run_platereader': Pr.run_platereader,
        'Pr.run_platereader': Pr.run_platereader,
        'Pr.check_wellplate_signal': Pr.check_platereader,
        'Pr.analyze_data': Pr.analyze_platereader,
        'Pr.preheat': Pr.preheat_platereader,
        'Pr.prepare_to_send_or_receive': Pr.prepare_to_receive_and_send,
        'Pr.go_to_initial_state': Pr.return_to_initial_state,
        'Pr.cancel_job': Pr.cancel_job,
        'is_spark_busy': Pr.is_spark_busy,
        'rebuild_sparkcontrol': Pr.rebuild_plate_reader,

        'debug.direct_lpx_cmd': Ss.lpx_direct,
        'debug.run_lpx_gui': Ss.open_lpx_ui,
        'debug.run_th_gui': Th.run_thermoreactor_gui,
        'debug.Th_door_open': Th.thermo_door_open,
        'debug.Th_door_close': Th.thermo_door_close,
        'debug.piston_up': Th.thermo_piston_up,
        'debug.piston_down': Th.thermo_piston_down,
        'debug.load_pfa_film': Th.pfa_film_load,
        'debug.unload_pfa_film': Th.pfa_film_unload,
        'debug.modify_spark_state': Pr.modify_spark_state,
        'debug.spark_ui': Pr.spark_ui_spawn,
        'debug.spark_lamp_on': Pr.lamp_on,
        'debug.spark_lamp_off': Pr.lamp_off,
        'debug.spark_raise': Pr.raise_forklift,
        'debug.spark_lower': Pr.lower_forklift,
        'debug.set_heater': Pr.set_photo_heater,
        'debug.move_tray': Pr.move_tray,
        'debug.read_anemometer': Pr.read_anemometer,
        'debug.calibrate_anemometer': Pr.calibrate_anemometer,
    }
    return function_list

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


def initialize(*_, **kwargs):
    _, internals, child_com = functionals.unpack_kwargs(kwargs, d_msg="")

    com_ports = internals.get([mcs.S_AUXILIS, "SP"], dict())
    system_log.info("Please confirm COM ports for startup...")

    while True:
        user_response = functionals.quick_gui(title="COM Port Confirmation",
                                              dialog=f"Please ensure that the following is correct:\n"
                                                     f"{pformat(com_ports, width=20)}",
                                              buttons={"No": lambda: None},
                                              ret_if_ok="OK")
        if user_response == "OK":
            break

        for port in com_ports:
            port_value = functionals.quick_gui(title="COM Port Update",
                                               dialog=f"Please Enter the COM port for {port}",
                                               has_entry=True)
            com_ports[port] = port_value.strip()

    with child_com as cc:
        cc['Ss'] = LpxController(com_ports["lpx_com"],
                                 r'..\Shell_SP\LPX_API\LPX.cfg',
                                 real_mode=True)
        cc['Fs'] = dict()
        cc['Th'] = ThermoController(door_com=com_ports["door_com"],
                                    pcb_com=com_ports["pcb_com"],
                                    omega_com=com_ports["omega_com"],
                                    valve_com=com_ports["valve_com"],
                                    real_mode=True)
        cc['Pr'] = dict()
        cc['Pr']['SparkControl'] = SparkController(lift_com=com_ports['lift_com'],
                                                   solar_button=com_ports['solar_button'],
                                                   solar_omega=com_ports['solar_omega'],
                                                   solar_valve=com_ports['solar_valve'],
                                                   anem_valve=com_ports['anem_valve'])

        cc['is_initialized'] = True  # TODO: use some sort of flag from LPX and ThermoController

    system_log.info(f"Testing connection to Databases: "
                    f"(Platform={dbi.test_connection()}, Results={dri.test_connection()})\n\t"
                    f"If you need to change DB settings, "
                    f"use the '__.Update_Database_Settings' command with no (kw)args")

    return mcs.RetObj.complete("Ss given controller and Fs, & Th given dictionaries")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


if __name__ == '__main__':
    from threading import Thread
    import system

    try:
        sp_system = system.System('SP')
        c4 = Thread(target=sp_system.run)
        c4.start()
        c4.join()
    except:  # noqa
        system_log.exception("[SP Post][Exception]")
