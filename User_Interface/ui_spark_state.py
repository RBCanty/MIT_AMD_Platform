""" A GUI for controlling the Spark and its associated equipment
@author: Ben C

Note on 'safety': While the use of a singleton pattern for the class would accomplish the goal of ensuring that there is
only one controller-GUI open at a time, for the sake of a nice UI where the button to spawn it also grays out, the
safety parameter is used instead.
"""

import tkinter as tk
from typing import Union, Dict

import message as msgm
import operations as oprtn
from gui_constants import *

GP3EW = {'sticky': "ew", **PADDING_3}


class SparkStateUI:
    """ A GUI for controlling the Spark and some state variables """
    def __init__(self, root: Union[tk.Tk, tk.Toplevel], spark, outbox=None, safety=None):
        """ Allows getting and setting instrument states as well as giving basic commands to the peripheral devices

        :param root: a tkinter root
        :param spark: either: (the controller for the Spark, needed for reading state information) or (None)
        :param outbox: If spark is None, the outbox is needed to tell the System to send the commands in the UI's stead
        :param safety: For use with the main GUI, to prevent spawning multiple instances
        """
        self.outbox = outbox
        self.safety = safety
        if self.safety:
            self.safety['state'] = tk.DISABLED

        self.root = root
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)
        master_frame = tk.Frame(root)
        master_frame.pack()
        self.ret_val: Dict[str, Union[str, float, bool, None]] = {}
        self.spark = spark if type(spark).__name__ == 'SparkAutomation' else None

        if spark is not None:
            state_dict = self.spark.state
            self.tray_position = tk.StringVar(value=state_dict['tray_position'])
            self.thermostat = str(state_dict['thermostat'])
            self.encumbered = tk.BooleanVar(value=state_dict['encumbered'])
            self.expose = tk.BooleanVar(value=state_dict['expose'])

            state_frame = tk.LabelFrame(master_frame, text="Set state")
            state_frame.pack(**TB33)

            tk.Label(state_frame, text="Tray Position").grid(row=0, column=0, columnspan=2)
            tk.Label(state_frame, text="Thermostat").grid(row=1, column=0, columnspan=2)
            tk.Label(state_frame, text="Is a plate on the tray").grid(row=2, column=0, columnspan=2)
            tk.Label(state_frame, text="Is a plate in the photoreactor").grid(row=3, column=0, columnspan=2)

            tk.Radiobutton(state_frame, text="In", variable=self.tray_position, value="In").grid(row=0, column=2)
            tk.Radiobutton(state_frame, text="Out Left", variable=self.tray_position, value="Left").grid(row=0, column=3)  # noqa
            tk.Radiobutton(state_frame, text="Out Right", variable=self.tray_position, value="Right").grid(row=0, column=4)  # noqa

            self.therm_entry = tk.Entry(state_frame)
            self.therm_entry.grid(row=1, column=2)
            self.therm_entry.insert(0, self.thermostat)

            tk.Radiobutton(state_frame, text="Yes", variable=self.encumbered, value=True).grid(row=2, column=2)
            tk.Radiobutton(state_frame, text="No", variable=self.encumbered, value=False).grid(row=2, column=3)

            tk.Radiobutton(state_frame, text="Yes", variable=self.expose, value=True).grid(row=3, column=2)
            tk.Radiobutton(state_frame, text="No", variable=self.expose, value=False).grid(row=3, column=3)

            tk.Button(state_frame, text="Submit", command=self.submit).grid(row=4, column=0)
            tk.Button(state_frame, text="Cancel", command=self.cancel).grid(row=4, column=1)

        if self.outbox is not None:
            ctrl_frame = tk.LabelFrame(master_frame, text="Send Commands")
            ctrl_frame.pack(**TB33)

            tk.Button(ctrl_frame, text="Plate In", command=lambda: self.cmd('debug.move_tray', mode="in")) \
                .grid(row=5, column=0, **GP3EW)
            tk.Button(ctrl_frame, text="Plate Out Left", command=lambda: self.cmd('debug.move_tray', mode="left")) \
                .grid(row=5, column=1, **GP3EW)
            tk.Button(ctrl_frame, text="Plate Out Right", command=lambda: self.cmd('debug.move_tray', mode="right")) \
                .grid(row=5, column=2, **GP3EW)

            tk.Button(ctrl_frame, text="Lamp On", command=lambda: self.cmd('debug.spark_lamp_on'))\
                .grid(row=6, column=0, **GP3EW)
            tk.Button(ctrl_frame, text="Lamp Off", command=lambda: self.cmd('debug.spark_lamp_off'))\
                .grid(row=6, column=1, **GP3EW)
            tk.Button(ctrl_frame, text="Heater Off", command=lambda: self.cmd('debug.set_heater')) \
                .grid(row=6, column=2, **GP3EW)
            tk.Button(ctrl_frame, text="Read Gas", command=lambda: self.cmd('debug.read_anemometer')) \
                .grid(row=6, column=3, **GP3EW)
            tk.Button(ctrl_frame, text="Calibrate Gas", command=lambda: self.cmd('debug.calibrate_anemometer')) \
                .grid(row=6, column=4, **GP3EW)
            tk.Button(ctrl_frame, text="Raise Lift", command=lambda: self.cmd('debug.spark_raise'))\
                .grid(row=6, column=5, **GP3EW)
            tk.Button(ctrl_frame, text="Lower Lift", command=lambda: self.cmd('debug.spark_lower'))\
                .grid(row=6, column=6, **GP3EW)

    def on_closing(self, *_):
        """ Restores the GUI ability to spawn another instance and closes the UI

        :return:
        """
        self.root.destroy()
        if self.safety:
            self.safety['state'] = tk.NORMAL

    def submit(self, *_):
        """ Loads int state-inputs and commits to the return value before closing

        :return:
        """
        tray_position = self.tray_position.get()
        thermostat = self.therm_entry.get()
        encumbered = bool(self.encumbered.get())
        expose = bool(self.expose.get())

        thermostat = None if thermostat == 'None' else float(thermostat)

        self.ret_val = {
            'tray_position': tray_position,
            'thermostat': thermostat,
            'encumbered': encumbered,
            'expose': expose,
        }
        self.on_closing()

    def cancel(self, *_):
        """ Clears state-inputs before closing

        :return:
        """
        self.ret_val = {}
        self.on_closing()

    def run(self):
        """ Executes popup

        If controller is present, will commit state-level changes

        :return: A dictionary of the state-inputs
        """
        self.root.mainloop()
        if self.spark and self.ret_val:
            self.spark.set_instrument_state(**self.ret_val)
        return self.ret_val

    def cmd(self, command, **kwargs):
        """ Converts a command (from Button press) into a Message that is sent to SP!Pr (if outbox is connected)

        :param command: The name of the function
        :param kwargs: Additional keyword arguments
        :return: None
        """

        if not self.outbox:
            return
        my_message = \
            msgm.Message.build_from_args(
                _type="CMD",
                _sender="SP",
                _recip="Pr",
                _contents="ib.Read_Operation",
                _data=oprtn.Operation.build_from_args(
                    func=command,
                    agent="Pr",
                    kwargs=kwargs
                ).package()
            )
        self.outbox.enqueue(my_message)


if __name__ == "__main__":
    tk_root = tk.Tk()
    from mcn_queues import MCNPriorityQueue

    spark_state = {
        "tray_position": "In",
        "thermostat": None,
        "encumbered": False,
        "expose": False,
    }
    pipes = MCNPriorityQueue()
    _safety = {'state': None}

    my_ui = SparkStateUI(tk_root, None, pipes, _safety)

    print(spark_state)
    spark_state.update(my_ui.run())
    print(spark_state)
    [m.print() for m in pipes.queue]
