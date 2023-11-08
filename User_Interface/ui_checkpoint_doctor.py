""" Checkpoint Doctor UI
@author: Ben C
"""

import tkinter as tk
from gui_constants import *
from mcn_status import Checkpoint


class Defibrillator:
    """ A GUI for reviving a Checkpoint """
    def __init__(self, root, task_key, checkpoint: Checkpoint, memory=None):
        """ Creates a basic input GUI to allow a user to create a fully specified Checkpoint

        :param root: a tkinter root
        :param task_key: The Checkpoints key value
        :param checkpoint: The Checkpoint
        :param memory: To pass in any retained information about the Checkpoint
        """
        # completion, location, level, data, queue, step, operation, task_key
        categories = ['Recipient', 'Checkpoint Key', 'Queue Name', 'Step #', 'Update Interval', 'Timeout']
        flag = False
        if memory:
            defaults = memory
        else:
            if checkpoint.queue == "Internal":
                flag = True
            defaults = [checkpoint.location, task_key, checkpoint.queue, checkpoint.step, "15", "None"]
        if len(categories) != len(defaults):
            raise ValueError(f"Lengths of categories {len(categories)} and defaults {len(defaults)} must be equal")

        self.root = root
        self.ret_val = defaults

        frame = tk.Frame(root)
        frame.winfo_toplevel().title('Checkpoint Doctor')
        frame.pack(side=tk.LEFT, fill=tk.BOTH, padx=10)
        prompt = tk.Label(frame, text="Please verify the following information", width=30)
        prompt.pack(**TB33)

        self.entries = {k: tk.StringVar(frame, value=defaults[i]) for i, k in enumerate(categories)}
        for k, entry in self.entries.items():
            tk.Label(frame, text=k, anchor=tk.W).pack(**TB33)
            tk.Entry(frame, textvariable=entry).pack(**TB33)

        tk.Button(frame, text="Submit", command=self.submit).pack(**TB33)
        tk.Button(frame, text="Cancel", command=self.cancel).pack(**TB33)

        self.info_box = tk.Label(frame, text="")
        self.info_box.pack(**TB33)
        if flag:
            self.info_box['text'] = f"Checkpoint is for Internal Operation\nResuscitation may fail"

    def submit(self):
        """ Compiles the input and commits them to the GUI's ret_val property

        :return: None
        """
        ret_val = list()
        for k, entry in self.entries.items():
            raw_entry = entry.get()
            if raw_entry == "":
                self.info_box['text'] = f"Missing {k}"
                return
            if raw_entry.lower() in ["none", "null"]:
                raw_entry = None
            ret_val += [raw_entry, ]

        self.ret_val = ret_val
        self.root.destroy()

    def cancel(self):
        """ Sets the return value to None and closes the GUI

        :return: None
        """
        self.ret_val = None
        self.root.destroy()

    def run(self):
        """ Mainloop for the GUI

        :return: [location, checkpoint_key, q_name, q_step#, update_interval(, timeout)] or None
        """
        self.root.wait_window(self.root)
        return self.ret_val


if __name__ == '__main__':
    my_gui = Defibrillator(tk.Tk(),
                           "abc123",
                           Checkpoint(None, "Loc", "Level", "info", "Internal", "step"),
                           memory=None)
    print(my_gui.run())
