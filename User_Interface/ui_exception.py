""" Plate Location Reconciliation UI
@author: Ben C
"""

import tkinter as tk
from gui_constants import TB33, LB33


class ReconcilePlateUI:
    """
    A quick UI that requests a human user move a plate
    """

    def __init__(self, root: tk.Tk, plate_name: str, plate_type: str, direction: int):
        """ Creates a GUI to prompt the User on whether a manual reconciliation of a plate location was successful

        :param root: a tkinter root
        :param plate_name: the name of the plate
        :param plate_type: the plate type
        :param direction: (-1) plate is leaving the LPX, (+1) plate is entering the LPX
        """
        direction = int(direction)
        if direction == -1:  # -1 for accessing a plate (plate leaves LPX)
            dialog = f"Please move plate {plate_name} ({plate_type}) to the transfer station"
        else:
            raise ValueError(f"Only direction=-1 has implementable user-recovery (direction={direction})")

        self.root = root
        self.ret_val = None

        frame = tk.Frame(root)
        frame.winfo_toplevel().title("Plate Reconciliation")
        root.attributes('-topmost', 1)
        frame.pack(side=tk.LEFT, fill=tk.BOTH, padx=10)
        prompt = tk.Label(frame, text=dialog)
        prompt.pack(**TB33)

        button_frame = tk.Frame(frame)
        button_frame.pack(**TB33)
        tk.Button(button_frame, text="Recovery Complete", command=self.ok).pack(**LB33)
        tk.Button(button_frame, text="Unable to Recover", command=self.cancel).pack(**LB33)

    def ok(self, *_):
        self.ret_val = True
        self.root.destroy()

    def cancel(self, *_):
        self.ret_val = False
        self.root.destroy()

    def run(self) -> bool:
        """ Main loop of GUI, generates prompt for user and awaits the response

        :return: True (ok), False (cancel), or None (not set)
        """
        self.root.mainloop()
        return self.ret_val


class LPXPopupDialog(LookupError):
    """ Exception for when a plate's location gets lost in the updating state and the LPX can't resolve transfers

    The Popup GUI is linked to an exception such that it can be raised and interrupt operation until handled.
    """
    def __init__(self, *args, **_):
        """ Constructor for Exception (for Popup constructor, see run() method) """
        super(LPXPopupDialog, self).__init__(*args)

    @staticmethod
    def run(plate_name: str, plate_type: str, direction: int):
        """ Creates popup from exception

        :param plate_name: The name of the plate being moved
        :param plate_type: The type of the plate
        :param direction: (only -1 is supported: "Move to transfer station")
        :return: True - Recovery success, False - Recovery failed, None - Recovery aborted
        """
        prompt = ReconcilePlateUI(tk.Tk(), plate_name=plate_name, plate_type=plate_type, direction=direction)
        ret_val = prompt.run()
        return ret_val


if __name__ == '__main__':
    codex = {True: "Success",
             False: "Explcit Fail",
             None: "Implicit Fail"}
    try:
        raise LPXPopupDialog
    except LPXPopupDialog as popup:
        result = popup.run("test_plate", "test_type", -1)
        print(codex.get(result, "???"))
        if not result:
            raise RuntimeError("Process Failed")
        print("Carrying on")
