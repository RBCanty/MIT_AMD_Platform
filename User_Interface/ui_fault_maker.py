""" Fault Creation/Reporting UI
@author: Ben C
"""

import tkinter

import mcn_status as mcs


class FaultMaker:
    """ A GUI for the Creation of Faults """
    def __init__(self, master):
        """ Creates a GUI for the creation and reporting of Fault objects

        :param master: a tkinter root
        """
        self.core = master
        self.ret_value = None

        all_keys = mcs.MCN_CFG[mcs.S_ALL_SYS]
        all_levels = list(mcs.LEVEL_MAP.keys())

        self.location = tkinter.StringVar()
        self.level = tkinter.StringVar()

        self.frame = tkinter.Frame(self.core)
        self.frame.winfo_toplevel().title("Fault Maker")
        self.frame.pack(side=tkinter.LEFT, fill=tkinter.BOTH, padx=10)

        tkinter.Label(self.frame, text="Location: ").grid(row=0, column=0, sticky=tkinter.W)
        tkinter.OptionMenu(self.frame, self.location, *all_keys).grid(row=0, column=1, sticky=tkinter.W)

        tkinter.Label(self.frame, text="Level: ").grid(row=1, column=0, sticky=tkinter.W)
        tkinter.OptionMenu(self.frame, self.level, *all_levels).grid(row=1, column=1, sticky=tkinter.W)

        tkinter.Label(self.frame, text="Info: ").grid(row=2, column=0, sticky=tkinter.W)
        self.info_entry = tkinter.Entry(self.frame, width=30)
        self.info_entry.grid(row=2, column=1, sticky=tkinter.W)

        tkinter.Button(self.frame, text="Create", command=self.report_fault, width=10).grid(row=3, column=0)
        tkinter.Button(self.frame, text="Cancel", command=self.core.destroy, width=10).grid(row=3, column=1)

    def report_fault(self):
        """ Validates and commits a user-reported Fault

        Closes on valid commit (remains open if unsuccessful)

        :return:
        """
        loc = self.location.get()
        level = self.level.get()
        info = self.info_entry.get()

        if not loc or not level or not info:
            return

        self.ret_value = mcs.Fault(loc, level, info)

        self.core.destroy()

    def wait(self):
        """ Mainloop of the FaultMaker

        Creates a blocking GUI element

        :return: the Fault object (or None if window exited via X or cancel)
        """
        self.core.wait_window(self.core)
        return self.ret_value


if __name__ == '__main__':
    tk = tkinter.Tk()
    my_ui = FaultMaker(tk)
    result = my_ui.wait()
    print(result)
