""" GUI for handling the Liquid Handler
@author: Ben C
"""

import tkinter as tk
from gui_constants import *
import queue
import json


class EvowareHandler:
    """ A GUI interface to the Liquid Handler communication queue """
    def __init__(self, master, child_com, safety=None):
        """ Creates GUI to input messages to or reset the System-LH communication pipe

        :param master: a tkinter root
        :param child_com: For communication with parent System
        :param safety: For use with System-level GUI (module_GUI.py) to prevent multiple copies from being spawned
        """
        self.core = master
        self.core.protocol("WM_DELETE_WINDOW", self.on_closing)
        self.child_com = child_com
        self.safety = safety
        if self.safety:
            self.safety['state'] = tk.DISABLED
        self.tk = dict()
        self.reg = dict()

        self.frame = tk.Frame(self.core)
        tk.Button(self.frame, text="Reset", command=self.reset)                    .grid(row=0, column=0, **GRID_B33)
        tk.Label(self.frame, text="Clears Lh comm. queue")                         .grid(row=0, column=1, **GRID_B33)
        tk.Button(self.frame, text="Stop", command=self.stop)                      .grid(row=1, column=0, **GRID_B33)
        tk.Label(self.frame, text="Sends stop signal to Lh comm. queue")           .grid(row=1, column=1, **GRID_B33)
        self.debug = tk.Entry(self.frame, width=32)
        self.debug                                                                 .grid(row=2, column=0, **GRID_B33)
        tk.Button(self.frame, text="Submit to Lh comm. queue", command=self.submit).grid(row=2, column=1, **GRID_B33)

        self.frame.pack()

    def on_closing(self, *_):
        """ If connected to larger GUI, releases button

        :return:
        """
        self.core.destroy()
        if self.safety:
            self.safety['state'] = tk.NORMAL

    def reset(self):
        """ Resets/Clears the communications pipe

        :return:
        """
        with self.child_com as cc:
            cc['Lh'] = queue.Queue()
            cc['Lh'].put(['IDLE', ])
        self.on_closing()

    def stop(self):
        """ Puts an emergency stop signal into the communications pipe

        :return:
        """
        with self.child_com as cc:
            cc['Lh'].queue.popleft()
            cc['Lh'].put(['EMERGENCYSTOP', ''])
        self.on_closing()

    def submit(self):
        """ Puts a message into the communications pipe

        If the message begins with "json:" the rest of the string will be parsed by json and submitted directly to
        the communications pipe

        Otherwise, the submission is interpreted as a string, wrapped by a list, then submitted to the communications
        pipe.

        :return:
        """
        submission = self.debug.get()
        if not submission:
            return
        if 'json:' == submission[:5]:
            try:
                signal = json.loads(submission[5:])
            except json.decoder.JSONDecodeError as jde:
                self.debug.select_range(jde.pos + 5, jde.pos + 6)
                return
        else:
            signal = [submission, ]
        with self.child_com as cc:
            cc['Lh'] = queue.Queue()
            cc['Lh'].put(signal)
        self.on_closing()


if __name__ == '__main__':
    from thread_safe_data import ThreadSafeDataContainer
    pipe = ThreadSafeDataContainer({'Lh': queue.Queue()})
    with pipe as _p:
        _p['Lh'].put(['Original queue', ])

    tk_root = tk.Tk()

    EvowareHandler(tk_root, pipe)

    tk_root.mainloop()

    with pipe as _p:
        print(_p['Lh'].queue)
