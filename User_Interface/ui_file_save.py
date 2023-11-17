""" A dialogue box for prompting about files
@author: Ben C
"""

import tkinter as tk
from tkinter import filedialog


class FileSaver:
    """ Opens a GUI which can spawn an OS-level explorer for finding files on a local machine """
    def __init__(self, root: tk.Tk, default_filename=None):
        """ Creates a small dialogue box that prompts the user to select a file using an OS-level file explorer.

        :param root: a tkinter root
        :param default_filename: The filepath to be returned if nothing is selected
        """
        self.core = root
        self.ret_val = default_filename
        frame = tk.Frame(self.core)
        frame.winfo_toplevel().title("File")
        frame.pack(side=tk.LEFT, fill=tk.BOTH, padx=10)
        prompt = tk.Label(frame, text="File received, please specify destination")
        prompt.pack(side=tk.TOP)
        save_as_button = tk.Button(frame, text="Save as", command=self.file_select)
        save_as_button.pack(side=tk.TOP)

    def file_select(self):
        """ Spawns the OS-level file explorer, updates the return value from what the explorer selected,
        closes after selection.

        :return:
        """
        self.ret_val = filedialog.asksaveasfilename(title='Save a file', initialdir='/')
        self.core.destroy()

    def run(self):
        """ Runs the manager UI

        :return: A file path
        """
        self.core.mainloop()
        return self.ret_val
