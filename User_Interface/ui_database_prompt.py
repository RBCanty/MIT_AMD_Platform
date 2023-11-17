""" A dialog box for getting database configuration details
@author: Ben C
"""

import tkinter as tk
from typing import Union, List, Tuple
from gui_constants import TB33
import database_interface as dbi


class DatabasePrompt:
    """ A quick GUI for getting Database-stye Network options """
    def __init__(self, root: Union[tk.Tk, tk.Toplevel], defaults: Union[List[str], Tuple[str, str]]):
        """
        Creates a simple UI with a prompt, two entries (url and port), and button

        :param root: A Tkinter object (should be tk.Tk() from the calling function, tkinter won't let me spawn it here)
        :param defaults: The default or initial values of the url and port
        """
        title = "Database Connection Error"
        dialog = "Please check Database server and connectivity options\n" \
                 "(E.g. Url: 'http://01.23.45.67:12345' & Port: '01234')"

        self.root = root
        if (defaults is None) or (not defaults):
            defaults = ("", "")
        self.ret_val = defaults

        frame = tk.Frame(root)
        frame.winfo_toplevel().title(title)
        root.attributes('-topmost', 1)
        frame.pack(side=tk.TOP, fill=tk.BOTH, padx=10, expand=True)

        prompt = tk.Label(frame, text=dialog)
        prompt.grid(column=0, row=0, columnspan=2, pady=6)
        tk.Label(frame, text="Url:  ", anchor="w").grid(column=0, row=1, sticky="W", pady=3)
        tk.Label(frame, text="Port: ", anchor="w").grid(column=0, row=2, sticky="W", pady=3)
        self.url_entry = tk.Entry(frame)
        self.url_entry.grid(column=1, row=1, sticky="EW", pady=3)
        self.port_entry = tk.Entry(frame)
        self.port_entry.grid(column=1, row=2, sticky="EW", pady=3)
        frame.columnconfigure(1, weight=1, minsize='42p')

        button_frame = tk.Frame(root)
        button_frame.pack(**TB33)
        tk.Button(button_frame, text="Confirm", command=self.ok).pack(**TB33)

        self.url_entry.insert(0, defaults[0])
        self.port_entry.insert(0, defaults[1])

    def ok(self, *_):
        """ Loads the entries into the return value and closes the GUI

        :return: None
        """
        self.ret_val = (self.url_entry.get(), self.port_entry.get())
        self.root.destroy()
        self.root.quit()

    def run(self) -> Tuple[str, str]:
        """ Executes the popup (blocking) and returns a tuple (or list)

        :return: (url, port)
        """
        self.root.mainloop()
        return self.ret_val


if __name__ == '__main__':
    new_database_settings = DatabasePrompt(
        tk.Tk(),
        defaults=[dbi.DATABASE_URL, dbi.DATABASE_PORT]
    ).run()

    print(f"Popup returned: '{new_database_settings}'")
