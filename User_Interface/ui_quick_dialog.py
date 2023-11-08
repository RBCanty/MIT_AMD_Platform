""" A quick and simple GUI prompt for asking basic questions
@author: Ben C
"""

import tkinter as tk
from gui_constants import *
from typing import Union


class QuickUI:
    """ A quick UI generator so you can create a dialog without needing to know tkinter

    Good if all input is Button or single-entry based.  For selecting between multiple options,
    see :class:`QuickSelectUI`
    """
    def __init__(self, root: Union[tk.Tk, tk.Toplevel],
                 *, title: str, dialog: str, buttons: dict, has_entry=False, ret_if_ok=""):
        """
        Creates a simple UI

        By default, all UI will have an OK button (you may overwrite it by adding your own {"OK": func} to the buttons
        input argument).  When has_entry is True, a Submit button is added by default and the return of run() will be
        the text in the entry field (the return is "" otherwise).  When has_entry is True, ALL functions in the buttons
        dictionary will receive the text of the entry widget as their args[0].  All buttons will close the message box.

        :param root: A Tkinter object (should be tk.Tk() from the calling function, tkinter won't let me spawn it here)
        :param title: The text shown on the header bar of the message box
        :param dialog: The text shown in the message box
        :param buttons: A dictionary of {"Button text": pyhton function handle}
        :param has_entry: (Default: False) True - adds a text entry field
        :param ret_if_ok: When using the default OK, the return value of run() will be "" or the value specified by
        ret_if_ok (note that only OK will generate this return value)
        """
        if title is None:
            title = "User Dialog"
        if dialog is None:
            dialog = "Default Message"
        if buttons is None:
            buttons = {"OK": self.ok}
        if "OK" not in buttons.keys() and not has_entry:
            buttons.update({"OK": self.ok})
        self.root = root
        self.has_entry = has_entry
        self.ret_val = ""
        self.ret_if_ok = ret_if_ok

        frame = tk.Frame(root)
        frame.winfo_toplevel().title(title)
        root.attributes('-topmost', 1)
        frame.pack(side=tk.LEFT, fill=tk.BOTH, padx=10)
        prompt = tk.Label(frame, text=dialog)
        prompt.pack(**TB33)

        if self.has_entry:
            buttons.update({"Submit": self.ok})
            self.entry_field = tk.Entry(frame)
            self.entry_field.pack(**TB33)

        button_frame = tk.Frame(frame)
        button_frame.pack(**TB33)

        for t, f in buttons.items():
            tk.Button(button_frame, text=t, command=lambda x=f: self.func_wrapper(x)).pack(**LB33)

    def ok(self, *_):
        """ the OK button

        If the GUI has an entry field, causes exit

        Otherwise, if sets the return value to the ret_if_ok (from constructor) if it is specified, otherwise leaves
        the return value unaltered (default: "")

        :return: None
        """
        if not self.has_entry:
            print("User Clicked OK")
            self.ret_val = self.ret_if_ok if self.ret_if_ok else self.ret_val

    def func_wrapper(self, f):
        """ Wraps functions calls for buttons such that they are given the contents of the Entry box (if present) and
        will exit the popup upon execution.  Also sets the return value to the contents of the Entry box (if present)

        :param f:
        :return:
        """
        if self.has_entry:
            self.ret_val = self.entry_field.get()
            f(self.ret_val)
        else:
            f()
        self.root.destroy()
        self.root.quit()

    def run(self) -> str:
        """ Executes the popup

        IF there is an Entry box:
          - The return will be the contents of said Entry

        IF there is no Entry box:
          - If the Popup was closed by pressing "OK", then the return is the value of ret_if_ok (specified in
            constructor and defaults to "")
          - Otherwise, returns ""

        :return: (See above)
        """
        self.root.mainloop()
        return self.ret_val


class QuickSelectUI:
    """ A simple Helper UI for making selections rather than giving inputs or pushing buttons
    (for that see :class:`QuickUI`) that doesn't require knowing how tkinter works
    """
    def __init__(self, root: Union[tk.Tk, tk.Toplevel], *, title: str, dialog: str, options: list, default=None):
        """ Creates a selector popup box

        :param root: a tkinter root
        :param title: the title of the popup
        :param dialog: the prompt
        :param options: a list of options available
        :param default: the default selection (added to options, if not already present)
        """

        if title is None:
            title = "User Dialog"
        if dialog is None:
            dialog = "Default Message"
        if options is None:
            options = list()
        if default is None:
            default = "None"
        if default not in options:
            options = [default, ] + options
        self.root = root
        self.ret_val = default

        frame = tk.Frame(root)
        frame.winfo_toplevel().title(title)
        frame.pack(side=tk.LEFT, fill=tk.BOTH, padx=10)
        prompt = tk.Label(frame, text=dialog)
        prompt.pack(**TB33)

        self.selection = tk.StringVar(frame, value=default)

        menu = tk.OptionMenu(frame, self.selection, *options)
        menu.pack(**TB33)

        button = tk.Button(frame, text="Submit", command=self.submit)
        button.pack(**TB33)

    def submit(self):
        """ Sets the return value to the selection and exits

        :return: None
        """
        self.ret_val = self.selection.get()
        self.root.destroy()
        self.root.quit()

    def run(self):
        """ Executes the popup

        :return: the value of the selection
        """
        self.root.mainloop()
        return self.ret_val


if __name__ == "__main__":

    def my_func(*args):
        print("I called my custom function")
        if args:
            print(f"User had entered: {args[0]}")

    my_list = [1, 2, 3]
    print(my_list)

    def add_four(*_):
        my_list.append(4)

    my_prompt = QuickUI(tk.Tk(),
                        title="My Title",
                        dialog="My custom dialog prompt",
                        buttons={"Do Something": my_func, "add 4": add_four},
                        has_entry=True,
                        ret_if_ok="user did click OK")
    ret = my_prompt.run()

    print(my_list)
    print(ret)

    prompt2 = QuickSelectUI(tk.Tk(),
                            title="My Selector",
                            dialog="Please select one of the following",
                            options=['Steak', 'Fish'])

    print(prompt2.run())
