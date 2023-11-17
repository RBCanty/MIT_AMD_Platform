""" GUI for editing HPLC solvent volumes
@author: Ben C
Based on HPLC code provided by Matt McDonald
"""

import tkinter as tk
from tkfire import TkFire, LAYOUT, TYPE, CHILDREN, grid_arg
from tkinter import filedialog
from gui_constants import *
import database_interface as dbi
from os import path
from typing import List, Union


def is_valid_volume(val):
    """ Helper method for determining if a number is both numerical and positive

    :param val: A volume
    :return: True - is valid, False - otherwise
    """
    try:
        float(val)
    except ValueError:
        return False
    else:
        return float(val) >= 0


class VolumeSet:
    """ A GUI for setting Solvent Volumes """
    def __init__(self, root: tk.Tk, title, dialog,
                 categories: List[str], defaults: List[Union[int, float]] = None):
        """ Creates a Dialogue box prompting for volumes for a given set of solvents

        :param root: a tkinter root
        :param title: The title of the dialog box
        :param dialog: The prompt
        :param categories: List of solvents
        :param defaults: (Optional, defaults to all "0") List of initial/default solvent volumes
        """
        if title is None:
            title = "User Dialog"
        if dialog is None:
            dialog = "Default Message"
        if categories is None:
            categories = list()
        if defaults is None:
            defaults = [0] * len(categories)
        defaults = [float(x) for x in defaults]

        if len(categories) != len(defaults):
            raise ValueError(f"Lengths of categories {len(categories)} and defaults {len(defaults)} must be equal")

        self.root = root
        self.ret_val = defaults

        frame = tk.Frame(root)
        frame.winfo_toplevel().title(title)
        frame.pack(side=tk.LEFT, fill=tk.BOTH, padx=10)
        prompt = tk.Label(frame, text=dialog, width=30)
        prompt.pack(**TB33)

        # Create (prompt: Label, input: Entry) entries in the GUI for each solvent
        self.entries = {k: tk.StringVar(frame, value=str(defaults[i])) for i, k in enumerate(categories)}
        for k, entry in self.entries.items():
            tk.Label(frame, text=k, anchor=tk.W).pack(**TB33)
            tk.Entry(frame, textvariable=entry).pack(**TB33)

        tk.Button(frame, text="Submit", command=self.submit).pack(**TB33)

        self.info_box = tk.Label(frame, text="")
        self.info_box.pack(**TB33)

    def submit(self):
        """ Validates all volumes, and if valid, sets the return value and exits the GUI

        :return:
        """
        raw_entries: List[str] = list()
        for k, entry in self.entries.items():
            raw_entry = entry.get()
            if raw_entry == "":
                self.info_box['text'] = f"Missing {k}"
                return
            if not is_valid_volume(raw_entry):
                self.info_box['text'] = f"{k} Volume '{raw_entry}' is not a valid volume"
                return
            raw_entries.append(raw_entry)

        self.ret_val = [float(x) for x in raw_entries]
        self.root.destroy()

    def run(self):
        """ Mainloop of UI

        :return: A list of volumes (float) in the same order as the solvents
        """
        self.root.mainloop()
        return self.ret_val


class EllseeManager:
    """ More in-depth GUI for the HPLC """

    def __init__(self, master, shimadzu_controller, safety=None):
        """ Creates a GUI for monitoring the solvent volumes, column, current wellplate, current folder, and view batch
        information.

        :param master: a tkinter root
        :param shimadzu_controller: A controller for the HPLC
        :param safety: For the System-level GUI (module_GUI.py) to manage buttons such that only one such controller can
          exist at a time.
        """
        self.core = master
        self.core.protocol("WM_DELETE_WINDOW", self.on_closing)
        self.shimadzu_controller = shimadzu_controller
        self.solvents = ['Aqueous', 'Organic', 'MS_makeup']
        self.sel_column = tk.StringVar(value=self.shimadzu_controller.lc_column)
        self.safety = safety
        if self.safety:
            self.safety['state'] = tk.DISABLED

        init_column = self.shimadzu_controller.lc_column

        self.memory = {'sel_column': tk.StringVar,
                       'column_options': ["analytical", "semiprep", "UPLC"],
                       'cmd_submit': self.submit,
                       'cmd_close': self.close,
                       'cmd_set_wellplate': self.set_wellplate,
                       'cmd_find_folder': self.sel_folder,
                       'cmd_set_folder': self.set_folder,
                       'cmd_get_status': self.get_status,
                       'cmd_stop_batch': self.stop_batch,
                       'cmd_resume_batch': self.resume_batch
                       }

        mother = {
            'mf': {
                TYPE: ['Frame', {}],
                LAYOUT: ['pack', {}],
                CHILDREN: {
                    'left': {
                        TYPE: ['Frame', {}],
                        LAYOUT: ['pack', LBA3],
                        CHILDREN: {
                            'volumes': {
                                TYPE: ['LabelFrame', {'text': "Solvent Volumes"}],
                                LAYOUT: ['pack', TB33],
                                CHILDREN: {
                                    'header_solvent': {
                                        TYPE: ['Label', {'text': "Solvent"}],
                                        LAYOUT: ['grid', grid_arg(0, 0)]
                                    },
                                    'header_volume': {
                                        TYPE: ['Label', {'text': "New Volume (uL)"}],
                                        LAYOUT: ['grid', grid_arg(0, 1)]
                                    },
                                    'header_previous': {
                                        TYPE: ['Label', {'text': "Previous"}],
                                        LAYOUT: ['grid', grid_arg(0, 2)]
                                    },
                                    # Defines the Label-Entry-Label triplet for solvents
                                    # There is probably a way to combine the comprehensions
                                    **{
                                        f"{solvent}_label": {
                                            TYPE: ['Label', {'text': solvent}],
                                            LAYOUT: ['grid', grid_arg(ri, 0)]
                                        }
                                        for ri, solvent in enumerate(self.solvents, 1)
                                    },
                                    **{
                                        f"{solvent}_entry": {
                                            TYPE: ['Entry', {'justify': tk.RIGHT}],  # initialize later
                                            LAYOUT: ['grid', grid_arg(ri, 1)]
                                        }
                                        for ri, solvent in enumerate(self.solvents, 1)
                                    },
                                    **{
                                        f"{solvent}_prev": {
                                            TYPE: ['Label', {'text': "-1.0"}],
                                            LAYOUT: ['grid', grid_arg(ri, 2)]
                                        }
                                        for ri, solvent in enumerate(self.solvents, 1)
                                    },
                                }
                            },
                            'columns': {
                                TYPE: ['LabelFrame', {'text': "Column"}],
                                LAYOUT: ['pack', TB33],
                                CHILDREN: {
                                    'column_selector': {
                                        TYPE: ['OptionMenu', {'variable': ['sel_column', {'value': init_column}],
                                                              'values': 'column_options'}],
                                        LAYOUT: ['pack', LBA3]
                                    },
                                    'prev_column': {
                                        TYPE: ['Label', {'text': f"Prev: {init_column}"}],
                                        LAYOUT: ['pack', LB33]
                                    }
                                }
                            },
                            'buttons': {
                                TYPE: ['Frame', {}],
                                LAYOUT: ['pack', TB33],
                                CHILDREN: {
                                    'b_submit': {
                                        TYPE: ['Button', {'text': "Submit", 'command': 'cmd_submit'}],
                                        LAYOUT: ['pack', LBA3]
                                    },
                                    'b_close': {
                                        TYPE: ['Button', {'text': "Close", 'command': 'cmd_close'}],
                                        LAYOUT: ['pack', LBA3]
                                    }
                                }
                            },
                            'info_box': {
                                TYPE: ['Text', {'height': 4, 'width': 32}],
                                LAYOUT: ['pack', LBA3]
                            }
                        }
                    },
                    'right': {
                        TYPE: ['Frame', {}],
                        LAYOUT: ['pack', LBA3],
                        CHILDREN: {
                            'wellplate_selector': {
                                TYPE: ['Frame', {}],
                                LAYOUT: ['pack', TB33],
                                CHILDREN: {
                                    'header': {
                                        TYPE: ['Label', {'text': "Current Wellplate:"}],
                                        LAYOUT: ['pack', TB33]
                                    },
                                    'e_wellplate': {
                                        TYPE: ['Entry', {}],
                                        LAYOUT: ['pack', TB33]
                                    },  # initialize later
                                    'b_select': {
                                        TYPE: ['Button', {'text': "Set", 'command': 'cmd_set_wellplate'}],
                                        LAYOUT: ['pack', TB33]
                                    }
                                }
                            },
                            'folder_selector': {
                                TYPE: ['Frame', {}],
                                LAYOUT: ['pack', TB33],
                                CHILDREN: {
                                    'header': {
                                        TYPE: ['Label', {'text': "Current Folder:"}],
                                        LAYOUT: ['pack', TB33]
                                    },
                                    'e_project': {
                                        TYPE: ['Entry', {}],
                                        LAYOUT: ['pack', TB33]
                                    },  # initialize later
                                    'b_find': {
                                        TYPE: ['Button', {'text': "Find", 'command': "cmd_find_folder"}],
                                        LAYOUT: ['pack', TB33]
                                    },
                                    'b_set': {
                                        TYPE: ['Button', {'text': "Set", 'command': "cmd_set_folder"}],
                                        LAYOUT: ['pack', TB33]
                                    },
                                }
                            },
                            'info_box_frame': {
                                TYPE: ['LabelFrame', {'text': "Info"}],
                                LAYOUT: ['pack', TB33],
                                CHILDREN: {
                                    'info_box': {
                                        TYPE: ['Text', {'height': 4, 'width': 24}],
                                        LAYOUT: ['pack', TB33]
                                    }
                                }
                            },
                            'buttons': {
                                TYPE: ['Frame', {}],
                                LAYOUT: ['pack', TB33],
                                CHILDREN: {
                                    'b_get_status': {
                                        TYPE: ['Button', {'text': "Get Status", 'command': "cmd_get_status"}],
                                        LAYOUT: ['pack', LBA3]
                                    },
                                    'b_batch_stop': {
                                        TYPE: ['Button', {'text': "Batch Stop", 'command': "cmd_stop_batch"}],
                                        LAYOUT: ['pack', LBA3]
                                    },
                                    'b_batch_resume': {
                                        TYPE: ['Button', {'text': "Batch Resume", 'command': "cmd_resume_batch"}],
                                        LAYOUT: ['pack', LBA3]
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        self.ui = TkFire(self.core, self.memory, mother)
        self.gui = self.ui.gui
        self.gui['mf'].winfo_toplevel().title("HPLC Manager")

        self.reset_volume_and_column()
        self.initialize_wellplate_and_project()

    def reset_volume_and_column(self):
        """

        :return:
        """
        for solvent in self.solvents:
            volume = self.shimadzu_controller.solvent_levels.get(solvent, -1)
            self.gui[f'mf!left!volumes!{solvent}_entry'].delete(0, tk.END)
            self.gui[f'mf!left!volumes!{solvent}_entry'].insert(0, volume)
            self.gui[f'mf!left!volumes!{solvent}_prev']['text'] = str(volume)
        self.gui['mf!left!columns!prev_column']['text'] = f"Prev: {self.shimadzu_controller.lc_column}"

    def initialize_wellplate_and_project(self):
        """

        :return:
        """
        wellplate_name = self.shimadzu_controller.wellplate.get('container_name', '<Not Found>')
        self.gui['mf!right!wellplate_selector!e_wellplate'].insert(0, wellplate_name)
        self.gui['mf!right!folder_selector!e_project'].insert(0, self.shimadzu_controller.folderName)

    def on_closing(self, *_):
        """ Used to note the closing of the window to release the button of the main GUI if present

        :return:
        """
        self.ui.destroy()
        if self.safety:
            self.safety['state'] = tk.NORMAL

    def update_infobox(self, side, message):
        """ Posts messages to the infobox

        :param side: Which infobox to use
        :param message: The message to be printed
        :return: None
        """
        if isinstance(message, str):
            pass
        elif isinstance(message, (list, tuple)):
            message = "\n".join([str(s) for s in message])
        else:
            message = str(message)
        if side == 'left':
            gui_path = "mf!left!info_box"
        else:
            gui_path = "mf!right!info_box_frame!info_box"
        self.gui[gui_path].delete("1.0", tk.END)
        self.gui[gui_path].insert(tk.END, message)

    def close(self, *_):
        """ Calls self.on_closing()

        :return:
        """
        self.on_closing()

    def submit(self, *_):
        """ Posts updates for solvent and column (prints Left)

        :return:
        """
        for solvent in self.solvents:
            sol_value = self.gui[f'mf!left!volumes!{solvent}_entry'].get()
            if not is_valid_volume(sol_value):
                self.update_infobox('left', f"Bad {solvent}")
                return

        col_sel = self.memory['sel_column'].get()  # noqa:
        # the StrVar is initialized in TkFire, and so will have a get method at runtime

        if col_sel not in self.memory['column_options']:
            self.update_infobox('left', f"Bad column")
            return

        for solvent in self.solvents:
            sol_value = self.gui[f'mf!left!volumes!{solvent}_entry'].get()
            self.shimadzu_controller.solvent_levels[solvent] = float(sol_value)

        self.shimadzu_controller.lc_mode = col_sel
        self.shimadzu_controller.lc_column = col_sel

        self.update_infobox('left', f"Updated")
        self.reset_volume_and_column()

    def set_wellplate(self):
        """ Updates the wellplate field of the controller to a new document from the database (prints Right)

        :return:
        """
        wellplate_name = self.gui['mf!right!wellplate_selector!e_wellplate'].get()
        if not wellplate_name or wellplate_name == "<Not Found>":
            self.update_infobox('right', 'Wellplate value not populated')
            return
        try:
            doc, _, _ = dbi.query_document('LC', dbi.COL_WELLPLATES, dbi.DBG_CONTAINER_NAME, wellplate_name)
        except dbi.cexc.DatabaseRequestError:
            self.update_infobox('right', 'Wellplate not found in database')
            return
        self.update_infobox('right', 'Updated wellplate')
        self.shimadzu_controller.wellplate = doc

    def sel_folder(self):
        """ Prompts user for selecting the project directory the controller should use (prints Right)

        :return:
        """
        init = r"C:\Shimadzu\Data\Automated" if path.isdir(r"C:\Shimadzu\Data\Automated") else '/'
        folder = filedialog.askdirectory(title='Find LC Project Folder', initialdir=init)
        if folder:
            self.gui['mf!right!folder_selector!e_project'].delete(0, tk.END)
            self.gui['mf!right!folder_selector!e_project'].insert(0, folder)

    def set_folder(self):
        """ Updates the project folder used by the controller (prints Right)

        :return:
        """
        folder = self.gui['mf!right!folder_selector!e_project'].get()
        if path.isdir(folder):
            self.update_infobox('right', 'Project folder updated')
            self.shimadzu_controller.folderName = folder
        else:
            self.update_infobox('right', 'Failed to locate project folder')

    def get_status(self):
        """ Requests status information from the controller (prints Right)

        :return:
        """
        ret_val = self.shimadzu_controller.status_ping()
        self.update_infobox('right', ret_val)

    def stop_batch(self):
        """ Asks the controller to stop the batch (prints Right)

        :return:
        """
        ret_val = self.shimadzu_controller.lc_stop()
        if ret_val[0]:
            self.shimadzu_controller.running = (False, "Batch Stopped")
            self.shimadzu_controller.update_solvent_levels(pumps_off=True)
        self.update_infobox('right', ret_val)

    def resume_batch(self):
        """ Asks the controller to resume a batch (prints Right)

        :return:
        """
        ret_val = self.shimadzu_controller.resume_batch()
        self.update_infobox('right', ret_val)


class DummyLC:
    """ Minimal class for testing GUI elements """
    def __init__(self):
        """ Creates a minimally viable mimic for GUI testing (see: Shell_LC/Shimadzu_API.py)  """
        self.lc_mode = "analytical"
        self.lc_column = "analytical"
        self.solvent_levels = {'Aqueous': 0.0, 'Organic': 0.0,
                               'MS_makeup': 0.0, 'last_update_time': -1.0,
                               'last_sample': "A0"}
        self.wellplate = {}
        self.running = (False, "No Batch")
        self.folderName = ""

    @staticmethod
    def status_ping():
        return True, "Dummy status message"

    @staticmethod
    def lc_stop():
        return True, "Dummy batch stop message"

    def resume_batch(self):
        self.running = (True, "Batch Resumed")
        return True, "B_Need to rearrange fraction collector and possibly fetch new plate"

    @staticmethod
    def update_solvent_levels(pumps_off=True):
        return True, f"Pumps set to {pumps_off}"


if __name__ == '__main__':
    # # VolumeSet
    # my_gui = VolumeSet(tk.Tk(), title="my title", dialog="specify!",
    #                    categories=["Apples", "Oranges", "Pears"], defaults=None)
    # print(my_gui.run())
    # print("Test")

    # # EllseeManager
    my_controller = DummyLC()
    tk_root = tk.Tk()
    EllseeManager(tk_root, my_controller)
    tk_root.mainloop()

    print(my_controller.lc_mode)
    print(my_controller.lc_column)
    print(my_controller.solvent_levels)
    print(my_controller.running)
    print(my_controller.folderName)
    print(my_controller.wellplate)
