""" A GUI for looking at and modifying the plates in the LPX
@author: Ben C
"""

from __future__ import annotations  # may no longer be needed, is typically used to change the order in which typing
# hints are processed, the classic example is "class A: def method() --> A: ..." will throw an error without this import
# because "A" hasn't been defined yet.  So this, along with the TYPE_CHECKING below, may be for the fact that the
# storage shell GUI references the LPX Controller and vice versa.

import tkinter as tk
from ttkwidgets.autocomplete import AutocompleteCombobox
from gui_constants import *

# Does nothing at Runtime, used to allow IDE to type-hint without circular import error
from typing import TYPE_CHECKING, Dict, Union
if TYPE_CHECKING:
    from AMD_LPX_Control import LpxController

ALLOWED_PLATETYPES = sorted(["96 Well Microplate", "96 Well Microplate PS", "96 Well Microplate Half Area",
                             "96 Well Filtration Plate", "96 Well PCR Plate",
                             "96 Well DeepWell", "54 Well Microplate", "DiTi SBS Waste", ])


class StorageShellUserInterface:
    """ A GUI which displays the contents of the Hotels """
    def __init__(self, root: tk.Tk, lpx_controller: LpxController):
        """ Creates a UI for viewing and modifying the contents of the LPX on a hotel basis.

        :param root: a tkinter root
        :param lpx_controller: the controller for the LPX
        """
        self.core = root
        master_frame = tk.Frame(self.core)
        master_frame.winfo_toplevel().title("LPX Storage Shell")
        master_frame.pack(side=tk.LEFT, fill=tk.BOTH, padx=10, expand=True)
        self.lpx_controller = lpx_controller
        self.mutables: Dict[str, Union[tk.Entry, tk.Button, dict]] = dict()
        self.ret_val = False

        # Part 1: Selector
        self.selector_frame = tk.LabelFrame(master_frame, text="Select")
        self.selector_frame.pack(**TB33)
        self.build_selector_frame()
        # Part 2: Fields
        self.slots_frame = tk.Frame(master_frame)
        self.slots_frame.pack(**TB33, expand=True)
        # Part 3: Button
        self.mutables['submit button'] = tk.Button(master_frame, text="Submit", command=self.submit, state=tk.DISABLED)
        self.mutables['submit button'].pack(**TB33)
        self.mutables['reset button'] = tk.Button(master_frame, text="Reset Form", command=self.reset)
        self.mutables['reset button'].pack(**TB33)

    def build_selector_frame(self):
        """ A header frame where the hotel number is provided such that the slots frame can be built

        :return:
        """
        prompt = tk.Frame(self.selector_frame)
        info = tk.Frame(self.selector_frame)
        prompt_label = tk.Label(prompt, text="Hotel: ")
        self.mutables['hotel selector'] = tk.Entry(prompt)
        self.mutables['hotel selector'].insert(tk.END, '0')
        info_details = tk.Label(info, text="Transfer Station: 0, Bottom: 1-10, Top: 11-20")
        self.mutables['search button'] = tk.Button(info, text="Search", command=self.build_slots_frame)
        prompt.pack(**TB33)
        info.pack(**TB33)
        prompt_label.pack(**LB33)
        self.mutables['hotel selector'].pack(**LB33)
        info_details.pack(**TB33)
        self.mutables['search button'].pack(**TB33)

    def build_slots_frame(self):
        """ Constructs a representation of the selected hotel for viewing and modification

        :return:
        """
        self.mutables['hotel selector']['state'] = tk.DISABLED

        for item in self.slots_frame.winfo_children():
            item.destroy()
        selected_hotel = self.mutables['hotel selector'].get()
        selected_hotel = selected_hotel if selected_hotel else '0'
        selected_hotel, _ = self.lpx_controller.validate_location(selected_hotel, 0)
        if selected_hotel == "updating":
            return
        size = self.lpx_controller.configuration[selected_hotel]['size']

        # create headers
        tk.Label(self.slots_frame, text="Slot").grid(row=0, column=0)
        tk.Label(self.slots_frame, text="Plate Name").grid(row=0, column=1)
        tk.Label(self.slots_frame, text="Plate Type").grid(row=0, column=2)

        self.slots_frame.columnconfigure(1, weight=3)
        self.slots_frame.columnconfigure(2, weight=1)

        # create the digital representation of the Hotel
        self.mutables['slot list'] = dict()
        for i in reversed(range(1, size+1)):
            row_index = size + 1 - i
            tk.Label(self.slots_frame, text=f"{i} ").grid(row=row_index, column=0)
            self.mutables['slot list'][i] = [tk.Entry(self.slots_frame),
                                             AutocompleteCombobox(self.slots_frame, completevalues=ALLOWED_PLATETYPES)]
            self.mutables['slot list'][i][0].grid(row=row_index, column=1, sticky='ew')
            self.mutables['slot list'][i][1].grid(row=row_index, column=2, sticky='ew')

        pm = self.lpx_controller.get_plates_from_db('partial_matches')
        for k in pm.keys():
            hotel, slot = pm[k]['location'][1]
            if hotel == selected_hotel:
                self.mutables['slot list'][slot][0].insert(tk.END, pm[k]['container_name'])
                self.mutables['slot list'][slot][1].insert(tk.END, pm[k]['labware_type'])

        self.mutables['submit button']['state'] = tk.ACTIVE

    def submit(self):
        """ Commit changes to the database (the LPX does not have a local copy, it uses the database)

        :return: None
        """
        manifest = list()  # [[name, type, [hotel, slot]], ...]
        for slot in self.mutables['slot list'].keys():
            item = [self.mutables['slot list'][slot][0].get(),
                    self.mutables['slot list'][slot][1].get(),
                    [self.mutables['hotel selector'].get(),
                    slot]]
            if item[0] != '' and item[1] != '':
                manifest.append(item)
        self.ret_val += self.lpx_controller.carousel_update_from_user(manifest, self.mutables['hotel selector'].get())
        self.build_slots_frame()
        self.mutables['submit button']['state'] = tk.DISABLED
        pass

    def reset(self):
        """ Clears the Hotel selection input and clears the display for the positions of the hotel.

        :return:
        """
        self.mutables['hotel selector'].delete(*INPUT_CONTENTS)
        self.mutables['hotel selector']['state'] = tk.NORMAL
        for item in self.slots_frame.winfo_children():
            item.destroy()

    def run(self):
        """ Runs the UI

        :return: True if a change has occurred, False otherwise
        """
        self.core.mainloop()
        return bool(self.ret_val)


if __name__ == '__main__':
    from AMD_LPX_Control import LpxController

    _root = tk.Tk()
    ui = StorageShellUserInterface(_root, LpxController(None, r'..\Shell_SP\LPX_API\LPX.cfg', real_mode=False))
    ui.run()
