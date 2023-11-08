""" Fault arbitrator

Mostly for allowing the User to take the blame (mea culpa) for Faults in a centralized form

@author: Ben C
"""

import logging
import tkinter as tk
from copy import deepcopy
from datetime import datetime

import custom_exceptions as cexc
import database_interface as dbi
import mcn_status as mcs
from custom_classes import Generic
from functionals import list_subtract
from gui_constants import *

MAX_LENGTH = 64
D_HEIGHT = 8

if __name__ == '__main__':
    # If debugging, just slap together an emitter that prints to the console
    emitter = Generic()
    emitter.info = print
    emitter.warning = print
else:
    emitter = logging.getLogger()


def fetch_data_from_db(down_select: list = None):
    """ Pulls queue documents from the Queue collection and Historical collection

    Collection keys are replaced to the queue name (were _id parameter)

    :param down_select: a list of queue names to be included
    :return: (queue collection queue documents: dict, historical collection queue documents: dict)
    """
    doc = dict()
    try:
        doc, _, _ = dbi.query_collection('MC', dbi.COL_QUEUE)
    except cexc.DatabaseRequestError as dre:
        emitter.warning(dre)
    queues = {v.get(dbi.Q_NAME, "?"): v for _, v in doc.items()}

    doc = dict()
    try:
        doc, _, _ = dbi.query_collection('MC', dbi.COL_HISTORICAL)
    except cexc.DatabaseRequestError as dre:
        emitter.warning(dre)
    historical = {v.get(dbi.Q_NAME, "?"): v
                  for _, v in doc.items()
                  if "?" not in v.get(dbi.Q_NAME, "?") + v.get(dbi.Q_STATUS, "?")}

    if down_select:
        for q_name in list_subtract(list(queues.keys()), down_select):
            del queues[q_name]
        for q_name in list_subtract(list(historical.keys()), down_select):
            del historical[q_name]

    return queues, historical


class QueueRepository:
    """ Container for queue documents, such that every single update does not require pulling the entire database over
    and over again
    """
    def __init__(self):
        """ Creates a container for queue documents """
        self.active = dict()
        self.historical = dict()
        self.update()

    def update(self):
        """ Updates the container

        :return: None
        """
        # print("update")
        try:
            doc, _, _ = dbi.query_collection('MC', dbi.COL_QUEUE)
        except cexc.DatabaseRequestError as dre:
            emitter.warning(dre)
        else:
            self.active = {v.get(dbi.Q_NAME, "?"): v for _, v in doc.items()}

        try:
            doc, _, _ = dbi.query_collection('MC', dbi.COL_HISTORICAL)
        except cexc.DatabaseRequestError as dre:
            emitter.warning(dre)
        else:
            self.historical = {v.get(dbi.Q_NAME, "?"): v
                               for _, v in doc.items()
                               if "?" not in v.get(dbi.Q_NAME, "?") + v.get(dbi.Q_STATUS, "?")}

    def get_queue_names(self):
        """ retrieves the names of all queues stored

        :return: a list of queue names
        """
        active = [k for k in self.active if self.active[k].get(dbi.Q_RECORD, [])]
        hist = [k for k in self.historical if self.historical[k].get(dbi.Q_RECORD, [])]
        return active + hist

    def get_relevant_instruments(self, sel_queue=None):
        """ Searches a queue/collection for all instruments mentioned in a Fault

        :param sel_queue: Reduces the search from both collections to a single queue
        :return: a list of instruments name
        """
        instrument_list = []
        for q_name in self.active:
            if sel_queue and (q_name != sel_queue):
                continue
            for fault_dict in self.active[q_name].get(dbi.Q_RECORD, []):
                loc = mcs.Fault.loadd(fault_dict).location
                if loc not in instrument_list:
                    instrument_list.append(loc)
        for q_name in self.historical:
            if sel_queue and (q_name != sel_queue):
                continue
            for fault_dict in self.historical[q_name].get(dbi.Q_RECORD, []):
                loc = mcs.Fault.loadd(fault_dict).location
                if loc not in instrument_list:
                    instrument_list.append(loc)
        return instrument_list

    def get_all_queues(self):
        """ retrieves all queues as a single document

        :return: A dictionary of queues keyed by name
        """
        return {**self.active, **self.historical}

    def pull_faults(self, q_name):
        """ Retrieves all Faults associated with a given queue

        :param q_name: a queue name
        :return: a list of Faults
        """
        if q_name in self.active:
            collection = self.active
        elif q_name in self.historical:
            collection = self.historical
        else:
            return
        return [mcs.Fault.loadd({**f, 'queue': q_name}) for f in collection[q_name].get(dbi.Q_RECORD, [])]

    def determine_collection(self, q_name):
        """ Determines the collection a queue name was pulled from

        :param q_name: a queue name
        :return: the name of the collection or None
        """
        if q_name in self.active:
            return dbi.COL_QUEUE
        elif q_name in self.historical:
            return dbi.COL_HISTORICAL
        else:
            return

    def _pre(self, sel_queue, sel_fault: mcs.Fault):
        """ Handle pre-edit requirements, such as moving historical documents into an editable mode and getting both the
        original document and making a local copy.

        :param sel_queue: The name of the queue being accessed
        :param sel_fault: The Fault being modified
        :return: (The original document, a new document for editing, the Fault) or None
        """
        collection = self.determine_collection(sel_queue)
        if collection == 'historical':
            try:
                temp, _, _ = dbi.query_document('MC', dbi.COL_HISTORICAL, dbi.Q_NAME, sel_queue)
                temp = temp[dbi.DBG_ID]
                dbi.move_database_document('MC', temp, dbi.COL_HISTORICAL, dbi.COL_QUEUE)
            except cexc.DatabaseRequestError as dre:
                emitter.warning(f"Failed to prepare a historical document for editing: {repr(dre)}")
                return
        elif collection is None:
            emitter.warning(f"Failed to determine collection for {sel_queue}")
            return

        try:
            original_doc, _, _ = dbi.query_document('MC', dbi.COL_QUEUE, dbi.Q_NAME, sel_queue)
        except cexc.DatabaseRequestError as dre:
            emitter.warning(f"Could not pull document from database: {repr(dre)}")
            return

        new_doc = deepcopy(original_doc)
        existing_faults = [mcs.Fault.loadd({**fault_dict, 'queue': sel_queue}) for fault_dict in new_doc[dbi.Q_RECORD]]
        if sel_fault not in existing_faults:
            emitter.warning("Could not find selected Fault in queue document")
            return

        return original_doc, new_doc, existing_faults.index(sel_fault)

    def _post(self, sel_queue, original_doc, new_doc):
        """ Handles post-editing details, such as moving historical queues back.

        :param sel_queue: The name of the queue being edits
        :param original_doc: The original document
        :param new_doc: The edited document
        :return: None if failed, True if successful
        """
        collection = self.determine_collection(sel_queue)
        if collection is None:
            return
        try:
            dbi.update_document('MC', dbi.COL_QUEUE, original_doc, new_doc)
        except cexc.DatabaseRequestError as dre:
            emitter.warning(f"Could not push to database: {repr(dre)}")
            return

        if collection == 'historical':
            try:
                temp, _, _ = dbi.query_document('MC', dbi.COL_QUEUE, dbi.Q_NAME, sel_queue)
                temp = temp[dbi.DBG_ID]
                dbi.move_database_document('MC', temp, dbi.COL_QUEUE, dbi.COL_HISTORICAL)
            except cexc.DatabaseRequestError as dre:
                emitter.warning(f"Failed to restore document to historical: {repr(dre)}")
                return
        return True

    def destroy_fault(self, sel_queue, sel_fault):
        """ Removes a Fault from a queue document

        :param sel_queue: The name of the queue
        :param sel_fault: The Fault being removed
        :return: None or True if successful
        """
        pre = self._pre(sel_queue, sel_fault)
        if not pre:
            return
        original_doc, new_doc, fault_index = pre

        fault_record = new_doc[dbi.Q_RECORD]
        fault_record.pop(fault_index)

        return self._post(sel_queue, original_doc, new_doc)

    def edit_fault(self, sel_queue, sel_fault, sel_loc):
        """ Reassigns the location (blame) of a a Fault

        :param sel_queue: The name of the queue being edited
        :param sel_fault: The Fault being modified
        :param sel_loc: The new Location associated with the Fault
        :return: None or True if successful
        """
        pre = self._pre(sel_queue, sel_fault)
        if not pre:
            return
        original_doc, new_doc, fault_index = pre

        fault_record = new_doc[dbi.Q_RECORD]
        origin = fault_record[fault_index]['location']
        fault_record[fault_index]['location'] = sel_loc
        stamp = f" (Edit: {origin} --> {sel_loc}; {datetime.now().strftime(mcs.TIME_FORMAT)})"
        if isinstance(fault_record[fault_index].get('data', None), str):
            fault_record[fault_index]['data'] += stamp
        else:
            fault_record[fault_index]['data'] = str(fault_record[fault_index].get('data', None)) + stamp

        post = self._post(sel_queue, original_doc, new_doc)
        return post


class Arbitrator:
    """ The GUI for modifying Faults """
    def __init__(self, master):
        """ Creates a GUI for deleting or reattributing blame for Faults

        :param master: a tkinter root
        """
        self.core = master
        self.queue_repo = QueueRepository()

        self.all_systems = mcs.MCN_CFG[mcs.S_ALL_SYS]
        self.constituents = ['queue', 'location', 'level', 'timestamp', 'data']
        self.variables = dict()
        self.forms = dict()
        self.options_menus = dict()
        self.update_queue_names()
        self.fault_list = None
        self.fault_lookup = dict()

        self.construct_selection_pane(tk.Frame(self.core))
        self.construct_edit_pane(tk.Frame(self.core))

    def update_queue_names(self):
        """ Updates the list of queue names in the selector

        :return: ['Select Queue', ...queue names]
        """
        self.variables['q_list'] = ['Select Queue', ] + self.queue_repo.get_queue_names()

    def get_relevant_instruments(self, q_name=None):
        """ Provides a list of afflicted instruments for a collection/queue

        :param q_name: If None, pulls for a collection; if a queue name, pulls for the specific queue
        :return: ['--All', ...instrument names]
        """
        return ["--All", ] + self.queue_repo.get_relevant_instruments(q_name)

    def construct_selection_pane(self, selection_pane: tk.Frame):
        """ Builds the pane where fauls are listed and selected

        :param selection_pane: a tkinter Frame
        :return: None
        """
        selection_pane.pack(**LB33, expand=True)

        # Create sub-panes
        filter_pane = tk.LabelFrame(selection_pane, text="Filter")
        filter_pane.pack(**TBA3)

        selector_pane = tk.LabelFrame(selection_pane, text="Select Fault")
        selector_pane.pack(**TB33, expand=True)

        # populate filter_pane
        tk.Label(filter_pane, text="Queue").grid(row=0, column=0, **GRID_B33)
        self.variables['filter_queue'] = tk.StringVar(self.core, value="Select Queue")
        self.options_menus['filter_queue'] = tk.OptionMenu(filter_pane, self.variables['filter_queue'],
                                                           *self.variables['q_list'], command=self.update_rel_inst)
        self.options_menus['filter_queue'].grid(row=0, column=1, columnspan=2, **GRID_B33)

        tk.Label(filter_pane, text="Instrument").grid(row=1, column=0, **GRID_B33)
        self.variables['filter_inst'] = tk.StringVar(self.core, value="--All")
        self.options_menus['filter_inst'] = tk.OptionMenu(filter_pane, self.variables['filter_inst'],
                                                          *self.get_relevant_instruments())
        self.options_menus['filter_inst'].grid(row=1, column=1, columnspan=2, **GRID_B33)

        tk.Label(filter_pane, text="Sorting").grid(row=2, column=0, **GRID_B33)
        self.variables['sort_mode'] = tk.StringVar(self.core, 'i')  # i - instrument, d - date
        tk.Radiobutton(filter_pane, text="Instrument", variable=self.variables['sort_mode'], value='i')\
            .grid(row=2, column=1, **GRID_B33)
        tk.Radiobutton(filter_pane, text="Date", variable=self.variables['sort_mode'], value='d')\
            .grid(row=2, column=2, **GRID_B33)

        tk.Button(filter_pane, text="Search", command=self.build_fault_list)\
            .grid(row=3, column=0, columnspan=3, **GRID_B33)

        # sketch selector_pane
        _pane = tk.Frame(selector_pane)
        scrolly = tk.Scrollbar(_pane)
        scrolly.pack(side=tk.RIGHT, fill=tk.Y)
        scrollx = tk.Scrollbar(_pane, orient=tk.HORIZONTAL)
        scrollx.pack(side=tk.BOTTOM, fill=tk.X)
        self.fault_list = tk.Listbox(_pane,
                                     yscrollcommand=scrolly.set,
                                     xscrollcommand=scrollx.set,
                                     width=MAX_LENGTH)
        self.fault_list.pack(fill=tk.BOTH, expand=True)
        scrolly.config(command=self.fault_list.yview)
        scrollx.config(command=self.fault_list.xview)
        self.fault_list.bind('<<ListboxSelect>>', self.update_fault_details)
        _pane.pack(**TB33, expand=True)

    def construct_edit_pane(self, edit_pane: tk.Frame):
        """ Creates the Pane in which Faults are edited

        :param edit_pane: a tkinter Frame
        :return: None
        """
        edit_pane.pack(**LB33, expand=True)

        # Create sub-panes
        display_pane = tk.LabelFrame(edit_pane)
        display_pane.pack(**LBA3, expand=True)

        button_pane = tk.LabelFrame(edit_pane)
        button_pane.pack(**LBA3)

        # populate display_pane
        i_data = 0
        for i, item in enumerate(self.constituents):
            if item == 'data':
                i_data = i
                _h = D_HEIGHT
                _s = 'nesw'
            else:
                _h = 1
                _s = 'ew'
            self.forms[f'{item}_label'] = tk.Label(display_pane, text=f'{item.capitalize()}: ')
            self.forms[f'{item}_label'].grid(row=i, column=0)
            self.forms[f'{item}_info'] = tk.Text(display_pane, height=_h, wrap=tk.WORD)
            self.forms[f'{item}_info'].insert(tk.END, ' ' * MAX_LENGTH)
            self.forms[f'{item}_info'].grid(row=i, column=1, sticky=_s)
            self.forms[f'{item}_info'].delete('1.0', tk.END)
        display_pane.rowconfigure(i_data, weight=1)
        display_pane.columnconfigure(1, weight=1)

        # populate button_pane
        self.variables['new_loc'] = tk.StringVar(button_pane, value="")
        tk.Button(button_pane, text='User Error', command=self.mea_culpa).grid(row=0, column=0, pady=12)
        tk.OptionMenu(button_pane, self.variables['new_loc'], *self.all_systems).grid(row=1, column=0)
        tk.Button(button_pane, text='Redirect Blame', command=self.reattribute).grid(row=2, column=0)
        tk.Button(button_pane, text='Erase Fault', command=self.destroy_fault).grid(row=3, column=0, pady=12)

    def update_rel_inst(self, *_):
        """ Updates the options for instruments when a queue is selected

        :return: None
        """
        if 'filter_inst' not in self.options_menus:
            return
        menu = self.options_menus['filter_inst']
        sel_queue = self.variables['filter_queue'].get()
        if sel_queue == "Select Queue":
            return
        inst_list = self.get_relevant_instruments(sel_queue)
        menu['menu'].delete(0, tk.END)
        for i in inst_list:
            menu['menu'].add_command(label=i, command=lambda v=i: self.variables['filter_inst'].set(v))

    def build_fault_list(self):
        """ Constructs a list of Faults when a queue is selected

        :return: None
        """
        self.clear_form()

        sel_queue = self.variables['filter_queue'].get()
        if sel_queue == "Select Queue":
            emitter.info("Select a queue first")
            return

        all_faults = self.queue_repo.pull_faults(sel_queue)
        self.fault_lookup = {f'{f.location} [{f.timestamp}]: {f.data[:64] if f.data else "..."}': f for f in all_faults}

        sel_location = self.variables['filter_inst'].get()
        if sel_location != '--All':
            self.fault_lookup = {k: v for k, v in self.fault_lookup.items() if v.location == sel_location}
        all_faults = [k for k in self.fault_lookup]

        i_sort_kwargs = {'key': lambda x: self.fault_lookup[x].location}
        d_sort_kwargs = {'key': lambda x: datetime.strptime(self.fault_lookup[x].timestamp, mcs.TIME_FORMAT),
                         'reverse': True}
        if self.variables['sort_mode'].get() == 'i':
            all_faults = sorted(all_faults, **d_sort_kwargs)
            all_faults = sorted(all_faults, **i_sort_kwargs)
        else:
            all_faults = sorted(all_faults, **i_sort_kwargs)
            all_faults = sorted(all_faults, **d_sort_kwargs)

        self.fault_list.delete(0, tk.END)
        for fault in all_faults:
            self.fault_list.insert(tk.END, fault)

    def clear_form(self):
        """ Clears the Fault report/edit frame when the Fault changes

        :return: None
        """
        if not self.forms:
            return
        for i, item in enumerate(self.constituents):
            self.forms[f'{item}_info'].delete('1.0', tk.END)
            self.forms[f'{item}_info'].insert(tk.END, ' ' * MAX_LENGTH)
            self.forms[f'{item}_info'].grid(row=i, column=1)
            self.forms[f'{item}_info'].delete('1.0', tk.END)

    def get_sel_fault(self):
        """ Retrieves a Fault from the list based on user selection

        :return: The Fault or None if nothing is selected or it cannot be retrieved
        """
        index = self.fault_list.curselection()[:1]
        if not index:
            return self.clear_form()
        sel_fault_key = self.fault_list.get(index)
        sel_fault = self.fault_lookup.get(sel_fault_key, None)
        return sel_fault

    def update_fault_details(self, *_):
        """ Updates the Fault viewing/editing pane with the details of a selected Fault

        :return: None
        """
        q_name, fault, _ = self.get_update_info()
        if None in [q_name, fault]:
            return self.clear_form()
        loc = fault.location
        if loc:
            self.variables['new_loc'].set(loc)

        for i, item in enumerate(self.constituents):
            self.forms[f'{item}_info'].delete('1.0', tk.END)
            information = getattr(fault, item)
            information = information if information else "NR"
            if (item == 'queue') and (information == "NR"):
                print("DEBUG 2")
                information = q_name
            # _m = D_HEIGHT if item == 'data' else 1
            # working_length = MAX_LENGTH * _m
            # if len(information) > working_length:
            #     information = information[:(working_length - 4)] + " ..."
            self.forms[f'{item}_info'].insert(tk.END, information)

    def mea_culpa(self):
        """ Shorthand for calling self.reattribute() with the new location being the User (_U)

        :return:
        """
        self.variables['new_loc'].set('_U')
        self.reattribute()

    def reattribute(self):
        """ Reassigns the location (blame) for a Fault

        :return:
        """
        q_name, fault, location = self.get_update_info()
        if None in [q_name, fault, location]:
            return
        # print('reattribute', q_name, location)
        self.queue_repo.edit_fault(q_name, fault, location)
        self.queue_repo.update()
        self.build_fault_list()

    def destroy_fault(self):
        """ Deletes a Fault from its associated Queue document

        :return:
        """
        q_name, fault, _ = self.get_update_info()
        if None in [q_name, fault]:
            return
        # print('reattribute', q_name)
        self.queue_repo.destroy_fault(q_name, fault)
        self.queue_repo.update()
        self.build_fault_list()

    def get_update_info(self):
        """ Collects information for a Fault update

        :return: (the queue name, the Fault, the new location)
        """
        sel_queue = self.variables['filter_queue'].get()
        sel_queue = sel_queue if sel_queue != "Select Queue" else None
        sel_fault = self.get_sel_fault()
        location = self.variables['new_loc'].get()
        location = location if location else None
        return sel_queue, sel_fault, location


if __name__ == '__main__':
    core = tk.Tk()
    all_queues = False

    if all_queues:
        sel_queues = None
    else:
        sel_queues = [
            "second_sulfur_explore_scaffold_round2_characterization_20220804",
            "sulfur_explore_scaffold_round2_characterization_20220804"
        ]

    my_manager = Arbitrator(core)

    core.mainloop()
