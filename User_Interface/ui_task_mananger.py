""" An interface for controlling the execution of tasks on the MCN
@author: Ben C
"""

import tkinter
from _tkinter import TclError

import custom_exceptions as cexc
import queue_builder
import functionals
import database_interface as dbi
from constants import *
from database_constants import *
from gui_constants import *
from mcn_logging_manager import system_log
import mcn_status as mcs
from thread_safe_data import ThreadSafeDataContainer as TSD  # noqa


class TaskManager:
    """ A GUI for controlling database-based operation (MC Only) """
    def __init__(self, master, data_pipes: dict, safety):
        """ Creates a GUI for controlling how the MC should behave and execute tasks as well as some helper methods
        to interface with the database

        :param master: a tkinter root
        :param data_pipes: for communication with the MC (parent System)
        :param safety: Used to ensure only one copy is open at a time
        """
        self.safety = safety
        self.safety['state'] = tkinter.DISABLED
        self.core = master
        self.core.protocol("WM_DELETE_WINDOW", self.on_closing)
        self._pipes = data_pipes
        self.child_com: TSD = data_pipes[CHILD_COM]
        internals: mcs.Status = data_pipes[INTERNAL]
        self.name = internals.name
        with self.child_com as cc:
            self.expired = cc['ui'][mcs.C_EXPIRED]

        self.master_frame = tkinter.Frame(self.core)
        self.master_frame.winfo_toplevel().title("MCN Task Manager")
        self.master_frame.pack(expand=True, **BOTH33)

        self.tk = dict()

        # Memory
        self.rx_mode = tkinter.BooleanVar(self.master_frame, False)  # True-Safe mode, False-Automatic mode
        self.rx_queue_collection = dict()  # Queue collection with names as keys
        self.rx_selected_queues = list()
        self.rx_blocked_queues = list()
        self.rx_sel_q = str()  # Queue selected for editing (name)
        self.rx_sel_s = str()  # Queue step selected for editing (repr)

        # Build GUI scaffold
        self.sketch()
        self.packer()

        self.update_queue_collection()

    def on_closing(self, *_):
        """ Restores the safety parameter upon exit

        :return:
        """
        self.core.destroy()
        self.safety['state'] = tkinter.NORMAL

    def sketch(self):
        """ Creates all tkinter objects

        :return:
        """
        self.tk = dict()

        # Structure
        self.tk['manager'] = dict()
        self.tk['manager']['control'] = dict()
        self.tk['manager']['settings'] = dict()
        self.tk['manager']['settings']['tors'] = dict()
        self.tk['manager']['settings']['tions'] = dict()
        self.tk['manager']['mode'] = dict()

        self.tk['editor'] = dict()
        self.tk['editor']['qselect'] = dict()
        self.tk['editor']['qselect']['selector'] = dict()
        self.tk['editor']['qselect']['mover'] = dict()
        self.tk['editor']['sselect'] = dict()
        self.tk['editor']['sselect']['selector'] = dict()
        self.tk['editor']['sselect']['marker'] = dict()

        # Level 1
        self.tk['manager.frame'] = tkinter.Frame(self.master_frame)
        self.tk['editor.frame'] = tkinter.Frame(self.master_frame)

        # Level 2 - manager
        self.tk['manager']['control.frame'] = tkinter.LabelFrame(self.tk['manager.frame'], text="Control")
        self.tk['manager']['settings.frame'] = tkinter.LabelFrame(self.tk['manager.frame'], text="Settings")
        self.tk['manager']['mode.frame'] = tkinter.LabelFrame(self.tk['manager.frame'], text="Modus Operandi")
        # Level 2 - editor
        self.tk['editor']['qselect.frame'] = tkinter.Frame(self.tk['editor.frame'])
        self.tk['editor']['sselect.frame'] = tkinter.Frame(self.tk['editor.frame'])

        # Level 3 - control
        self.tk['manager']['control']['db control'] = tkinter.Button(self.tk['manager']['control.frame'], text="Start Database Control", command=self.change_mode, width=18)
        self.tk['manager']['control']['db pull'] = tkinter.Button(self.tk['manager']['control.frame'], text="Pull Database", command=self.update_queue_collection, width=18)
        # Level 3 - settings
        self.tk['manager']['settings']['select'] = tkinter.Label(self.tk['manager']['settings.frame'], text="Selected")
        self.tk['manager']['settings']['block'] = tkinter.Label(self.tk['manager']['settings.frame'], text="Blocked")
        self.tk['manager']['settings']['tors.frame'] = tkinter.Frame(self.tk['manager']['settings.frame'])
        self.tk['manager']['settings']['tions.frame1'] = tkinter.Frame(self.tk['manager']['settings.frame'])
        self.tk['manager']['settings']['tions.frame2'] = tkinter.Frame(self.tk['manager']['settings.frame'])
        self.tk['manager']['settings']['allow all'] = tkinter.Button(self.tk['manager']['settings.frame'], text="Allow All", command=self.allow_all)
        self.tk['manager']['settings']['block all'] = tkinter.Button(self.tk['manager']['settings.frame'], text="Block All", command=self.block_all)
        # Level 3 - mode
        self.tk['manager']['mode']['auto'] = tkinter.Radiobutton(self.tk['manager']['mode.frame'], text="Automatic", variable=self.rx_mode, value=False, command=self.toggle_safe_mode)
        self.tk['manager']['mode']['safe'] = tkinter.Radiobutton(self.tk['manager']['mode.frame'], text="Safe", variable=self.rx_mode, value=True, command=self.toggle_safe_mode)
        # Level 3 - qselect
        self.tk['editor']['qselect']['selector.frame'] = tkinter.LabelFrame(self.tk['editor']['qselect.frame'], text="Queue Select")
        self.tk['editor']['qselect']['mover.frame'] = tkinter.LabelFrame(self.tk['editor']['qselect.frame'], text="Move")
        # Level 3 - sselect
        self.tk['editor']['sselect']['selector.frame'] = tkinter.LabelFrame(self.tk['editor']['sselect.frame'], text="Step Select")
        self.tk['editor']['sselect']['marker.frame'] = tkinter.LabelFrame(self.tk['editor']['sselect.frame'], text="Edit")

        # Level 4 - tors
        self.tk['manager']['settings']['tors']['selector'] = tkinter.Button(self.tk['manager']['settings']['tors.frame'], text="<", command=lambda x="sel": self.tor(x))
        self.tk['manager']['settings']['tors']['deselector'] = tkinter.Button(self.tk['manager']['settings']['tors.frame'], text=">", command=lambda x="block": self.tor(x))
        # Level 4 - tions
        self.tk['manager']['settings']['tions']['scrolly1'] = tkinter.Scrollbar(self.tk['manager']['settings']['tions.frame1'])
        self.tk['manager']['settings']['tions']['scrolly2'] = tkinter.Scrollbar(self.tk['manager']['settings']['tions.frame2'])
        self.tk['manager']['settings']['tions']['scrollx1'] = tkinter.Scrollbar(self.tk['manager']['settings']['tions.frame1'], orient=tkinter.HORIZONTAL)
        self.tk['manager']['settings']['tions']['scrollx2'] = tkinter.Scrollbar(self.tk['manager']['settings']['tions.frame2'], orient=tkinter.HORIZONTAL)
        self.tk['manager']['settings']['tions']['selected'] = tkinter.Listbox(self.tk['manager']['settings']['tions.frame1'], selectmode=tkinter.MULTIPLE, yscrollcommand=self.tk['manager']['settings']['tions']['scrolly1'].set, xscrollcommand=self.tk['manager']['settings']['tions']['scrollx1'].set, width=56)
        self.tk['manager']['settings']['tions']['blocked'] = tkinter.Listbox(self.tk['manager']['settings']['tions.frame2'], selectmode=tkinter.MULTIPLE, yscrollcommand=self.tk['manager']['settings']['tions']['scrolly2'].set, xscrollcommand=self.tk['manager']['settings']['tions']['scrollx2'].set, width=56)
        self.tk['manager']['settings']['tions']['scrolly1'].config(command=self.tk['manager']['settings']['tions']['selected'].yview)
        self.tk['manager']['settings']['tions']['scrolly2'].config(command=self.tk['manager']['settings']['tions']['blocked'].yview)
        self.tk['manager']['settings']['tions']['scrollx1'].config(command=self.tk['manager']['settings']['tions']['selected'].xview)
        self.tk['manager']['settings']['tions']['scrollx2'].config(command=self.tk['manager']['settings']['tions']['blocked'].xview)
        # Level 4 - Q-selector
        self.tk['editor']['qselect']['selector.scrolly1'] = tkinter.Scrollbar(self.tk['editor']['qselect']['selector.frame'])
        self.tk['editor']['qselect']['selector.scrollx1'] = tkinter.Scrollbar(self.tk['editor']['qselect']['selector.frame'], orient=tkinter.HORIZONTAL)
        self.tk['editor']['qselect']['selector']['qselect'] = tkinter.Listbox(self.tk['editor']['qselect']['selector.frame'],
                                                                              selectmode=tkinter.SINGLE,
                                                                              yscrollcommand=self.tk['editor']['qselect']['selector.scrolly1'].set,
                                                                              xscrollcommand=self.tk['editor']['qselect']['selector.scrollx1'].set,
                                                                              width=56)
        self.tk['editor']['qselect']['selector.scrolly1'].config(command=self.tk['editor']['qselect']['selector']['qselect'].yview)
        self.tk['editor']['qselect']['selector.scrollx1'].config(command=self.tk['editor']['qselect']['selector']['qselect'].xview)
        self.tk['editor']['qselect']['selector']['qselect'].bind('<<ListboxSelect>>', self.select_queue_for_editing)
        # Level 4 - mover
        self.tk['editor']['qselect']['mover']['move2idle'] = tkinter.Button(self.tk['editor']['qselect']['mover.frame'], text="Move to Idle", command=lambda x="idle": self.move_group(x), width=18)
        self.tk['editor']['qselect']['mover']['move2working'] = tkinter.Button(self.tk['editor']['qselect']['mover.frame'], text="Move to Working", command=lambda x="in_progress": self.move_group(x), width=18)
        self.tk['editor']['qselect']['mover']['move2error'] = tkinter.Button(self.tk['editor']['qselect']['mover.frame'], text="Move to Error", command=lambda x="error": self.move_group(x), width=18)
        # Level 4 - S-selector
        self.tk['editor']['sselect']['selector.scrolly1'] = tkinter.Scrollbar(self.tk['editor']['sselect']['selector.frame'])
        self.tk['editor']['sselect']['selector']['sselect'] = tkinter.Listbox(self.tk['editor']['sselect']['selector.frame'], selectmode=tkinter.SINGLE, yscrollcommand=self.tk['editor']['sselect']['selector.scrolly1'].set, width=48)
        self.tk['editor']['sselect']['selector.scrolly1'].config(command=self.tk['editor']['sselect']['selector']['sselect'].yview)
        self.tk['editor']['sselect']['selector']['sselect'].bind('<<ListboxSelect>>', self.select_step_for_editing)
        # Level 4 - marker
        self.tk['editor']['sselect']['marker']['set2complete'] = tkinter.Button(self.tk['editor']['sselect']['marker.frame'], text="Mark Complete", command=lambda x="yes": self.mark_step(x), width=18)
        self.tk['editor']['sselect']['marker']['set2incomplete'] = tkinter.Button(self.tk['editor']['sselect']['marker.frame'], text="Mark Incomplete", command=lambda x="no": self.mark_step(x), width=18)
        self.tk['editor']['sselect']['marker']['soft_reset_queue'] = tkinter.Button(self.tk['editor']['sselect']['marker.frame'], text="Soft Reset", command=self.soft_reset_queue, width=18)

    def packer(self):
        """ Packs/Grids all the widgets created by sketch()

        :return:
        """
        # Level 1
        self.tk['manager.frame'].pack(expand=True, side=tkinter.LEFT, **BOTH33)
        self.tk['editor.frame'].pack(expand=True, side=tkinter.LEFT, **BOTH33)

        # Level 2 - manager
        self.tk['manager']['control.frame'].pack(side=tkinter.TOP, **BOTH33)
        self.tk['manager']['settings.frame'].pack(expand=True, side=tkinter.TOP, **BOTH33)
        self.tk['manager']['mode.frame'].pack(side=tkinter.TOP, **BOTH33)
        # Level 2 - editor
        self.tk['editor']['qselect.frame'].pack(expand=True, side=tkinter.TOP, **BOTH33)
        self.tk['editor']['sselect.frame'].pack(expand=True, side=tkinter.TOP, **BOTH33)

        # Level 3 - control
        self.tk['manager']['control']['db control'].pack(side=tkinter.LEFT, **BOTH33)
        self.tk['manager']['control']['db pull'].pack(side=tkinter.LEFT, **BOTH33)
        # Level 3 - settings
        self.tk['manager']['settings.frame'].rowconfigure(1, weight=1)
        self.tk['manager']['settings.frame'].columnconfigure((0, 2), weight=1)
        self.tk['manager']['settings']['select'].grid(row=0, column=0)
        self.tk['manager']['settings']['block'].grid(row=0, column=2)
        self.tk['manager']['settings']['tions.frame1'].grid(row=1, column=0, sticky='nsew')
        self.tk['manager']['settings']['tors.frame'].grid(row=1, column=1,)
        self.tk['manager']['settings']['tions.frame2'].grid(row=1, column=2, sticky='nsew')
        self.tk['manager']['settings']['allow all'].grid(row=2, column=0)
        self.tk['manager']['settings']['block all'].grid(row=2, column=2)
        # Level 3 - mode
        self.tk['manager']['mode']['auto'].pack(side=tkinter.LEFT, **BOTH33)
        self.tk['manager']['mode']['safe'].pack(side=tkinter.LEFT, **BOTH33)
        # Level 3 - qselect
        self.tk['editor']['qselect']['selector.frame'].pack(expand=True, side=tkinter.LEFT, **BOTH33)
        self.tk['editor']['qselect']['mover.frame'].pack(side=tkinter.LEFT, **BOTH33)
        self.tk['editor']['qselect']['selector.scrolly1'].pack(side=tkinter.RIGHT, fill=tkinter.Y)
        # level 3 - sselect
        self.tk['editor']['sselect']['selector.frame'].pack(expand=True, side=tkinter.LEFT, **BOTH33)
        self.tk['editor']['sselect']['marker.frame'].pack(side=tkinter.LEFT, **BOTH33)
        self.tk['editor']['sselect']['selector.scrolly1'].pack(side=tkinter.RIGHT, fill=tkinter.Y)
        self.tk['editor']['sselect']['selector']['sselect'].pack(expand=True, side=tkinter.LEFT, **BOTH33)
        self.tk['editor']['sselect']['marker']['set2complete'].pack(side=tkinter.TOP, **BOTH33)
        self.tk['editor']['sselect']['marker']['set2incomplete'].pack(side=tkinter.TOP, **BOTH33)
        self.tk['editor']['sselect']['marker']['soft_reset_queue'].pack(side=tkinter.TOP, **BOTH33)

        # Level 4 - tors
        self.tk['manager']['settings']['tors']['selector'].pack(side=tkinter.LEFT)
        self.tk['manager']['settings']['tors']['deselector'].pack(side=tkinter.LEFT)
        # Level 4 - tions
        self.tk['manager']['settings']['tions']['scrolly1'].pack(side=tkinter.RIGHT, fill=tkinter.Y)
        self.tk['manager']['settings']['tions']['scrolly2'].pack(side=tkinter.RIGHT, fill=tkinter.Y)
        self.tk['manager']['settings']['tions']['scrollx1'].pack(side=tkinter.BOTTOM, fill=tkinter.X)
        self.tk['manager']['settings']['tions']['scrollx2'].pack(side=tkinter.BOTTOM, fill=tkinter.X)
        self.tk['manager']['settings']['tions']['selected'].pack(expand=True, side=tkinter.LEFT, **BOTH33)
        self.tk['manager']['settings']['tions']['blocked'].pack(expand=True, side=tkinter.LEFT, **BOTH33)
        # Level 4 - Q-selector
        self.tk['editor']['qselect']['selector.scrolly1'].pack(side=tkinter.RIGHT, fill=tkinter.Y)
        self.tk['editor']['qselect']['selector.scrollx1'].pack(side=tkinter.BOTTOM, fill=tkinter.X)
        self.tk['editor']['qselect']['selector']['qselect'].pack(expand=True, side=tkinter.LEFT, **BOTH33)
        # Level 4 - mover
        self.tk['editor']['qselect']['mover']['move2idle'].pack(side=tkinter.TOP, **BOTH33)
        self.tk['editor']['qselect']['mover']['move2working'].pack(side=tkinter.TOP, **BOTH33)
        self.tk['editor']['qselect']['mover']['move2error'].pack(side=tkinter.TOP, **BOTH33)

    def change_mode(self, *_):
        """ Changes the mode (DB control on/off) and enables/disables GUI elements accordingly

        :return:
        """
        with self.child_com as cc:
            cc['ui'][mcs.C_ACTIVE_MODE] = not cc['ui'][mcs.C_ACTIVE_MODE]
            current_mode = cc['ui'][mcs.C_ACTIVE_MODE]
        if current_mode:
            system_log.mark_up()  # .info("$TASK_MANAGER$ : DB CONTROL *ON*")
            self.tk['manager']['control']['db control']['text'] = "Pause DB Control"
            self.recursive_status_change(self.tk['manager']['settings'], 'state', tkinter.DISABLED)
            self.recursive_status_change(self.tk['manager']['mode'], 'state', tkinter.DISABLED)
            self.recursive_status_change(self.tk['editor'], 'state', tkinter.DISABLED)
        else:
            system_log.mark_down()  # .info("$TASK_MANAGER$ : DB CONTROL *OFF*")
            self.tk['manager']['control']['db control']['text'] = "Resume DB Control"
            self.recursive_status_change(self.tk['manager']['settings'], 'state', tkinter.NORMAL)
            self.recursive_status_change(self.tk['manager']['mode'], 'state', tkinter.NORMAL)
            self.recursive_status_change(self.tk['editor'], 'state', tkinter.NORMAL)
            self.rx_sel_q = str()
            self.rx_sel_s = str()
        self.update_queue_collection()

    def recursive_status_change(self, domain, field, state):
        """ Goes through a domain and disables/enables all widgets therein

        :param domain: A hierarchical dictionary of tkinter widgets
        :param field: the key (widget[field] = state)
        :param state: the value  (widget[field] = state)
        :return:
        """
        if isinstance(domain, dict):
            for _, item in domain.items():
                self.recursive_status_change(item, field, state)
        else:
            try:
                domain[field] = state
            except KeyError:
                pass
            except TclError:
                pass

    def update_queue_collection(self, *_):
        """ Pulls DB for an update

        :return:
        """
        try:
            doc, err_details, resp_code = dbi.query_collection("MC", "queue")
        except cexc.DatabaseRequestError:
            system_log.exception("Task Manager - DB request failed - Queue lookup")
            return
        self.rx_queue_collection = {v[Q_NAME]: v for k, v in doc.items()}
        self.update_queue_lists()

    def allow_all(self, *_):
        """ Removes all queues from the block list

        :return:
        """
        with self.child_com as cc:
            cc['ui'][mcs.C_QUEUE_BLOCK] = list()
        self.update_queue_collection()

    def block_all(self, *_):
        """ Blocks all queues

        :return:
        """
        with self.child_com as cc:
            cc['ui'][mcs.C_QUEUE_BLOCK] = list(self.rx_queue_collection.keys())
        self.update_queue_collection()

    def toggle_safe_mode(self, *_):
        """ Changes the state of the safe-mode parameter (makes Scheduler ask the user to confirm
        before a command is sent)

        :return:
        """
        mode = self.rx_mode.get()
        with self.child_com as cc:
            cc['ui'][mcs.C_SAFE_MODE] = mode

    def tor(self, x, *_):
        """ Command that blocks or unblocks a queue (or set of queues)

        :param x: "sel" or unblocking or "block" for blocking
        :return:
        """
        if x == "sel":
            selection = list()
            cursor = self.tk['manager']['settings']['tions']['blocked'].curselection()
            if not cursor:
                return
            for q_index in cursor:
                selection.append(self.tk['manager']['settings']['tions']['blocked'].get(q_index))
            with self.child_com as cc:
                cc['ui'][mcs.C_QUEUE_BLOCK] = functionals.list_subtract(cc['ui'][mcs.C_QUEUE_BLOCK], selection)
        elif x == "block":
            selection = list()
            cursor = self.tk['manager']['settings']['tions']['selected'].curselection()
            if not cursor:
                return
            for q_index in cursor:
                selection.append(self.tk['manager']['settings']['tions']['selected'].get(q_index))
            with self.child_com as cc:
                cc['ui'][mcs.C_QUEUE_BLOCK] = cc['ui'][mcs.C_QUEUE_BLOCK] + selection
        else:
            pass
        self.update_queue_collection()

    def update_queue_lists(self, *_):
        """ Updates the display of available queues in the manager and editor

        :return:
        """
        self.tk['manager']['settings']['tions']['selected'].delete(0, tkinter.END)
        self.tk['manager']['settings']['tions']['blocked'].delete(0, tkinter.END)
        self.tk['editor']['qselect']['selector']['qselect'].delete(0, tkinter.END)
        self.tk['editor']['sselect']['selector']['sselect'].delete(0, tkinter.END)

        with self.child_com as cc:
            cc_q_block: list = cc['ui'][mcs.C_QUEUE_BLOCK]

        self.rx_blocked_queues = cc_q_block
        self.rx_selected_queues = [k for k in self.rx_queue_collection.keys() if k not in self.rx_blocked_queues]

        try:
            all_qs_long_name = [f"{k} ({v.get(DBQ_STATUS, 'na')}, {dbi.get_step_number(v)})"
                                for k, v in self.rx_queue_collection.items()]
        except cexc.DatabaseRequestError:
            all_qs_long_name = []

        for aq in self.rx_selected_queues:
            self.tk['manager']['settings']['tions']['selected'].insert(tkinter.END, aq)
        for bq in self.rx_blocked_queues:
            self.tk['manager']['settings']['tions']['blocked'].insert(tkinter.END, bq)
        for ql in all_qs_long_name:
            self.tk['editor']['qselect']['selector']['qselect'].insert(tkinter.END, ql)
        # no update for self.tk['editor']['sselect']['selector']['sselect']
        #     it get's updated when a queue is selected for editing

    def select_queue_for_editing(self, *_):
        """ Updates the editing panel when a queue is selected

        :return:
        """
        # Get what queue is selected
        index = self.tk['editor']['qselect']['selector']['qselect'].curselection()[:1]
        if not index:
            return
        q_repr: str = self.tk['editor']['qselect']['selector']['qselect'].get(index)
        new_rx_sel_q, *_ = q_repr.split(" ", 1)

        if (self.rx_sel_q != new_rx_sel_q) or (self.tk['editor']['sselect']['selector']['sselect'].size() == 0):
            # Clear steps
            self.tk['editor']['sselect']['selector']['sselect'].delete(0, tkinter.END)
            # Update the steps
            operations = self.rx_queue_collection[new_rx_sel_q][dbi.Q_OPERATIONS_LIST]
            operation_long_name = [f"{k}: {v[QOP_OPERATION]} ({v[QOP_COMPLETED]})" for k, v in
                                   operations.items()]
            for op in operation_long_name:
                self.tk['editor']['sselect']['selector']['sselect'].insert(tkinter.END, op)
        # Update
        self.rx_sel_q = new_rx_sel_q

    def move_group(self, x, *_):
        """ Changes the state of a queue (idle, in-progress, error)

        :param x: Destination
        :return:
        """
        # Grab the name of the queue
        sel_q_doc = self.rx_queue_collection.get(self.rx_sel_q, None)
        if not sel_q_doc:
            self.update_queue_collection()
            return

        # Change the status
        try:
            dbi.move(self.rx_sel_q, _from="Lookup", _to=x)
        except cexc.DatabaseRequestError:
            pass
        else:
            self.update_queue_collection()

    def select_step_for_editing(self, *_):
        """ Loads a selected step into the register

        :return:
        """
        # Get what step is selected
        index = self.tk['editor']['sselect']['selector']['sselect'].curselection()[:1]
        if not index:
            return
        s_repr: str = self.tk['editor']['sselect']['selector']['sselect'].get(index)
        self.rx_sel_s, *_ = s_repr.split(": ", 1)

    def mark_step(self, x, *_):
        """ Sets the completion of a step

        :param x: yes/no
        :return:
        """
        if not self.rx_sel_q:
            return
        if not self.rx_sel_s:
            return
        try:
            dbi.update_mongo_queue(self.name, self.rx_sel_q, self.rx_sel_s, x)
            if x == 'yes':
                dbi.mark_time(self.rx_sel_q, self.rx_sel_s, stop=True, logger=system_log)
            elif x == 'no':
                dbi.mark_time(self.rx_sel_q, self.rx_sel_s, reset=True, logger=system_log)
        except cexc.DatabaseRequestError:
            pass
        else:
            self.update_queue_collection()

    def soft_reset_queue(self, *_):
        """ invokes reset_queue from queue_builder.py

        :return:
        """
        if not self.rx_sel_q:
            return
        print(queue_builder.reset_queue(self.rx_sel_q))
        self.update_queue_collection()


if __name__ == '__main__':
    from custom_classes import Nop
    root = tkinter.Tk()
    my_state = mcs.Status("MC")
    my_pipes = dict()
    my_pipes[INTERNAL] = my_state
    my_pipes[CHILD_COM] = TSD(
        {'ui': {mcs.C_ACTIVE_MODE: False, mcs.C_QUEUE_BLOCK: list(), mcs.C_EXPIRED: [False, ]},
         'is_initialized': True})

    ui = TaskManager(root, my_pipes, Nop())
    root.mainloop()
