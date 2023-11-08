""" The Tab in the main GUI which manages status elements
@author: Ben C
"""

import mcn_queues
import tkfire as tkf
import mcn_status as mcs
import synchronization as sync
import ui_fault_maker
from mcn_logging_manager import system_log
from functools import wraps
from constants import *
from functionals import list_subtract
import message as message_manager
from ui_checkpoint_doctor import Defibrillator
import database_interface as dbi
from custom_exceptions import DatabaseRequestError
from typing import Tuple, Union

DEBUG = False


# This module makes use of tkfire, a custom tkinter manager.  The only other use is in the queue builder GUI (which
# may not be included in the final release), so additional explanations will be given in this module.  TkFire can
# be found within this project as well as on GitHub at https://github.com/RBCanty/TkFire .

# 'mother' defines the layout of the GUI in dictionary form.  Elements of the GUI are given names (the key) and will
# possess up to three values: (Type) the Widget type (e.g. Frame, Button), (Layout) how the Widget is displayed (e.g.
# pack or grid), and (Children, optional) if the Widget is a Frame or other object which can support having children.
# The Type and Layout fields map to lists which contain [class or method, {kwargs}].
# If the IDE supports collapsing sub-dictionaries, this allows the specification of the GUI to be collapsed and for
# specific elements to be expanded for inspection---hopefully making it easier to read.
# The names used in mother will allow elements to be referenced by name using bang-path notation.  Tkinter does already
# support referencing UI elements by name implicitly, this module makes these references explicit and mandatory.
# Custom Widgets and Tkinter Variables can be included as well, this is done below with 'sel_checkpoint'
# (low_Section!Checkpoints!Selector/type/variable) and with 'checkpoints' (same location/type/values).
# Widgets such as buttons can be "plugged in" (or bound) later by assigning their 'command' field once the parent
# class has been created.  This is done below in the "StatusTab.plug_in()" method.
mother = {
    'top_section': {
        tkf.TYPE: ['Frame', {}],
        tkf.LAYOUT: ['pack', tkf.TB33],
        tkf.CHILDREN: {
            'Network': {
                tkf.TYPE: ['LabelFrame', {'text': "Network"}],
                tkf.LAYOUT: ['pack', tkf.LB33],
                tkf.CHILDREN: {
                    'Connectivity': {
                        tkf.TYPE: ['Label', {'text': "unknown"}],
                        tkf.LAYOUT: ['pack', tkf.LB33]
                    },
                    'IP_and_Port': {
                        tkf.TYPE: ['Label', {'text': "unknown"}],
                        tkf.LAYOUT: ['pack', tkf.LB33]
                    }
                }
            },
            'Members': {
                tkf.TYPE: ['LabelFrame', {'text': "Members"}],
                tkf.LAYOUT: ['pack', tkf.LB33],
                tkf.CHILDREN: {
                    'Area_label': {
                        tkf.TYPE: ['Label', {'text': "Network"}],
                        tkf.LAYOUT: ['grid', tkf.grid_arg(0, 0, tkf.P33)],
                    },
                    'Local_label': {
                        tkf.TYPE: ['Label', {'text': "Local"}],
                        tkf.LAYOUT: ['grid', tkf.grid_arg(1, 0, tkf.P33)],
                    },
                    'Area_members': {
                        tkf.TYPE: ['Label', {'text': "..."}],
                        tkf.LAYOUT: ['grid', tkf.grid_arg(0, 1, tkf.P33)],
                    },
                    'Local_members': {
                        tkf.TYPE: ['Label', {'text': "..."}],
                        tkf.LAYOUT: ['grid', tkf.grid_arg(1, 1, tkf.P33)],
                    }
                }
            }
        }
    },
    'low_section': {
        tkf.TYPE: ['LabelFrame', {'text': "State"}],
        tkf.LAYOUT: ['pack', {'expand': True, **tkf.TB33}],
        tkf.CHILDREN: {
            'Checkpoints': {
                tkf.TYPE: ['Frame', {}],
                tkf.LAYOUT: ['pack', tkf.TB33],
                tkf.CHILDREN: {
                    'Label': {
                        tkf.TYPE: ['Label', {'text': 'Checkpoints'}],
                        tkf.LAYOUT: ['pack', tkf.LB33],
                    },
                    'Selector': {
                        tkf.TYPE: ['OptionMenu',
                                   {'variable': ['sel_checkpoint', {'value': 'None'}], 'values': 'checkpoints'}],
                        tkf.LAYOUT: ['pack', tkf.LB33],
                    },
                    'B_Bypass': {
                        tkf.TYPE: ['Button', {'text': 'Bypass'}],
                        tkf.LAYOUT: ['pack', tkf.LB33],
                    },
                    'B_Retry': {
                        tkf.TYPE: ['Button', {'text': 'Retry'}],
                        tkf.LAYOUT: ['pack', tkf.LB33],
                    },
                    'B_Fail': {
                        tkf.TYPE: ['Button', {'text': 'Fail'}],
                        tkf.LAYOUT: ['pack', tkf.LB33],
                    },
                    'B_Remove': {
                        tkf.TYPE: ['Button', {'text': 'Remove'}],
                        tkf.LAYOUT: ['pack', tkf.LB33],
                    },
                    'B_Revive': {
                        tkf.TYPE: ['Button', {'text': 'Revive'}],
                        tkf.LAYOUT: ['pack', tkf.LB33],
                    }
                }
            },
            'Fault_list_frame': {
                tkf.TYPE: ['Frame', {}],
                tkf.LAYOUT: ['pack', {'expand': True, **tkf.TB33}],
                tkf.CHILDREN: {
                    'Fault_list': {
                        tkf.TYPE: ['Listbox',
                                   {'selectmode': tkf.tk.EXTENDED, 'width': 160, 'height': 12, 'scrolly': True}],
                        tkf.LAYOUT: ['pack', {'expand': True, 'side': 'left', 'fill': 'both'}]
                    }
                }
            },
            'Fault_manager': {
                tkf.TYPE: ['Frame', {}],
                tkf.LAYOUT: ['pack', tkf.TB33],
                tkf.CHILDREN: {
                    'B_Remove': {
                        tkf.TYPE: ['Button', {'text': 'Remove', 'width': 18}],
                        tkf.LAYOUT: ['pack', tkf.L33]
                    },
                    'B_Report': {
                        tkf.TYPE: ['Button', {'text': 'Report Fault', 'width': 18}],
                        tkf.LAYOUT: ['pack', tkf.L33]
                    }
                }
            }
        }
    }
}


class StatusTab:
    """ GUI for the Status Tab in the System-level GUI (module_GUI.py)

    Shows Status-level information such as network information, Checkpoint information, and Fault information
    """
    def __init__(self, core, data_pipes, child_com):
        """ Widget for network, Checkpoint, and Fault interfacing

        :param core: a tkinter root
        :param data_pipes: For communication with parent System
        :param child_com: For communication with parent System
        """
        self.core = core
        self.internals: mcs.Status = data_pipes[INTERNAL]
        self.data_pipes = data_pipes
        self.child_com = child_com

        # Create dynamic memory which will be bound to the TkFire GUI
        self.memory = dict()
        self.memory['sel_checkpoint'] = tkf.tk.StringVar
        self.memory['checkpoints'] = ["None", ]

        # Create GUI
        self.display = tkf.TkFire(self.core, self.memory, mother)
        self.gui = self.display.gui

        # Populate
        self.plug_in()
        self.update()

    def update_post(func):
        """ Calls update after method

        (Relieves need to put a "self.update()" before all return statements)

        :param func: Wrapped function
        :return: Return of Wrapped function
        """
        @wraps(func)
        def wrapper_function(self, *args, **kwargs):
            try:
                return func(self, *args, **kwargs)
            finally:
                self.update()
        return wrapper_function

    def plug_in(self):
        """ Binds commands to buttons

        :return: None
        """
        self.gui['low_section!Checkpoints!B_Bypass']['command'] = lambda x='bypass': self.edit_checkpoint(x)
        self.gui['low_section!Checkpoints!B_Retry']['command'] = lambda x='retry': self.edit_checkpoint(x)
        self.gui['low_section!Checkpoints!B_Fail']['command'] = lambda x='fail': self.edit_checkpoint(x)
        self.gui['low_section!Checkpoints!B_Remove']['command'] = self.remove_checkpoint
        self.gui['low_section!Checkpoints!B_Revive']['command'] = self.revive_checkpoint

        self.gui['low_section!Fault_manager!B_Remove']['command'] = self.remove_fault
        self.gui['low_section!Fault_manager!B_Report']['command'] = self.report_fault

    def _extract_checkpoint(self, checkpoint_repr: str) -> Union[None, Tuple[str, mcs.Checkpoint]]:
        """ Converts a GUI's Checkpoint representation to an actual checkpoint

        The GUI does not display all Checkpoint information, so instead of converting the string to Checkpoint directly,
        the Checkpoint is looked up in the System's Status object.

        :param checkpoint_repr: The GUI representation of a Checkpoint --- only used to confirm that it exists
        :return: either (the task key from Status, the Checkpoint) or None
        """
        if (not checkpoint_repr) or (not isinstance(checkpoint_repr, str)) or (checkpoint_repr == "None"):
            return None
        sel_checkpoint = self.memory['sel_checkpoint'].get()  # noqa
        if (not sel_checkpoint) or (not isinstance(sel_checkpoint, str)) or (sel_checkpoint == "None"):
            return None
        task_key, *_ = sel_checkpoint.split(":", 1)
        sel_checkpoint = self.internals.get_checkpoints(_all=True).get(task_key, None)
        return task_key, sel_checkpoint

    @update_post
    def remove_checkpoint(self, *_):
        """ Removes a Checkpoint

        :return: None
        """
        sel_checkpoint = self.memory['sel_checkpoint'].get()  # noqa
        sel_checkpoint = self._extract_checkpoint(sel_checkpoint)
        if sel_checkpoint is None:
            return
        task_key, checkpoint = sel_checkpoint
        self.internals.release_checkpoint(task_key)
        self.memory['sel_checkpoint'].set("None")  # noqa
        if DEBUG:
            return
        try:
            dbi.mark_time(checkpoint.queue, checkpoint.step, checkpoint.operation, reset=True, logger=system_log)
        except (DatabaseRequestError, IndexError, TypeError) as _ce_:
            if isinstance(_ce_, DatabaseRequestError):
                _msg = "Database request bad"
            elif isinstance(_ce_, IndexError):
                _msg = "Fault poorly formatted, indices 4 and 5 should contain a queue name and step number"
            else:
                _msg = "Checkpoint not found in internals"
            system_log.warning(f"Failed to reset checkpoint {checkpoint} in the Database ({_msg})")

    @update_post
    def revive_checkpoint(self, *_):
        """ Sends a resuscitation message to the System to continue awaiting a given Checkpoint

        :return: None
        """
        sel_checkpoint = self.memory['sel_checkpoint'].get()  # noqa
        sel_checkpoint = self._extract_checkpoint(sel_checkpoint)
        if sel_checkpoint is None:
            return
        task_key, checkpoint = sel_checkpoint
        if checkpoint is None:
            return

        msg_outbox: mcn_queues.MCNPriorityQueue = self.data_pipes[MSG_Q_OUT]
        data = None
        safe_to_send = False
        while not safe_to_send:
            prompt = Defibrillator(tkf.tk.Toplevel(self.core), task_key, checkpoint, data)
            data = prompt.run()  # [location, checkpoint_key, q_name, q_step#, update_interval(, timeout)]
            if data is None:
                return
            safe_to_send = all([d != "" for d in data])
            try:
                int(data[3])
                int(data[4])
            except (ValueError, TypeError):
                safe_to_send = False
        rx_message = message_manager.Message.resuscitation_message(self.internals.name, *data)
        msg_outbox.enqueue(rx_message)

        self.memory['sel_checkpoint'].set("None")  # noqa

    @update_post
    def edit_checkpoint(self, mode: str, *_):
        """ Modifies the completion (and data) field of a Checkpoint (and consequentially triggers the Event)

        - bypass: Completion --> None
        - retry: Completion --> False, Level --> Busy
        - fail: Completion --> False, Level --> Fatal

        :param mode: How the Checkpoint is to be changed: 'bypass', 'retry', or 'fail'
        :return: None
        """
        sel_checkpoint = self.memory['sel_checkpoint'].get()  # noqa
        sel_checkpoint = self._extract_checkpoint(sel_checkpoint)
        if sel_checkpoint is None:
            return
        task_key, checkpoint = sel_checkpoint
        if checkpoint is None:
            return

        if mode == 'bypass':
            checkpoint.update_checkpoint(completion=None,
                                         data=f"Checkpoint bypassed by user in GUI")
        elif mode == 'retry':
            checkpoint.update_checkpoint(completion=False, level=mcs.V_BUSY,
                                         data=f"Checkpoint set to Busy by user to trigger retry in GUI")
        elif mode == 'fail':
            checkpoint.update_checkpoint(completion=False, level=mcs.V_FATAL,
                                         data=f"Checkpoint failed by user in GUI")

        # Let the task manager handle updating the stop/start/end times

        # Deselect the checkpoint
        self.memory['sel_checkpoint'].set("None")  # noqa

    @update_post
    def report_fault(self, *_):
        """ Opens FaultMaker UI for the creation of a Fault

        :return: None
        """
        popup = ui_fault_maker.FaultMaker(tkf.tk.Toplevel(self.core))
        reported_fault = popup.wait()

        if reported_fault:
            # Let synchronize add the Fault (or in Debug mode, add it directly)
            if DEBUG:
                self.internals.add_fault(reported_fault)
            else:
                sync.synchronized_fault_add(self.data_pipes, reported_fault)

    @update_post
    def remove_fault(self, *_):
        """ Removes a Fault

        :return:
        """
        cursor = self.gui['low_section!Fault_list_frame!Fault_list'].curselection()
        if not cursor:
            return
        for index in cursor:
            selected_fault_repr: str = self.gui['low_section!Fault_list_frame!Fault_list'].get(index)
            if not selected_fault_repr:
                continue
            selected_fault = mcs.Fault.loads(selected_fault_repr)
            system_log.debug(f"(GUI) Attempting to remove: {selected_fault}")

            # Let synchronize handle the Fault removal, unless in Debug mode
            if DEBUG:
                self.internals.remove_fault(selected_fault)
            else:
                sync.synchronized_remove_fault(self.data_pipes, selected_fault)

    def update(self):
        """ Updates GUI elements to reflect their current values

        :return: None
        """
        # Network
        conn, mode = self.internals.get_network_state()
        if (mode == mcs.V_ONLINE) and (conn == mcs.V_ONLINE):
            net_state = "Online"
        elif (mode == mcs.V_OFFLINE) and (conn == mcs.V_OFFLINE):
            net_state = "Offline"
        elif (mode == mcs.V_ONLINE) and (conn == mcs.V_OFFLINE):
            net_state = "Disconnected"
        elif (mode == mcs.V_OFFLINE) and (conn == mcs.V_ONLINE):
            net_state = "Misconnected"
        else:
            net_state = "Unknown"
        self.gui['top_section!Network!Connectivity']['text'] = net_state
        ip, port = self.internals.get_ip_and_port()
        self.gui['top_section!Network!IP_and_Port']['text'] = f"{ip} @{port}"

        # Members
        area_members = self.internals.get(mcs.DN_MEMBERS, list())
        with self.child_com as cc:  # noqa
            local_members = list(cc.keys())
        local_members = list_subtract(local_members, ["ui", 'is_initialized'])

        self.gui['top_section!Members!Area_members']['text'] = " ".join(area_members)
        self.gui['top_section!Members!Local_members']['text'] = " ".join(local_members)

        # Checkpoints
        self.memory['checkpoints'] = [f"{tk}: {ch.location}, {ch.level}, {ch.data}"
                                      for tk, ch in self.internals.get_checkpoints(_all=True).items()]
        if not self.memory['checkpoints']:
            self.memory['checkpoints'] = ["None", ]
        checkpoint_optionmenu_widget = self.gui['low_section!Checkpoints!Selector']
        menu = checkpoint_optionmenu_widget['menu']
        menu.delete(*tkf.INPUT_CONTENTS)
        for ch in self.memory['checkpoints']:
            menu.add_command(label=ch,
                             command=lambda value=ch: self.memory['sel_checkpoint'].set(value))  # noqa

        # Faults
        self.gui['low_section!Fault_list_frame!Fault_list'].delete(0, tkf.tk.END)
        all_faults = [repr(f) for f in sorted(self.internals.get_faults(_all=True))]
        for f in all_faults:
            self.gui['low_section!Fault_list_frame!Fault_list'].insert(tkf.tk.END, f)

    update_post = staticmethod(update_post)


if __name__ == '__main__':
    from thread_safe_data import ThreadSafeDataContainer
    from mcn_queues import MCNPriorityQueue

    DEBUG = True

    root = tkf.tk.Tk()

    mailbox = MCNPriorityQueue()
    my_state = mcs.Status("MC")
    _child_com = ThreadSafeDataContainer(dict())
    with _child_com as cc:
        cc['ui'] = {mcs.C_ACTIVE_MODE: False,
                    mcs.C_QUEUE_BLOCK: list(),
                    mcs.C_SAFE_MODE: False,
                    mcs.C_EXPIRED: [False, ]}
        cc['is_initialized'] = True
        cc['subsystem 1'] = dict()
        cc['subsystem 2'] = dict()
        cc['subsystem 3'] = dict()

    my_state.add_fault(mcs.Fault("_U", mcs.V_BUSY, "Test fault 1"))
    my_state.add_fault(mcs.Fault("_U", mcs.V_PROBLEM, "Test fault 2"))
    my_state.add_fault(mcs.Fault("_U", mcs.V_FATAL, "Test fault 3"))
    my_state.add_fault(mcs.Fault("_U", mcs.V_BUSY, "Test fault 4"))

    my_state.add_checkpoint("123abc", mcs.Checkpoint(None, "MC", None, "Test Checkpoint 1", "Test_queue_1", "5"))
    my_state.add_checkpoint("456ijk", mcs.Checkpoint(None, "AH", None, "Test Checkpoint 2", "Test_queue_2", "13"))
    my_state.add_checkpoint("789lmn", mcs.Checkpoint(None, "SP", None, "Test Checkpoint 3", "Test_queue_3", "7"))
    my_state.add_checkpoint("123xyz", mcs.Checkpoint(None, "LC", None, "Test Checkpoint 4", "Internal", "Sync"))
    my_state.add_checkpoint("000pqr", mcs.Checkpoint(None, "__", None, "Test Checkpoint 5", "Internal", "Proxy"))

    my_gui = StatusTab(root, {INTERNAL: my_state, MSG_Q_OUT: mailbox, MSG_Q_IN: mailbox}, _child_com)

    root.mainloop()

    print(my_state)
    while not mailbox.empty():
        print(mailbox.get_nowait().print())
