""" Main GUI elements for a System

@author: Ben C
"""

import queue
import time
import tkinter as tk
from copy import deepcopy
from io import StringIO
from itertools import groupby
from tkinter import scrolledtext as st
from tkinter import ttk

import yaml

import mcn_queues
import mcn_status as mcs
import message as msgm
import operations as oprtn
import ui_database_manager
import ui_message_maker
import ui_status_tab
from constants import *
from gui_constants import *
from mcn_logging import level_translator, LOGGER_LEVELS
from mcn_logging_manager import system_log
from ui_quick_dialog import QuickUI

# Import system-specific UI elements
try:
    import ui_task_mananger
    import ui_hplc_volumes
    import ui_lh_manager
    import ui_spark_state
except (ImportError, OSError, ):
    system_log.exception(f"{__name__}: Failed to load a module/component!")


def generate_main_menu_text(system_name=""):
    leaflet = ""
    if system_name == "MC":
        leaflet = """
  * "Open Task Manager" will open a dialog to control how the MC schedules and runs
    - Only active on the MC
    - To make edits, database operation must be paused
    - Control:
      > Database control can be started, paused, or resumed
      > Pulling the database populates queues for selection
    - Settings:
      > By selecting a queue, you restrict operation to only that queue (use allow all to allow all queues to be run)
    - Modus Operandi:
      > Automatic mode - Scheduler will move freely from one task to another
      > Safe mode - Scheduler will ask the user before starting a task
    - Select Group:
      > You may move queues between groups, first select which group to view and the queues within this group will
        be presented in Select Queue
    - Select Queue:
      > With a group selected, you may select a queue, this will enable the buttons in Move Queue
    - Move Queue:
      > Select to which group you wish to move the selected queue"""
    elif system_name == "AH":
        leaflet = """
  * "Open Lh Queue Manager" will open a dialog to control the liquid handler
    - Allows the user to Reset or Stop the instrument
      > Useful for error recovery
    - Also allows for custom signals to be sent (use at own risk)"""
    elif system_name == "LC":
        leaflet = """
  * "Open Instrument Manager" opens a dialog to control the column and volume (left)
    - Additional controls for wellplate and project selection and Batch stop/resume commands (right)
    - Feedback provided in text boxes below the left and right sections, respectively"""
    elif system_name == "NM":
        leaflet = """
  * "Open Spinner Manager" will open a dialog to control the NMR
    - Not implemented"""
    elif system_name == "SP":
        leaflet = """
  * "Open LPX Manager" will open a dialog to control the LPX
    - Only active on SP
    - Provides the user with a form to lookup hotels on the LPX and to edit occupancy
  * "Open Spark Manager" will open a dialog to control the Spark
    - Only active on SP
    - Presents the state of the Spark and allows the user to change the state as needed (functionality is questionable)
    - Also presents controls for the mechanical operation of the Spark and its peripheral devices
      > Tray movement
      > Lamp
      > Heater
      > Gas-flow
      > Photo-reactor lift
  * "Open Thermoreactor Manager" will open a dialog to control the thermo-reactor
    - Only active on SP
    - Provides the ability to mutate the state object
    - Provides control over all subsystem functions"""

    return f"""Welcome to the Master Controller Network (MCN) User Interface (UI)
This tool is designed to help with the controlling and monitoring of systems on the AMD platform.
This main tab will provide an overview of functionality for this UI

The Controls tab will allow you to issue command and log notes.
  * "Create Message" will open a dialog to create a custom message which can be sent over the MCN
    - You may spoof an alternate sender if you so wish by entering the name of the system in the sender field
    - You may select the recipient from a list of (sub)systems (Note "_A" means "to all")
    - All messages have ID keys, these can be automatically generated or a special bypass key can be used
    - Message Contents specify how the message is to be read by the recipient
    - Since the UI has no checkpoints, the priority can only be normal or high
    - The data field will be populated depending on the message contents
      > For Args (almost never used) and Kwargs, enter the data as a json dictionary (note this requires double quotes)
  * "(Dis)connect" is in reference to the network (it is fine to "disconnect" from a network multiple times)
  * "Network Settings" allows the user to change the IP and port of the MCN{leaflet}
  * "Shutdown" exits the GUI and will close out the system
  * "Kill Server" will shut down the server
  * "Emergency Shutdown All" calls Shutdown and Kill Server
  * "Open DB Manager" opens a dialog to allow quick edits and reads from the database
    - In the normal manager, "Pull DB" refreshes the list of queues available and "Select Queue" selects a queue for
      abbreviated listing.  By selecting an item in the list, further details will be presented.
    - You may mark steps as complete or incomplete manually (the edit will appear in the details but not in the list)
    - The "Advanced Editor" button opens a dialog which allows you direct access to the database.
      > The large text edit box can be edited and used to "save" documents to registers--used for updates and adds
      > The "Reopen Document" button will pull the (add), (new), or (previous) document (in that hierarchy)
      > The "Unload Saved Documents" clears all saved documents
      > When adding or editing documents, don't forget to trim the ID header (but not the _id field)
      > "Reformat" will clean the display (document-->yaml-->document)
      > When using "query_document", the "What's here?" has two behaviors depending on the contents of "search_field"
        * If "search_field" is empty, the top-level of the collection will be printed (e.g. queue names for the queue 
          collection)
        * If "search_field" contains an expression "key::criterion", the top-level of the collection which meets the
          criterion will be printed (e.g. queue names which contain the word 'sulfur' for the queue collection)
          -  "key" is a dot-flattened dictionary key (e.g. "toplevel.level_one.field")
          -  "criterion" is an eval()-compliant string which returns a boolean (e.g. "_val_ < 100") and which can use 
             the keyword "_val_" to represent the value keyed by "key" in the document
    - The "Fault Editor" allows the user to view and modify faults
      > The user may also delete a Fault
      > The user may re-attribute (change location field of) a Faults (incuding claiming "User Error")
  * "Shutdown" exits the GUI and will close out the system
  * "Kill Server" will shut down the server
  * "Emergency Shutdown All" calls Shutdown and Kill Server
  * "Open DB Manager" opens a dialog to allow quick edits and reads from the database
    - In the normal manager, "Pull DB" refreshes the list of queues available and "Select Queue" selects a queue for
      abbreviated listing.  By selecting an item in the list, further details will be presented.
    - You may mark steps as complete or incomplete manually (the edit will appear in the details but not in the list)
    - The "Advanced Editor" button opens a dialog which allows you direct access to the database.
      > The large text edit box can be edited and used to "save" documents to registers--used for updates and adds
      > The "Reopen Document" button will pull the (add), (new), or (previous) document (in that hierarchy)
      > The "Unload Saved Documents" clears all saved documents
      > When adding or editing documents, don't forget to trim the ID header (but not the _id field)
      > "Reformat" will clean the display (document-->yaml-->document)
  
The Status tab contains information on the (sub)systems.
  * Network connectivity will show the IP and PORT of the connection
  * The State section will allow you to view faults and checkpoints for a system and its subsystems
    - Checkpoints can be selected and bypassed, retried, failed, or removed
      > Bypass - marked as successful
      > Retry - marked as busy
      > Fail - marked as fatal
      > Remove - deleted
      > Revive - restart monitoring thread
  * The Fault manager can select a fault (regardless of which (sub)system is selected in State) and remove it.
    It can also create faults ("Report Fault") which will open a dialog that asks for the location, level, and Info.
      > Level should be one of the pre-defined values ('Good', 'Busy', 'Problem', 'Fatal'), but it allowed to be set 
        to any value.
"""


class MCN_GUI:  # noqa # This is an acronym
    def __init__(self, data_pipes, internals, cmd_list_keys):
        """ Constructor for a GUI element

        :param data_pipes: Used to access properties of the associated System
        :param internals: Used to access properties of the associated System
        :param cmd_list_keys: Used to know what functionality the associated System has
        """
        self._data_pipes = data_pipes
        self.inbox = data_pipes[MSG_Q_IN]
        self.outbox = data_pipes[MSG_Q_OUT]
        self.child_com = data_pipes[CHILD_COM]
        self.internals: mcs.Status = internals
        self.name = self.internals.name
        self.cmd_list_str = list(cmd_list_keys)
        self.core = tk.Tk()
        self.core.withdraw()
        self.tabs = dict()
        self.frames = dict()
        self.widgets = dict()
        self.ui_states = dict()
        self.root = None
        self.nb = None
        self.dialog = dict()
        self.exit_by_shutdown = False
        # Each System has some functionality/subsystem unique to it, so add GUI controllers if present
        self.special_cmd = {"MC": [("Open Task Manager", self.open_task_manager), ],
                            "SP": [("Open LPX Manager", self.open_ss_instrument_manager),
                                   ("Open Spark Manager", self.open_pr_instrument_manager),
                                   ("Open Thermoreactor Manager", self.open_th_instrument_manager), ],
                            "AH": [("Open Lh Queue Manager", self.open_liquid_handler_manager), ],
                            "LC": [("Open Instrument Manager", self.lc_advanced_gui), ]}

        self.build_gui()
        self.core.after(0, self._refresh)

    def build_gui(self):
        """ Calls constructors for the notebook widget and child frames

        :return: None
        """
        # Build the notebook super-structure
        self.build_notebook()
        # Build frames for each tab
        self.build_frames()
        # Lock until initialized
        while True:
            with self.child_com as cc:
                if cc['is_initialized']:
                    break
                time.sleep(1)
        self.widgets['ctrls.direct_cmds']["(dis)connect"]['state'] = tk.ACTIVE
        if self.internals.get_network_state()[1] == mcs.V_ONLINE:
            self.widgets['ctrls.direct_cmds']["(dis)connect"]['text'] = "Disconnect"
        else:
            self.widgets['ctrls.direct_cmds']["(dis)connect"]['text'] = "Connect"
        self.nb.enable_traversal()

    def _refresh(self):
        """ Periodic refresh for state elements

        :return: None
        """

        try:
            self.core.winfo_exists()
        except tk.TclError:
            return

        if self.internals.get_network_state()[1] == mcs.V_ONLINE:
            self.widgets['ctrls.direct_cmds']["(dis)connect"]['text'] = "Disconnect"
        else:
            self.widgets['ctrls.direct_cmds']["(dis)connect"]['text'] = "Connect"

        self.widgets['status.edit'].update()

        self.core.after(5140, self._refresh)

    def build_notebook(self):
        """ Constructs the notebook widget and its tabs, does not populate them

        :return: None
        """
        self.root = tk.Toplevel()
        self.root.title(f"{self.name}: Control Panel")
        self.nb = ttk.Notebook(self.root)
        self.nb.pack(expand=True, fill='both')
        # Main tab
        self.tabs["main"] = tk.Frame(self.nb)
        self.nb.add(self.tabs["main"], text="Main")
        # Tab 2: Controls
        self.tabs["ctrls"] = tk.Frame(self.nb)
        self.nb.add(self.tabs["ctrls"], text="Controls")
        # Tab 3: Status
        self.tabs["status"] = tk.Frame(self.nb)
        self.nb.add(self.tabs["status"], text="Status")
        # Tab 4: Check
        self.tabs["debug"] = tk.Frame(self.nb)
        self.nb.add(self.tabs["debug"], text="Debug")
        #
        self.nb.select(self.tabs["main"])
        # self.nb.enable_traversal()

    def build_frames(self):
        """ Populates frames for each notebook tab

        :return: None
        """
        # Main Frame                                                                                               ~~~~#
        self.frames['main'] = dict()
        self.frames['main']['help'] = st.ScrolledText(self.tabs["main"],
                                                      wrap=tk.WORD,
                                                      width=120,
                                                      height=29)
        self.frames['main']['help']['state'] = tk.NORMAL
        self.frames['main']['help'].insert(tk.INSERT, generate_main_menu_text(self.name))
        self.frames['main']['help']['state'] = tk.DISABLED
        self.frames['main']['help'].pack(expand=True, fill='both')

        # Controls Frame                                                                                           ~~~~#
        self.frames['ctrls'] = dict()
        self.frames['ctrls']['direct_cmds'] = tk.LabelFrame(self.tabs["ctrls"], text="Direct Commands")
        self.frames['ctrls']['direct_cmds'].pack(**LB33)
        self.frames['ctrls']['logger'] = tk.LabelFrame(self.tabs['ctrls'], text="Logger")
        self.frames['ctrls']['logger'].pack(expand=True, **LB33)

        self.widgets['ctrls.direct_cmds'] = dict()
        self.widgets['ctrls.direct_cmds']["create_msg"] = tk.Button(self.frames['ctrls']['direct_cmds'],
                                                                    text="Create Message",
                                                                    command=self.open_create_message_dialog)
        self.widgets['ctrls.direct_cmds']["(dis)connect"] = tk.Button(self.frames['ctrls']['direct_cmds'],
                                                                      text="Connect",
                                                                      command=self.change_connectivity,
                                                                      state=tk.DISABLED)
        self.widgets['ctrls.direct_cmds']["network"] = tk.Button(self.frames['ctrls']['direct_cmds'],
                                                                 text="Network Settings",
                                                                 command=self.edit_network)
        self.widgets['ctrls.direct_cmds']["spec_control"] = dict()
        for i, w in enumerate(self.special_cmd.get(self.name, [])):
            try:
                self.widgets['ctrls.direct_cmds']["spec_control"][i] = tk.Button(self.frames['ctrls']['direct_cmds'],
                                                                                 text=w[0],
                                                                                 command=w[1])
            except (KeyError, IndexError):
                self.widgets['ctrls.direct_cmds']["spec_control"][i] = tk.Button(self.frames['ctrls']['direct_cmds'],
                                                                                 text="<Not Implemented>")
                self.widgets['ctrls.direct_cmds']["spec_control"][i]['state'] = tk.DISABLED

        self.widgets['ctrls.direct_cmds']["shutdown"] = tk.Button(self.frames['ctrls']['direct_cmds'],
                                                                  text="Shutdown",
                                                                  command=self.shutdown)
        self.widgets['ctrls.direct_cmds']["kill_server"] = tk.Button(self.frames['ctrls']['direct_cmds'],
                                                                     text="Kill Server",
                                                                     command=self.kill_server)
        self.widgets['ctrls.direct_cmds']["emergency_stop"] = tk.Button(self.frames['ctrls']['direct_cmds'],
                                                                        text="Emergency Shutdown All",
                                                                        command=self.emergency_stop)
        self.widgets['ctrls.direct_cmds']["database_mngr"] = tk.Button(self.frames['ctrls']['direct_cmds'],
                                                                       text="Open DB Manager",
                                                                       command=self.manage_db)
        for button in self.widgets['ctrls.direct_cmds'].keys():
            selection = self.widgets['ctrls.direct_cmds'][button]
            if isinstance(selection, dict):
                for _button in selection:
                    selection[_button].pack(pady=1, **TB33)
            else:
                selection.pack(pady=1, **TB33)

        self.widgets['ctrls.logger'] = dict()
        self.ui_states['logger'] = dict()
        self.ui_states['logger']['log_level'] = tk.StringVar(value="notset")
        self.widgets['ctrls.logger']['log_level'] = tk.OptionMenu(self.frames['ctrls']['logger'],
                                                                  self.ui_states['logger']['log_level'],
                                                                  *LOGGER_LEVELS)
        self.widgets['ctrls.logger']['log_level'].pack(**TB33)
        self.widgets['ctrls.logger']['log_entry'] = st.ScrolledText(self.frames['ctrls']['logger'],
                                                                    wrap=tk.WORD,
                                                                    width=100,
                                                                    height=9)
        self.widgets['ctrls.logger']['log_entry'].pack(expand=True, **TB33)
        self.widgets['ctrls.logger']['log_submit'] = tk.Button(self.frames['ctrls']['logger'],
                                                               text="Submit",
                                                               command=self.user_log)
        self.widgets['ctrls.logger']['log_submit'].pack(**TB33)

        # Status Frame                                                                                             ~~~~#
        self.frames['status'] = dict()
        self.frames['status']['state'] = tk.LabelFrame(self.tabs['status'], text="State")
        self.frames['status']['state'].pack(expand=True, **TB33)

        self.widgets['status.edit'] = ui_status_tab.StatusTab(self.frames['status']['state'],
                                                              self._data_pipes,
                                                              self.child_com)

        self.frames['debug'] = dict()
        self.widgets['debug'] = dict()
        self.frames['debug']['text'] = st.ScrolledText(self.tabs["debug"],
                                                       wrap=tk.WORD,
                                                       width=120,
                                                       height=29)
        self.frames['debug']['text'].pack(expand=True, **TB33)
        if self.name == 'MC':
            tk.Button(self.tabs["debug"], text="Refresh", command=self._debug).pack(**LB33)
            tk.Button(self.tabs["debug"], text="Scheduler", command=self.schedule).pack(**LB33)
        else:
            tk.Button(self.tabs["debug"], text="Refresh", command=self._debug).pack(**TB33)

    def run(self):
        """ Calls mainloop and handles shutdown on mainloop escape

        :return: A boolean representing if the GUI is exiting due to a shutdown call (opposed to crash/interrupt)
        """
        self.core.protocol("WM_DELETE_WINDOW", self.on_closing)
        self.core.mainloop()
        return self.exit_by_shutdown

    def on_closing(self):
        """ Log shutdown and call shutdown

        :return: None
        """
        # Disconnect, Send a signal that this system is shutting down, cleanup
        system_log.debug(f"MCN_GUI({self.name}).on_closing() called :)")
        with self.child_com as cc:
            if cc['ui'][mcs.C_ACTIVE_MODE]:
                system_log.mark_down()
        self.shutdown()

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    def open_create_message_dialog(self):
        """ Spawns the MessageMaker GUI to allow the user to create and send their own messages

        :return: None
        """
        self.dialog["msg"] = ui_message_maker.MessageMaker(tk.Toplevel(self.core),
                                                           self._data_pipes,
                                                           self.cmd_list_str,
                                                           self.name)

    def change_connectivity(self):
        """ Attempts to change the network mode, and updates the GUI to respect the connectivity.

        :return: None
        """
        conn, mode = self.internals.get_network_state()
        if mode == mcs.V_ONLINE:
            system_log.debug("Network mode set FALSE by user")
            self.internals.network_mode_is(False)
        else:
            system_log.debug("Network mode set TRUE by user")
            self.internals.network_mode_is(True)
        self.widgets['ctrls.direct_cmds']["(dis)connect"]['state'] = tk.DISABLED
        time.sleep(5)
        conn, mode = self.internals.get_network_state()
        if mode == mcs.V_ONLINE:
            self.widgets['ctrls.direct_cmds']["(dis)connect"]['text'] = "Disconnect"
        else:
            self.widgets['ctrls.direct_cmds']["(dis)connect"]['text'] = "Connect"
        self.widgets['status.edit'].update()
        self.widgets['ctrls.direct_cmds']["(dis)connect"]['state'] = tk.ACTIVE

    def edit_network(self):
        """ Opens a QuickUI to allow the user to change the IP and PORT of the MCN

        :return: None
        """
        popup = QuickUI(tk.Toplevel(self.core),
                        title="Network Settings",
                        dialog="Please enter the IP and PORT\nFormat 'IP:PORT'",
                        buttons={},
                        has_entry=True)
        _ip, _port = self.internals.get_ip_and_port()
        popup.entry_field.insert(0, f"{_ip}:{_port}")
        new_network_settings = popup.run()
        print(f"Popup returned: '{new_network_settings}'")
        if not new_network_settings:
            return
        try:
            new_network_settings = new_network_settings.replace(" ", "")
            new_ip, new_port = new_network_settings.split(":")
            new_port = int(new_port)
        except ValueError:
            system_log.info(f"User attempted to change network settings to '{new_network_settings}'. "
                            f"Could not be parsed as 'IP:PORT'")
        else:
            self.internals.modify_network_settings(new_ip, new_port)

    def open_task_manager(self):
        """ (MC only) Opens a TaskManager GUI to allow Database Control

        :return:
        """
        self.dialog["task"] = ui_task_mananger.TaskManager(tk.Toplevel(self.core),
                                                           self._data_pipes,
                                                           self.widgets['ctrls.direct_cmds']["spec_control"][0])

    def open_ss_instrument_manager(self):
        """ (SP only) Opens the LPX (storage carousel) GUI to allow resource management

        :return:
        """
        self.widgets['ctrls.direct_cmds']["spec_control"][0]['state'] = tk.DISABLED
        try:
            inbox: MCNPriorityQueue = self._data_pipes[MSG_Q_IN]
            inbox.enqueue(msgm.Message.build_from_args(_type="CMD",
                                                       _sender="SP",
                                                       _recip="Ss",
                                                       _contents="ib.Read_Operation",
                                                       _data=oprtn.Operation({
                                                           oprtn.FUNC: "debug.run_lpx_gui",
                                                           oprtn.AGENT: "Ss"
                                                       }).package(),
                                                       _priority=msgm.P_CMD_NR
                                                       )
                          )
        finally:
            self.widgets['ctrls.direct_cmds']["spec_control"][0]['state'] = tk.NORMAL

    def open_pr_instrument_manager(self):
        """ (SP only) Opens the SparkStateUI to allow the user to change the state of the plate reader and associated
        photo-reaction tools

        :return: None
        """
        try:
            with self.child_com as cc:
                spark = cc['Pr']['SparkControl'].spark
        except (KeyError, AttributeError, TypeError):
            spark = None
        self.dialog["pr"] = ui_spark_state.SparkStateUI(tk.Toplevel(self.core),
                                                        spark,
                                                        self._data_pipes[MSG_Q_OUT],
                                                        self.widgets['ctrls.direct_cmds']["spec_control"][1],
                                                        )

    def open_th_instrument_manager(self):
        """ (SP only) Opens the Thermo-reactor UI to allow the user to change the state of the thermoreactor and
        associated thermo-reactor tools

        :return: None
        """
        self.widgets['ctrls.direct_cmds']["spec_control"][2]['state'] = tk.DISABLED
        try:
            inbox: MCNPriorityQueue = self._data_pipes[MSG_Q_IN]
            inbox.enqueue(msgm.Message.build_from_args(_type="CMD",
                                                       _sender="SP",
                                                       _recip="Th",
                                                       _contents="ib.Read_Operation",
                                                       _data=oprtn.Operation({
                                                           oprtn.FUNC: "debug.run_th_gui",
                                                           oprtn.AGENT: "Th"
                                                       }).package(),
                                                       _priority=msgm.P_CMD_NR
                                                       )
                          )
        finally:
            self.widgets['ctrls.direct_cmds']["spec_control"][2]['state'] = tk.NORMAL

    def open_liquid_handler_manager(self):
        """ (AH only) Opens an EvowareHandler GUI to allow modifications to the communications pipe (queue)

        :return: None
        """
        self.dialog["db"] = ui_lh_manager.EvowareHandler(tk.Toplevel(self.core),
                                                         self._data_pipes[CHILD_COM],
                                                         self.widgets['ctrls.direct_cmds']["spec_control"][0])

    def lc_advanced_gui(self):
        """ (LC only) opens the EllseeManager GUI to allow the user to update column and solvent details

        :return:
        """
        child_com = self._data_pipes[CHILD_COM]
        with child_com as cc:
            ls_object = cc['shimadzu_wrapper']
        self.dialog['lc'] = ui_hplc_volumes.EllseeManager(tk.Toplevel(self.core),
                                                          ls_object,
                                                          self.widgets['ctrls.direct_cmds']["spec_control"][0])

    def shutdown(self):
        """ Destroys GUI and sets shutdown flag

        :return:
        """
        system_log.info("User has issued shutdown")
        self.core.destroy()
        self.exit_by_shutdown = True

    def kill_server(self):
        """ Sends a message to the Server to shut down

        :return:
        """
        self.outbox.enqueue(msgm.Message.close_server())

    def emergency_stop(self):
        """ (Not Implemented)  Currently same as kill_server(), but also closes GUI without setting shutdown flag.

        :return:
        """
        self.outbox.enqueue(msgm.Message.close_server())
        self.core.destroy()

    def manage_db(self):
        """ Opens the DatabaseManager GUI to allow the user to view and edit the database

        :return:
        """
        self.dialog["db"] = ui_database_manager.DatabaseManager(tk.Toplevel(self.core),
                                                                self.name)

    def user_log(self):
        """ Commits a user-made message to the log

        :return:
        """
        level = self.ui_states['logger']['log_level'].get()
        message = self.widgets['ctrls.logger']['log_entry'].get(*ST_CONTENTS)
        system_log.log(level_translator(level), message)
        self.widgets['ctrls.logger']['log_entry'].delete(*ST_CONTENTS)

    def _debug(self):
        """ Updates the debug page with state data to aid in troubleshooting and debugging

        :return: None
        """
        # This is all to get YAML to behave
        proxy = dict()
        proxy[MSG_Q_IN] = len(list(self._data_pipes[MSG_Q_IN].queue))
        proxy[MSG_Q_OUT] = len(list(self._data_pipes[MSG_Q_OUT].queue))
        proxy[INTERNAL] = deepcopy(self._data_pipes[INTERNAL].get_short_copy())
        proxy[INTERNAL][mcs.S_FAULTS] = [f.dumps() for f in proxy[INTERNAL][mcs.S_FAULTS]]
        proxy[INTERNAL][mcs.S_CHECKPOINTS] = {k: repr(v) for k, v in proxy[INTERNAL][mcs.S_CHECKPOINTS].items()}

        cc = self._data_pipes[CHILD_COM].__dict__()["ThreadSafeDataContainer"]
        proxy[CHILD_COM] = dict()
        proxy['ui'] = dict()
        proxy['ui'][mcs.C_ACTIVE_MODE] = bool(cc['ui'][mcs.C_ACTIVE_MODE])
        proxy['ui'][mcs.C_QUEUE_BLOCK] = list(cc['ui'][mcs.C_QUEUE_BLOCK])
        proxy['ui'][mcs.C_SAFE_MODE] = bool(cc['ui'][mcs.C_SAFE_MODE])
        proxy['ui'][mcs.C_EXPIRED] = [bool(cc['ui'][mcs.C_EXPIRED][0]), ]
        proxy['is_initialized'] = bool(cc['is_initialized'])

        stream = StringIO()
        yaml.safe_dump(proxy, stream)
        text = stream.getvalue()
        del stream

        self.frames['debug']['text']['state'] = tk.NORMAL
        self.frames['debug']['text'].delete(*ST_CONTENTS)
        self.frames['debug']['text'].insert(tk.INSERT, text)
        self.frames['debug']['text']['state'] = tk.DISABLED

    def schedule(self):
        """ Similar to _debug() but with task-management/scheduling details

        :return:
        """
        child_com = self._data_pipes[CHILD_COM]
        with child_com as cc:
            q: mcn_queues.LogQueue = cc['ui'][mcs.C_TASK_QUEUE]
        q_list = list(q.queue)
        q_list = [key for key, _group in groupby(q_list)]  # remove consecutive duplicates
        if q_list:
            q_list = [str(i) for i in q_list]
            text = "\n".join(q_list)
        else:
            text = "<Scheduling Log moved to task_manager.log in Logs folder>"

        self.frames['debug']['text']['state'] = tk.NORMAL
        self.frames['debug']['text'].delete(*ST_CONTENTS)
        self.frames['debug']['text'].insert(tk.INSERT, text)
        self.frames['debug']['text']['state'] = tk.DISABLED
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing


if __name__ == '__main__':
    from mcn_queues import MCNPriorityQueue
    from thread_safe_data import ThreadSafeDataContainer

    my_state = mcs.Status("MC")
    my_state.add_checkpoint("test_ch", mcs.Checkpoint(False, "Pr", None, "waiting for checkpoint", "queue_name"))
    boxes = {MSG_Q_OUT: MCNPriorityQueue(),
             MSG_Q_IN: MCNPriorityQueue(),
             CHILD_COM: ThreadSafeDataContainer({'ui': {mcs.C_ACTIVE_MODE: False, mcs.C_QUEUE_BLOCK: list(),
                                                        mcs.C_SAFE_MODE: False, mcs.C_EXPIRED: [False, ],
                                                        mcs.C_TASK_QUEUE: mcn_queues.LogQueue(29)},
                                                 'is_initialized': True}),
             INTERNAL: my_state}
    with boxes[CHILD_COM] as _cc:
        for j in range(0, 100):
            print(f"putting {j}")
            _cc['ui'][mcs.C_TASK_QUEUE].put(j)
    with boxes[CHILD_COM] as _cc:
        _cc['Lh'] = queue.Queue()
        _cc['Lh'].put(["IDLE", ])
        _cc['shimadzu_wrapper'] = ui_hplc_volumes.DummyLC()
        _cc['Pr'] = dict()


        class Dum:
            def __getattr__(self, item):
                if item == 'state':
                    return {
                        "tray_position": "In",
                        "thermostat": None,
                        "encumbered": False,
                        "expose": False,
                    }
                return self

        _cc['Pr']['SparkControl'] = Dum()
    my_gui = MCN_GUI(boxes, my_state, ["1abc", "2xyz"])
    my_state.modify_network_settings('localhost', 8080)

    my_gui.run()

    try:
        print(boxes[MSG_Q_OUT].get(block=False).print())
    except queue.Empty:
        pass
    try:
        print(boxes[MSG_Q_IN].get(block=False).print())
    except queue.Empty:
        pass

    with boxes[CHILD_COM] as chcm:
        print(chcm['ui'])
        print(chcm['Lh'].queue)
