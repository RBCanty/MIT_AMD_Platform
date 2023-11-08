""" System

Workhorse of the MCN

Executor of all the sub-processes which execute operations on the computer and interface with the platform

Class Fields
  - name                   : Name of the system
  - data_pipes             : A collection of thread-safe data pipes for shared resources
  - status                 : A special data pipe for system status/state information
  - core_process_threads   : Holds thread objects spawned by the System which are required for operation
  - core_process_functions : Holds the functions which give the System its core functionality
  - child_processes        : Holds thread objects spawned by the System for message-based tasks
  - itinerary              : (MC only) A checklist of tasks to be performed on the platform
  - kill_signal_received   : A poison-pill flag which indicates if the system needs to shut down
  - received_keys          : A list of idempotency keys, incoming messages checked against this list
                             (cleaned of keys older than 5 days)
  - command_list           : A dictionary mapping keyword strings to python functions

Class Threads
  - user_com()          : Manages the GUI and user I/O
  - network_com()       : Manages network I/O via a Client

  - processor_intake()  : Spawn child processes to complete the tasks issued by the incoming messages
  - processor_cleanup() : Delete completed processes

  - monitor()           : (MC only) Routinely generates requests for status/state updates to send to all systems
  - scheduler()         : (MC only) Calls scheduling.Coordinator to try to coordinate tasks on the platform

Class Methods
  - build_command_list()     : Helper function, imports the System's methods and adds them tot the command list

  - run()                    : Spawns the core threads and joins the sentinel
  - sentinel()               : Watches threads and restarts them when they terminate/crash

  - validate_key()           : Checks an incoming idempotency key against the list of keys already received
  - clean_keyring()          : Clears keys that are ~5 days old

@author: Ben C
"""

import importlib
import time
from datetime import date, timedelta
from pprint import pformat
from threading import Thread

from mcn_logging_manager import system_log

from mc_client import MCClient
import functionals
import mcn_queues
import message as msgm
import module_GUI
import aceso
import task_management
from constants import *
import mcn_status as mcs
from thread_safe_data import ThreadSafeDataContainer


class System:
    """ Workhorse of the MCN, the main processor of MCN information

    A System spawns threads to perform tasks specified by messages from the MCN
    """

    def __init__(self, name, *, skip_initialization=False, debug_mode=False):
        """ Creates a new System instance

        :param name: Name of the system
        :param skip_initialization: (Default False, optional parameter, keyword only) If True, the initialize method
        will not be run
        :param debug_mode: (Default False, optional parameter, keyword only) If True, ALL system processes will ask for
        user confirmation before execution (not just queue/database operations)
        """
        #
        self.name = name
        self.status = mcs.Status(name)
        #
        self.data_pipes = dict()
        self.data_pipes[MSG_Q_IN] = mcn_queues.MCNPriorityQueue()
        self.data_pipes[MSG_Q_OUT] = mcn_queues.MCNPriorityQueue()
        self.data_pipes[CHILD_COM] = ThreadSafeDataContainer(dict())
        self.data_pipes[INTERNAL] = self.status
        #
        self.core_process_threads = []
        self.core_process_functions = [self.user_com,
                                       self.network_com,
                                       self.processor_intake,
                                       self.processor_cleanup,
                                       self.monitor]
        self.child_processes = ThreadSafeDataContainer(list())
        if self.name == 'MC':
            self.core_process_functions.append(self.scheduler)
        with self.child_processes as cp:
            cp += aceso.anabasis(r"./orpheus.json", self.data_pipes, functionals.resuscitate_queue_operation)
        #
        self.client = MCClient(self.name, self.data_pipes, *self.status.get_ip_and_port())
        self.gui = None
        self.received_keys = ThreadSafeDataContainer(list())
        self.command_list = dict()
        self.build_command_list()
        #
        self.kill_signal_received = False
        self.debug_mode = debug_mode
        self.debug_memory = dict()
        #
        with self.data_pipes[CHILD_COM] as cc:
            cc['ui'] = {mcs.C_ACTIVE_MODE: False,
                        mcs.C_QUEUE_BLOCK: list(),
                        mcs.C_SAFE_MODE: False,
                        mcs.C_EXPIRED: [False, ],
                        mcs.C_TASK_QUEUE: mcn_queues.LogQueue(29)}
            cc['is_initialized'] = skip_initialization

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    def build_command_list(self):
        """ Constructor for the System's command list

        Imports from module specified by the system's name

        :return: A dictionary of {'function name': function handle}
        """
        # Add basic commands
        self.command_list.update({
            # Connections
            "__.Connect": self.set_network_mode_conn,
            "__.Disconnect": self.set_network_mode_disconn,
            "__.Note_Network_Change": functionals.note_network_update,
            "__.Update_Database_Settings": functionals.update_database,
            # Execution
            "__.Run_From_DB": functionals.run_from_db,
            "ib.Confirmation": functionals.confirmation_of_action,
            "ib.Read_String": functionals.read_as_string,
            "ib.Read_File": functionals.read_as_file,
            "ib.Read_Operation": functionals.read_as_operation,
            "ib.Run_Local": functionals.read_as_operation,
            "ib.Resuscitate": functionals.resuscitate_queue_operation,
            "complete_queue": functionals.complete_queue,
            # Statuses
            "__.Resolve_Fault": functionals.resolve_fault,
            "__.Add_Fault": functionals.add_fault,
            "__.Send_State_to": functionals.send_state_to,
            "__.Read_Status_Report": functionals.read_status,
            # Misc.
            "__.Pause_for_User": functionals.wait_for_user,
            "__.Exit": self.exit,
            "__.void": functionals.nop,
            "void": functionals.db_nop,
            # Debug
            "__.Wait": functionals.wait_debug,
            "__.Echo": functionals.echo,
            "__.Except": functionals.exception_echo,
            "db_print": functionals.print_from_database,
        })
        # Load methods into here from the *_module.py files
        try:
            module = importlib.import_module(self.name + '_module')
        except ModuleNotFoundError:
            system_log.error(f"Could not build command list, could not find requisite module: '{self.name}_module'")
            return
        self.command_list.update(module.get_functions())

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    def run(self):
        """ Core loop of System

        :return: None
        """
        for func in self.core_process_functions:
            self.core_process_threads.append(Thread(target=func, daemon=True))

        for thread in self.core_process_threads:
            thread.start()

        sentinel = Thread(target=self.sentinel)
        sentinel.start()
        sentinel.join()

        system_log.info("Exiting...")

        aceso.katabasis(r"./orpheus.json", self.status)

        system_log.info("Calling CHILD_COM::__del__ to shut down")
        with self.data_pipes[CHILD_COM] as cc:
            cc.clear()
        del self.data_pipes[CHILD_COM]
        del self.client
        time.sleep(5)
        system_log.info("Shutdown signals sent")

    def sentinel(self):
        """ Manages active threads

        :return: None
        """
        counter = [0] * len(self.core_process_threads)
        thread_func_names = [x.__name__ for x in self.core_process_functions]
        while not self.kill_signal_received:
            for index, thread in enumerate(self.core_process_threads):
                if not thread.is_alive():
                    if index == 1 and self.status.get_network_state()[1] == mcs.V_OFFLINE:
                        continue
                    system_log.critical(f"Sentinel is reviving something ({index}, {thread_func_names[index]})")
                    self.core_process_threads[index] = Thread(target=self.core_process_functions[index], daemon=True)
                    self.core_process_threads[index].start()
                    counter[index] += 1
                else:
                    counter[index] = 0
            max_counter = max(counter)
            max_index = counter.index(max_counter)
            if max_counter > 3:
                if max_index == 1:
                    self.status.network_mode_is(False)
                else:
                    system_log.critical(f"Sentinel has tried {counter[max_index]} consecutive times "
                                        f"to revive process {max_index}.\n"
                                        f"Sentinel is letting the process die.")
                    break
            with self.data_pipes[CHILD_COM] as cc:
                initialized = cc['is_initialized']
            if not initialized:
                time.sleep(3)
                self.data_pipes[MSG_Q_IN].enqueue(msgm.Message.initialization_message(self.name))
                while not initialized:
                    with self.data_pipes[CHILD_COM] as cc:
                        initialized = cc['is_initialized']
                    if initialized:
                        break
                    time.sleep(1)
                system_log.info(f"{self.name} initialized!")
            time.sleep(2)
        system_log.info("Sentinel has terminated")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    def monitor(self):
        """ Whereas Sentinel manages Threads, Monitor manages system behaviors

        Right now that's just network connectivity

        :return: None
        """
        print(banner(self.name))

        while True:
            with self.data_pipes[CHILD_COM] as cc:
                initialized = cc['is_initialized']
            if not initialized:
                time.sleep(1)
            else:
                break

        if self.name == 'MC':
            monitor_ping = msgm.Message.request_of_state(requester='MC', requestee="_A")
            self.data_pipes[MSG_Q_OUT].enqueue(monitor_ping)

        while not self.kill_signal_received:
            conn, mode = self.status.get_network_state()
            if (mode == mcs.V_ONLINE) and (conn == mcs.V_ONLINE):
                if not self.data_pipes[MSG_Q_OUT].contains("ping",
                                                           _eq=lambda x, y: x == y[msgm.CONTENTS]):
                    self.data_pipes[MSG_Q_OUT].put(msgm.Message.ping())
            time.sleep(5)

        system_log.info("monitor has terminated")

    def scheduler(self):
        """ (MC Only) Tries to coordinate actions on the platform

        :return: None
        """
        my_coordinator = task_management.Coordinator(self.data_pipes, self.status)
        my_coordinator.run()
        try:
            if self.gui and self.gui.dialog and "task" in self.gui.dialog:
                self.gui.dialog["task"].on_closing()
        except (AttributeError, TypeError, KeyError):
            pass

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    def network_com(self):
        """ Manages I/O with the MCN

        :return: None
        """
        system_log.info("Starting up network communication")
        while True:
            ip_address, port_address = self.status.get_ip_and_port()
            if not isinstance(ip_address, str) or not isinstance(port_address, int):
                system_log.debug("Network mode&conn set FALSE by System.network_com() startup")
                self.status.network_is(False, False)
                time.sleep(5)
                continue
            else:
                break
        try:
            self.client.restart(ip_address, port_address)  # Runs until exception kills it (idles when mode is offline)
        except:
            system_log.exception(f"Client encountered unhandled exception")
        system_log.info("Client closed")

    def user_com(self):
        """ Manages I/O with the user

        :return: None
        """
        system_log.info("Starting up user interface")
        self.gui = module_GUI.MCN_GUI(self.data_pipes, self.status, self.command_list.keys())
        clean_shutdown = self.gui.run()
        if not clean_shutdown:
            system_log.info("GUI has died")
        self.kill_signal_received = True

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    def processor_intake(self):
        """ Spawner of threads for incoming messages

        :return: None
        """
        while not self.kill_signal_received:
            current_message = self.data_pipes[MSG_Q_IN].get_next()

            if isinstance(current_message, msgm.Message) and current_message[msgm.CONTENTS] == "ping":
                continue

            if self.debug_mode:
                options = ["Run", "Destroy", "Skip", "Skip & Demote", "Exit Debug Mode"]
                reference = pformat(current_message.print(False, True), width=40)
                user_response = functionals.quick_select(title="Debug Mode",
                                                         dialog=f"Next message:\n"
                                                                f"{reference}\n\n"
                                                                f"--Select Action--",
                                                         options=options,
                                                         default="Run")
                if user_response == "Run":
                    if reference in self.debug_memory:
                        current_message[msgm.PRIORITY] = self.debug_memory[reference]
                        del self.debug_memory[reference]
                elif user_response == "Destroy":
                    continue
                elif user_response == "Skip":
                    self.data_pipes[MSG_Q_IN].enqueue(current_message)
                elif user_response == "Skip & Demote":
                    self.debug_memory[reference] = current_message[msgm.PRIORITY]
                    current_message[msgm.PRIORITY] = 8
                    self.data_pipes[MSG_Q_IN].enqueue(current_message)
                elif user_response == "Exit Debug Mode":
                    self.debug_mode = False
                else:
                    raise RuntimeError(f"In Debug Mode, user response in unknown state: "
                                       f"{user_response}.\nNeeds be one of: {options}")

            system_log.info(f"Processing: {current_message.prints()}")

            if not (self.name == current_message[msgm.RECIP]):
                # This wasn't for this system, put it back on the message queue and the communicator will handle it
                self.data_pipes[MSG_Q_OUT].enqueue(current_message)
            else:
                # Message is for this system
                if self.validate_key(current_message[msgm.IDKEY]):
                    child = Thread(target=functionals.execute_message,
                                   args=(current_message, self.data_pipes),
                                   kwargs={'cmd_list': self.command_list, 'status': self.status},
                                   daemon=True)
                    with self.child_processes as cp:
                        cp.append(child)
                    child.start()
                else:
                    pass

    def processor_cleanup(self):
        """ Garbage collector for spent threads

        :return: None
        """
        while not self.kill_signal_received:
            with self.child_processes as cp:
                for index, child in enumerate(cp):
                    if not child.is_alive:
                        del cp[index]
            time.sleep(5*60)  # every 5 minutes

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    def validate_key(self, idemp_key):
        """ Checks the idempotency key of a message against those already received

        :param idemp_key: idempotency key of incoming message
        :return: True/False - Message is new, run message
        """
        if idemp_key == msgm.ID_KEY_BYPASS:
            return True
        if idemp_key in self.received_keys:
            system_log.debug(f"Rejecting key: {idemp_key}")
            return False
        else:
            self.received_keys.callattr('append', idemp_key)
            return True

    def clean_keyring(self):
        """ Routinely removes old idempotency keys

        :return: None
        """
        while True:
            expired_key_date = (date.today() - timedelta(days=5)).strftime("%Y%m%d")
            # I have no clue if this with statement will actually work (python-console analogous test worked so...)
            with self.received_keys as keychain:
                for index, idkey in enumerate(keychain):
                    if idkey[2:10] == expired_key_date:
                        del keychain[index]
            time.sleep(12*60*60)

    def set_network_mode_conn(self, *_, **__):
        """ Command for setting the network mode to connected (that the System *should* be connected)

        Calls Status.network_mode_is(True) and MCClient.connect()

        :return: A Return Object with data populated by the return from MCClient.connect()
        """
        self.status.network_mode_is(True)
        return mcs.RetObj.complete(self.client.connect())

    def set_network_mode_disconn(self, *_, **__):
        """ Command for setting the network mode to disconnected (that the System *should not* be connected)

        Calls MCClient.disconnect()

        :return: A complete Return Object
        """
        self.client.disconnect()
        return mcs.RetObj.complete()

    def exit(self):
        """ Sets the Kill Signal flag to True

        :return: A complete Return Object
        """
        self.kill_signal_received = True
        return mcs.RetObj.complete("Kill signal received")

# ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## ######## #
