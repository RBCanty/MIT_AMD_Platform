""" mc_server.py

A server network node

Use of "..." may require Python v3.10 or later (could replace with NARG from custom_classes.py if this causes problems)

@author: Ben Canty
"""

import json
import socket
import time
import tkinter as tk
from datetime import datetime
from queue import Queue, Empty
from threading import Thread, Event
from traceback import print_exc
from typing import Tuple, Dict

import mcn_status
import message as msgm
import operations as oprtn
from custom_classes import safe_open, pathfinder
from mcn_status import MCN_CFG
from thread_safe_data import ThreadSafeDataContainer
from gui_constants import T33


class ServerLog:
    """ Used to log actions on the server independently of the built-in logging module

    Use of the built-in logging module resulted in writes permissions errors (the server and MC run on the same
    computer)--the solution is probably tied to implementing a hierarchical logging system.
    """
    def __init__(self):
        """ Constructs a logger for the server

        Messages are of the form::

        "yyyy-mm-dd hh:mm:ss - level - module - message"
        """
        self.format = "%Y-%m-%d %H:%M:%S"

    @property
    def logfile(self):
        """ Locates and safe-opens the logfile """
        filename = pathfinder(r"Logs", f"server_logfile ({datetime.now().strftime('%Y-%m-%d')}).log")
        return safe_open(filename, 'a+')

    @property
    def date(self):
        """ Generates the current timestamp in the logger's format """
        return datetime.now().strftime(self.format)

    def debug(self, module, message):
        """ Commits a debug-level message

        :param module: The module from which the call is made
        :param message: The message to be committed
        :return: None
        """
        msg = f"{self.date} - DEBUG - {module} - {message}"
        print(msg)
        self.log_to_file(msg)

    def info(self, module, message):
        """ Commits an info-level message

        :param module: The module from which the call is made
        :param message: The message to be committed
        :return: None
        """
        msg = f"{self.date} - INFO - {module} - {message}"
        print(msg)
        self.log_to_file(msg)

    def error(self, module, message):
        """ Commits an error-level message

        :param module: The module from which the call is made
        :param message: The message to be committed
        :return: None
        """
        msg = f"{self.date} - ERROR - {module} - {message}"
        print(msg)
        self.log_to_file(msg)

    def exception(self, module, message):
        """ Commits a warning-level message with exception traceback

        :param module: The module from which the call is made
        :param message: The message to be committed
        :return: None
        """
        msg = f"{self.date} - WARN - {module} - {message}\n"
        print(msg)
        print_exc()
        self.log_to_file(msg, True)  # This will take care of printing the traceback to the log-file

    def warning(self, module, message):
        """ Commits a warning-level message

        :param module: The module from which the call is made
        :param message: The message to be committed
        :return: None
        """
        msg = f"{self.date} - WARN - {module} - {message}"
        print(msg)
        self.log_to_file(msg)

    def log_to_file(self, msg, is_exc=False):
        """ Writes log messages to the log file

        :param msg: The message to be written to the log file
        :param is_exc: False (default) - do not look for traceback, True - write traceback to file as well
        :return: None
        """
        try:
            with self.logfile as log:
                log.write(msg + "\n")
                if is_exc:
                    print_exc(file=log)
                    log.write("\n")
        except (FileNotFoundError, FileExistsError, PermissionError) as fe:
            print(f"Log file could not be written to: {repr(fe)}")


class ServerMessageBuffer(Queue):
    """ A queue for messages in the Server which have backed up """
    def put(self, item: msgm.Message, block=..., timeout=...):
        """ Adds a message to the queue

        :param item: A Message
        :param block: Block argument to Queue constructor
        :param timeout: Timeout argument to Queue constructor
        :return: None
        """
        sender = item[msgm.SENDER]
        recip = item[msgm.RECIP]
        contents = item[msgm.CONTENTS]

        # Only add a ping if one does not already exist in the queue
        if contents == "ping":
            with self.mutex:
                for msg in self.queue:
                    if recip != msg[msgm.RECIP]:
                        continue
                    if msg[msgm.CONTENTS] == "ping":
                        return

        # Safety checks for ...
        if contents == "ib.Read_Operation":
            operation = oprtn.Operation.build_from_json(json.loads(item.read_data()))
            function = operation[oprtn.FUNC]
            # ... if the message is a Note of the network members changing, only keep the most recent version
            if function == "__.Note_Network_Change":
                with self.mutex:
                    for msg in self.queue:
                        if recip != msg[msgm.RECIP]:
                            continue
                        _function = oprtn.Operation.build_from_json(json.loads(msg.read_data()))[oprtn.FUNC]
                        if _function == "__.Note_Network_Change":
                            msg[msgm.DATA] = item[msgm.DATA]
                            return
            # ... if the message is an update of current state, only keep the oldest(?) version
            elif function == "__.Send_State_to":
                with self.mutex:
                    for msg in self.queue:
                        if (recip != msg[msgm.RECIP]) or (sender != msg[msgm.SENDER]):
                            continue
                        _function = oprtn.Operation.build_from_json(json.loads(msg.read_data()))[oprtn.FUNC]
                        if _function == "__.Send_State_to":
                            # Should there be a line like "msg[msgm.DATA] = item[msgm.DATA]" here?
                            return

        # Add message to queue
        super(ServerMessageBuffer, self).put(item, block, timeout)


class ServerUI:
    """ A GUI for the Server """
    def __init__(self, root: tk.Tk):
        """ Creates a popup for the server, such that it can be controlled without an Agent connected

        :param root: A tkinter root for the popup
        """
        self.root = root
        self.root.protocol("WM_DELETE_WINDOW", self.kill)
        self.flag: Event = None  # noqa
        self.buffer: ServerMessageBuffer = None  # noqa

        master = tk.Frame(root)
        master.winfo_toplevel().title("Server")
        master.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        tk.Label(master, text="MCN Server\nMember Systems:").pack(**T33)

        self.sys_list = tk.Text(master, height=1, width=24)
        self.sys_list.pack(**T33)

        self.sel_client = tk.StringVar(master, "--")
        ban_window = tk.LabelFrame(master, text="Kick from server")
        ban_window.pack(fill=tk.BOTH, **T33)
        self.option_menu = tk.OptionMenu(ban_window, self.sel_client, *["--", ])
        self.option_menu.pack(**T33)
        self.option_menu.bind('<Button-1>', self.update_opt_menu)
        tk.Button(ban_window, text="Kick", command=self.kick).pack(**T33)

        tk.Button(master, text="Shut Down", command=self.kill).pack(side=tk.TOP, ipadx=3, ipady=3)
        self.kill_switch = False

        self.root.after(500, self.check)

    def attach(self, flag: Event, buffer: ServerMessageBuffer):
        """ Connects the GUI to the server by linking the flag and buffer between the objects

        :param flag: An event which signals the server to close down
        :param buffer: the message buffer of the server
        :return: None
        """
        self.flag = flag
        self.buffer = buffer

    def check(self):
        """ Checks the kill flag periodically

        :return: None
        """
        if self.flag and self.flag.is_set():
            self.kill()
        self.root.after(500, self.check)

    def kill(self):
        """  Destroys the window (tkinter destroy) and the server (by setting the kill flag)

        :return: None
        """
        try:
            if self.flag:
                self.flag.set()
        finally:
            try:
                self.root.destroy()
            except:  # noqa
                pass

    def update_sys_list(self, systems: list):
        """ Updates the list of systems online

        :param systems: An Iterable of system names
        :return: None
        """
        self.sys_list['state'] = tk.NORMAL
        self.sys_list.delete("1.0", tk.END)
        self.sys_list.insert(tk.INSERT, ", ".join(systems))
        self.sys_list['state'] = tk.DISABLED

    def update_opt_menu(self, *_):
        """ Updates the option menu for booting members from the server

        :param _: Consumes arguments passed by tkinter
        :return: None
        """
        sys_list = self.sys_list.get("1.0", tk.END).replace("\n", "")
        sys_list = ["--", *[v for v in sys_list.split(", ") if v]]
        self.option_menu["menu"].delete(0, tk.END)  # update_opt_menu
        for _sys in sys_list:
            self.option_menu["menu"].add_command(label=_sys, command=lambda v=_sys: self.sel_client.set(v))

    def kick(self):
        """ Removes a member from the network

        :return: None
        """
        sel_client = self.sel_client.get()
        if self.buffer and sel_client and (sel_client != "--"):
            self.buffer.put(msgm.Message.leave_server(sel_client))
            print(f"Kicking '{sel_client}'")
        else:
            print("Ignored")


# Create the server log
server_log = ServerLog()


class MCServer:
    """ The server for the MCN """
    def __init__(self, ip_address=MCN_CFG['network']['IP'], port_address=MCN_CFG['network']['PORT'], ui=None):
        """ Creates a Server

        :param ip_address: IP address of the server as a string
        :param port_address: Port number of the server as an integer
        :param ui: A ServerUI object, if desired
        """

        self.name = "_S"
        self.ip_address = ip_address
        self.port_address = port_address

        self.buffer = ServerMessageBuffer()
        self.flag = Event()

        self._build_server_socket()

        self._workers: Dict[str, Tuple[Thread, socket.socket]] = dict()
        self.workers = ThreadSafeDataContainer(self._workers)
        self.server_thread = Thread(target=self.server_worker, daemon=True, name="S_WORKER")
        self.main_thread = Thread(target=self.processor, daemon=True, name="S_PROCESS")

        self.ui: ServerUI = ui
        if self.ui:
            self.ui.attach(self.flag, self.buffer)

    def __del__(self):
        """ Soft destructor for the GUI

        Python does not reliably call __del__ (use a context manager next time)

        :return: None
        """
        if self.ui:
            self.ui.kill()

    @property
    def kill_switch(self):
        """ A boolean for if the kill flag is set """
        return self.flag.is_set()

    def _build_server_socket(self):
        """ Creates, binds, and initiates listening to a socket

        :return:
        """
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        self.sock.settimeout(1)
        self.sock.bind((self.ip_address, self.port_address))
        self.sock.listen()

    def _add_worker(self, name, _socket):
        """ Adds a worker for a given system to handle IO

        :param name: Name of the Agent
        :param _socket: Socket of the Agent
        :return: False - no addition, None - otherwise
        """
        with self.workers as workers:
            if (name in workers) and (workers[name][0].is_alive()):
                server_log.warning("_add_worker", f"Could not add worker {name}, already registered")
                return False
            workers[name] = (
                Thread(target=self.read_worker, args=(name, _socket, ), daemon=True, name=f"{name}_WORKER"),
                _socket
            )
            workers[name][0].start()
        self.print_member_list()

    def _remove_worker(self, name):
        """ Removes a worker

        If worker is the server, triggers server shutdown

        :param name:
        :return:
        """
        if name == "_S":
            self.flag.set()
        with self.workers as workers:
            workers.pop(name, None)
        self.print_member_list()

    @staticmethod
    def listen(name, sock: socket.socket):
        """ Watches the socket for a binary header, decodes the message, and returns the Message

        :param name: The name associated with the socket being monitored
        :param sock: The socket being monitored
        :return: A message object or None if failed
        """

        # await a header
        binary_header = sock.recv(64)
        if not binary_header:
            return None

        # Decode the header
        _header = msgm.Message.parse_header(binary_header)

        # Clear socket if header is bad
        if not (len(_header) == 4):
            server_log.error("listen", f"Server received a bad message from {name}, attempting to clear socket")
            accumulator = 0
            try:
                bad_data = sock.recv(4096)
                accumulator += len(bad_data)
                while len(bad_data) == 4096:
                    bad_data = sock.recv(4096)
                    accumulator += len(bad_data)
            except KeyboardInterrupt as ki:  # allow keyboard interrupt as a way to close the server if GUI absent
                raise ki
            except Exception:  # noqa
                server_log.exception("listen", f"Server failed to clear socket ({name}) of bad data")
            else:
                server_log.info("listen", f"Server flushed {accumulator} bytes to clear socket ({name})")
            finally:
                return None

        # Read in message
        try:
            binary_body = sock.recv(int(_header[3]))
        except ConnectionError:
            server_log.exception("listen", f"Server encountered a connection error ({name})")
            return None

        # Translate message
        my_message = msgm.Message.load_message(binary_body)
        return my_message

    def monitor(self):
        """ Watches all worker's sockets

        :return: None
        :raises KeyboardInterrupt: On KeyboardInterrupt, otherwise exceptions mark a worker to be removed
        """
        prev_member_list = []
        while not self.kill_switch:
            try:
                members = sorted(list(self.workers.callattr('keys')))
                if members == prev_member_list:
                    hit_list = []
                    with self.workers as workers:
                        for _r, (_, _s) in workers.items():  # name, (thread, socket)
                            temp_sock: socket.socket = _s
                            try:
                                temp_sock.sendall(msgm.Message.ping(_r).encode())
                            except KeyboardInterrupt as ki:
                                raise ki
                            except:  # noqa
                                hit_list.append(_r)
                    for hit in hit_list:
                        self._remove_worker(hit)
                    time.sleep(2.5)
                else:
                    print("DEBUG MEMBER LIST CHANGED")
                    self.buffer.put(msgm.Message.note_network_change(members))
                    prev_member_list = members
            except KeyboardInterrupt as ki:
                raise ki
            except:  # noqa
                pass
            finally:
                time.sleep(0.5)

    def read_worker(self, name, _socket: socket.socket):
        """ Worker which listens to incoming messages

        :param name: The associated system
        :param _socket: The associated socket
        :return: None
        """
        server_log.info("read_worker", f"New reade_worker created for {name}:{_socket}")
        _socket.settimeout(1)
        while not self.kill_switch:
            try:
                new_message = self.listen(name, _socket)
            except socket.timeout:
                continue
            except OSError:
                return
            if new_message:
                if new_message[msgm.CONTENTS] != "ping":
                    server_log.info("read_worker", f"Received new message from {name}: {new_message.prints()}")
                self.buffer.put(new_message)
            time.sleep(0.1)

    def server_worker(self):
        """ Worker which reads messages for the server itself

        :return: None
        """
        while not self.kill_switch:
            try:
                client_sock, _ = self.sock.accept()
            except socket.timeout:
                continue
            except OSError:
                if not self.kill_switch:
                    server_log.exception("server_worker", f"self.sock is {type(self.sock)} not 'socket'")
                return
            client_sock.setblocking(True)
            server_log.info("server_worker", "New Socket accepted!")
            new_message = self.listen('new', client_sock)
            if not new_message:
                return
            server_log.info("server_worker", f"Received new message: {new_message.prints()}")

            if new_message[msgm.CONTENTS] == "connect_signal":
                server_log.debug("server_worker", "SERVER RECEIVED CONNECT SIGNAL")
                user = new_message[msgm.SENDER]
                self._add_worker(user, client_sock)
                # Given that something just connected, we should tell it to report its status to the MC...
                if user not in ["MC", "", "_S"]:
                    startup_ping = msgm.Message.request_of_state(requester='MC', requestee=user)
                    self.buffer.put(startup_ping)
                # ...unless it was the MC, in which case, all other systems should greet the MC with their state
                elif user == "MC":
                    resync_ping = msgm.Message.request_of_state(requester='MC', requestee="_A")
                    self.buffer.put(resync_ping)
            elif new_message[msgm.CONTENTS] == "disconnect_signal":
                server_log.debug("server_worker", "SERVER RECEIVED *DIS*CONNECT SIGNAL")
                self.buffer.put(new_message)
            else:  # A message to the server that isn't a (dis)connect signal is just a ping
                png_response = {msgm.TYPE: 'COM',
                                msgm.SENDER: "_S",
                                msgm.RECIP: new_message[msgm.SENDER],
                                msgm.IDKEY: "",
                                msgm.CONTENTS: new_message[msgm.CONTENTS], msgm.DATA: None,
                                msgm.PRIORITY: msgm.P_PNG}
                self.buffer.put(msgm.Message(png_response))

    def processor(self):
        """ Translate messages from the buffer into action

        :return: None
        """
        while not self.kill_switch:
            # Pull a message
            try:
                new_message: msgm.Message = self.buffer.get(timeout=1)
            except Empty:
                continue
            if not new_message:
                continue
            if new_message[msgm.CONTENTS] == "ping":
                continue
            # Collect basic info (who sent it and to whom it is sent--both system and subsystem)
            recipient_sub = new_message[msgm.RECIP]
            recipient_maj = mcn_status.rs_dir(recipient_sub)
            sender_sub = new_message[msgm.SENDER]
            sender_maj = mcn_status.rs_dir(sender_sub)

            # Server Stuff
            if new_message[msgm.CONTENTS] == "connect_signal":
                user = new_message[msgm.SENDER]
                server_log.debug("processor", f"PROCESSOR RECEIVED CONNECT SIGNAL from {user}")
                continue
            if new_message[msgm.CONTENTS] == "disconnect_signal":
                user = new_message[msgm.SENDER]
                server_log.debug("processor", f"PROCESSOR RECEIVED *DIS*CONNECT SIGNAL from {user}")
                self._remove_worker(user)
                continue
            if new_message[msgm.CONTENTS] == "close_server":
                server_log.info("processor", "Server has been commanded to shut down")
                with self.workers as workers:
                    for _r, (_, _s) in workers.items():
                        temp_sock: socket.socket = _s
                        try:
                            temp_sock.sendall(msgm.Message.leave_server(_r).encode())
                        except KeyboardInterrupt as ki:
                            raise ki
                        except:  # noqa
                            pass
                self.flag.set()
                continue

            # Client Stuff
            if recipient_maj == "_A":
                server_log.info("processor", "Reprocessing @All Message")
                mailing_list = [k for k in mcn_status.MCN_CFG[mcn_status.S_MAJ_SYS]
                                if k not in [sender_maj, "_S", "_U", "_A"]]
                for recip in mailing_list:
                    self.buffer.put(new_message.modified_copy({msgm.RECIP: recip}))
                continue
            if recipient_maj not in self.workers:
                if recipient_maj in mcn_status.MCN_CFG[mcn_status.S_MAJ_SYS]:
                    self.buffer.put(new_message)
                continue
            # Misc. Stuff
            if recipient_maj in ["__", "_U", ""]:
                # burn the message
                server_log.info("processor (void)", new_message.print(sysio=False))
                continue

            # Message Forwarding
            bad_flag = False
            with self.workers as workers:
                server_log.info("processor",
                                f"Server is attempting to forward a message '{new_message.prints()}'...")
                try:
                    client_socket: socket.socket = workers[recipient_maj][1]
                    client_socket.sendall(new_message.encode())
                    server_log.info("processor",
                                    f"... message forwarded: '{new_message.print(sysio=False)}' ")
                except KeyError:
                    server_log.exception("processor",
                                         f"...forward of message to '{recipient_maj}' failed, "
                                         f"{recipient_maj} not connected")
                    self.buffer.put(new_message)
                except (ConnectionError, ConnectionResetError, ConnectionAbortedError, ConnectionRefusedError, OSError):
                    server_log.exception("processor",
                                         f"...forward of message to '{recipient_maj}' failed with exception:")
                    self.buffer.put(new_message)
                    bad_flag = True
            if bad_flag:
                self._remove_worker(recipient_maj)

    def host_server(self):
        """ Perpetuates threads for workers (this is the __main__ of the server)

        :return: None
        :raises KeyboardInterrupt: On KeyboardInterrupt
        """
        self.server_thread.start()
        self.main_thread.start()
        Thread(target=self.monitor, daemon=True, name="MONITOR").start()

        while not self.kill_switch:
            if not self.server_thread.is_alive():
                self._build_server_socket()
                self.server_thread = Thread(target=self.server_worker, daemon=True, name="S_WORKER")
                self.server_thread.start()
            time.sleep(0.01)

        try:
            server_log.info("host_server", "Attempting to close all sockets gracefully")
            for _, (_, sock) in self._workers.items():
                try:
                    sock.close()
                except KeyboardInterrupt as ki:
                    raise ki
                except:  # noqa
                    pass
            if self.buffer.qsize() > 0:
                server_log.warning("host_server", f"Server closed with {self.buffer.qsize()} items in its buffer")
                while self.buffer.qsize() > 0:
                    try:
                        item: msgm.Message = self.buffer.get_nowait()
                        server_log.info("host_server", f"Item {self.buffer.qsize()}: {item.prints()}")
                    except KeyboardInterrupt as ki:
                        raise ki
                    except:  # noqa
                        pass
            else:
                server_log.info("host_server", "Server closed with no items in its buffer")
        except KeyboardInterrupt as ki:
            raise ki
        except:  # noqa
            server_log.exception("host_server", "Exception when closing out Server")
        else:
            server_log.info("host_server", "Closing out Server")

    def print_member_list(self):
        """ Reports to the logger the current members of the network

        :return: None
        """
        # create a str and a list representation of the member list (unfortunate nomenclature)
        member_list_str = ""
        member_list_list = []
        for key, (_, val) in self._workers.items():
            member_list_str += f"{' '*47}{key} {val}\n"  # formatting to look nice with the log headers
            member_list_list.append(key)

        server_log.info("print_member_list", f"MCN-Members:\n{member_list_str}")
        if self.ui:
            self.ui.update_sys_list(member_list_list)


if __name__ == '__main__':
    tk_root = tk.Tk()
    my_server = MCServer(ui=ServerUI(tk_root))
    Thread(target=my_server.host_server, name="__SERVERHOST__").start()
    tk_root.mainloop()
