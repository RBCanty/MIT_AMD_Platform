""" mc_client.py

A client network node

Historical note: This and the server were the first parts of the entre project to be coded, which is why they are at
the socket level rather than using an established module such as requests.  It should have been updated to such a model
during the v0 --> v1 update, but that update was focused on operating from a database (rather than GUI). And while the
network component was updated at that time to handle network interruptions, overhauling it to use existing Python
modules wasn't a priority at the time.

@author: Ben Canty
"""

import socket
import time
from threading import Thread

import mcn_status as mcs
import message
from constants import *
from mcn_logging_manager import system_log
from mcn_queues import MCNPriorityQueue
import message as msgm


class MCClient:
    """ System-side network node for MCN communication """
    def __init__(self, name, pipes: dict, ip_address, port_address, ):
        """ Creates a client node on the MCN

        :param name: The name of the system (e.g. "LC", "AH", etc.)
        :param pipes: A dictionary of data pipes for message passing with the System
        :param ip_address: the IP address of the server
        :param port_address: the Port address of the server
        """
        self.name = name

        self.ip_address = str(ip_address)
        self.port_address = int(port_address)
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.sock.settimeout(1)

        self.status: mcs.Status = pipes[INTERNAL]
        self.inbox: MCNPriorityQueue = pipes[MSG_Q_IN]
        self.outbox: MCNPriorityQueue = pipes[MSG_Q_OUT]

    def __del__(self):
        """ Weak close-out procedure

        Note: Python does not reliably call __del__

        :return: None
        """
        try:
            self.sock.close()
        except:  # noqa
            pass

    def connect(self):
        """ Connects to the network

        Sets 'mode' to True and updates 'conn' to reflect the success

        :return: Either an Online or Offline code
        """
        self.status.network_mode_is(True)  # This means "system *should* be online", but the state for "*is* online"
        # is a different attribute (see "status.network_conn_is()" below)

        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.sock.settimeout(1)
        try:
            self.sock.connect((self.ip_address, self.port_address))
            self.sock.send(message.Message.connect(self.name).encode())
            self.status.network_conn_is(True)  # connection successful
            return mcs.V_ONLINE
        except (socket.timeout, OSError, ConnectionError):
            system_log.exception("Network conn set FALSE by Client.connect()")
            self.status.network_conn_is(False)  # connection failed
            return mcs.V_OFFLINE

    def disconnect(self):
        """ Disconnects from the network

        Sets 'mode' and 'conn' to False

        :return: None
        """
        conn, _ = self.status.get_network_state()
        self.status.network_mode_is(False)
        if conn == mcs.V_ONLINE:
            try:
                self.sock.sendall(message.Message.disconnect(self.name).encode())
                self.sock.close()
            except:  # noqa
                system_log.exception(f"From {self.name} Client, disconnect was not clean")
            else:
                system_log.info(f"From {self.name} Client, clean disconnect")
        elif conn == mcs.V_OFFLINE:
            try:
                self.sock.close()
            except:  # noqa
                pass
            system_log.info(f"From {self.name} Client, confirmed disconnect")
        self.status.network_conn_is(False)

    def send_msg(self, binary_msg):
        """ Calls socket.sendall() from the Communicator class rather than the child (Client/Server) class

        :param binary_msg: A binary object which can be sent over a socket
        :return: None
        :raises ConnectionError: On OSError from socket.sendall()
        """
        try:
            self.sock.sendall(binary_msg)
        except OSError as e:
            raise ConnectionError(e)

    def flush_socket(self):
        """ Cleans socket of content

        :return: None
        """
        system_log.error(f"{self.name} received a bad message, attempting to cleat socket")
        accumulator = 0
        try:
            while True:
                bad_data = len(self.sock.recv(4096))
                accumulator += bad_data
                if bad_data < 4096:
                    break
        except socket.timeout:
            pass
        except Exception:  # noqa
            system_log.exception(f"{self.name} failed to clear socket of bad data")
        system_log.info(f"{self.name} flushed {accumulator} bytes to clear socket")

    def listen(self):
        """ Watches the socket for a binary header, decodes the message, and returns the Message

        :return: A message object or None if failed
        :raises ConnectionError: On a socket.recv() exception (other than timeout, which returns None)
        """
        # Await a message header
        try:
            binary_header = self.sock.recv(64)
        # If timeout, return None (and wait to be called again)
        except socket.timeout:
            return None
        # If error, raise exception
        except Exception as e:
            system_log.exception(f"{self.name} encountered a connection error awaiting a header")
            raise ConnectionError from e
        # If the header is somehow null
        if not binary_header:
            return None

        # Convert header into readable form [Type, sender, idempotency key, length of message]
        _header = msgm.Message.parse_header(binary_header)
        if not (len(_header) == 4):
            self.flush_socket()
            return None

        # Lookup message length and read in that many bytes
        try:
            binary_body = self.sock.recv(int(_header[3]))
        except socket.timeout:
            system_log.exception(f"{self.name} encountered a timeout, message header gave wrong message length")
            return None
        except ConnectionError as ce:
            system_log.exception(f"{self.name} encountered a connection error reading message contents")
            raise ce

        # Convert message to a Message object and hand off to caller
        return msgm.Message.load_message(binary_body)

    def run_in_loop(self, method, exec_interval=0.1, pause_interval=3):
        """ Wrapper to execute a method perpetually but to pause when the mode or conn are offline

        The loop is only broken when an exception is raised to this level

        :param method: A method to be executed when the client is online and connected
        :param exec_interval: The interval on which the loop is executed
        :param pause_interval: The period to wait if either connection or mode is False
        :return: None
        :raises Exception: (From the given method to halt the loop)
        """
        while True:
            time.sleep(exec_interval)
            conn, mode = self.status.get_network_state()
            if (mode == mcs.V_OFFLINE) or (conn == mcs.V_OFFLINE):
                time.sleep(pause_interval)
                continue
            method()

    def writer(self):
        """ Sends messages to server (executes once)

        Catches ConnectionError, other exceptions will be raised to caller

        :return: None
        """
        msg = self.outbox.get_next()  # Is blocking, so loop does not need a time.sleep() call
        if (not msg) or (not isinstance(msg, msgm.Message)):
            return
        try:
            self.send_msg(msg.encode())
        except ConnectionError:
            self.outbox.enqueue(msg)
            system_log.exception("Network conn set FALSE by Client.writer() via exception")
            self.status.network_conn_is(False)

    def reader(self):
        """ Receives messages from server (executes once)

        Catches ConnectionError, other exceptions will be raised to caller

        :return: None
        """
        try:
            msg = self.listen()
        except ConnectionError:
            system_log.exception("Network conn set FALSE by Client.reader() via exception")
            self.status.network_conn_is(False)
            return
        if (not msg) or (not isinstance(msg, msgm.Message)):
            return
        self.inbox.enqueue(msg)

    def run(self):
        """ Run method for activating the Client

        :return: None upon reader or writer dying
        """
        reader_thread = Thread(target=self.run_in_loop, args=(self.reader, ), daemon=True)
        writer_thread = Thread(target=self.run_in_loop, args=(self.writer, ), daemon=True)
        threads = [writer_thread, reader_thread]

        for thread in threads:
            thread.start()

        while True:
            # If a thread has died, restart the client
            if not reader_thread.is_alive() or not writer_thread.is_alive():
                self.status.network_conn_is(False)
                return

            # Check the network state (reconnect/disconnect as needed)
            conn, mode = self.status.get_network_state()
            if (conn == mcs.V_OFFLINE) and (mode == mcs.V_ONLINE):
                system_log.info("Trying to reconnect to server")
                self.connect()
            elif (conn == mcs.V_ONLINE) and (mode == mcs.V_OFFLINE):
                system_log.info("Trying to disconnect from server")
                self.disconnect()

            # Sleep, then check again
            time.sleep(5)

    def restart(self, ip=None, port=None):
        """ Updates the IP and Port addresses, and invokes self.run() if self.connect() is successful

        :param ip: The new IP to use
        :param port: The new port to use
        :return:
        """
        if ip and port:
            self.ip_address = str(ip)
            self.port_address = int(port)
        if self.connect() == mcs.V_ONLINE:
            self.run()
