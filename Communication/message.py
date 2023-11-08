""" Message class

Proxies Python Dictionary

Basic unit of network communication

Class Fields
  - TYPE      : 'COM' for networking commands, 'CMD' for platform commands
  - SENDER    : System who sent the message
  - RECIP     : Intended recipient of the message ('_A' for 'to all')
  - IDKEY     : Idempotency key, used to prevent repeat messages and to track completion of jobs triggered by a message
  - CONTENTS  : A network command (COM) or a command string instructing how DATA is to be interpreted (CMD)
  - DATA      : String representation of whatever data type is specified by CONTENTS
  - PRIORITY  : An integer rank (0-8; 0-highest) which controls priority and if a response
                is needed or not for CMD Messages

Class Methods
  - __setitem__()  - Allows constrained mutation of Message's values with bracket indexing, like a Python Dictionary
  - __getitem__()  - Allows access of Message's values with [] indexing
  - encode()       - Converts a Message into a Bytes object for socket traversal
  - add_data()     - kwarg=data <- a string representation of an object, or kwarg=filename <- Full path of file
                     being encoded into the message (uses utf-8 encoding of base64 literal encoding to ensure
                     traversal over string objects is safe)
  - read_data()    - If a kwarg=filename is provided, it writes the string in DATA to the specified file,
                     otherwise, it returns DATA as a string
  - prints()       - Returns a string summary of the Message: "TYPE (SENDER>RECIP) IDKEY[16:]: 'CONTENTS'"
  - print()        - Prints a string representation (IDKEY and first 100 characters of CONTENTS as a string) of
                     the Message to the console or returns a list(str) version to the caller

Static Methods
  - build_from_args()        - Allows a Message to be constructed using kwargs rather than supplying a dictionary
  - parse_header()           - Takes the Binary header associated with an encoded Message and returns a list of
                               [TYPE, RECIP, IDKEY, Length of encoded Message]
  - load_message()           - Takes the encoded Message associated with a Binary Header and returns a Message object
  - gen_idemp_key()          - Generates an idempotency key (2:sender name)(14:timestamp)(16:unique characters)

Quick Message Maker
  -   connect()                - Creates a COM message to connect to the server
  -   disconnect()             - Creates a COM message to disconnect from the server
  -   close_server()           - Creates a COM message which will make the server shut down
  -   request_of_state()       - Creates a CMD message asks one system to send its state to another
  -   status_report_message()  - Creates a CMD message which contains a summary of a system's state
  -   note_network_change()    - Creates a CMD message for updating the MC that a system has disconnected form the MCN
  -   confirm_ro()             - Creates a CMD message which is sent to the issuer of a previous message that
                                 said message has been completed
  -   kill()                   - TODO: Creates a CMD message which is sent to have the recipient shutdown its subsystems

@author: Ben Canty
"""

import base64  # used for File Transfer, but isn't used in practice given databases
import random
from copy import deepcopy
from datetime import datetime

import operations as oprtn
from mcn_status import rs_dir
import json
from args_and_kwargs import *
from mcn_logging_manager import system_log

MSG_TYPES = ['CMD', 'COM']
TYPE = "type"
SENDER = "sender"
RECIP = "recipient"
IDKEY = "idempotency-key"
CONTENTS = "contents"
DATA = "data"
PRIORITY = "priority"
MESSAGE_KEYS = [TYPE, SENDER, RECIP, IDKEY, CONTENTS, DATA, PRIORITY]
FORMAT = 'Ascii'
P_EMERGENCY = 0
P_ERROR = 1
P_WARN = 2
P_SERVER = 3
P_CMD_RR = 4
""" Command that requests a confirmation of receipt and execution """
P_CMD_NR = 5
""" Command that does not request a confirmation of receipt and execution """
P_RSP = 6
P_7 = 7
P_PNG = 8
ID_KEY_SOURCE = '#0123456789abcdefghijklmnopqrstuvwxyz-ABCDEFGHIJKLMNOPQRSTUVWXYZ+'
ID_KEY_BYPASS = "0" * 32  # python: a string of 32 zero-characters, i.e. "000...0"


class Message(object):
    """ A container class for the creation and formatting of messages within the MCN

    A Message is just a dictionary object with constrained keys and values
    Type - CMD (issuing a command) or COM (networking)
    Sender - The name of the message's sender
    Recipient - The name of the message's recipient (by System)
    Idempotency Key - (2:sender name)(14:timestamp)(16:unique characters)
    Contents - (CMD) A string which specifies how Message[DATA] is to be interpreted, (COM) a network command string
    Data - A message, file, or object
    Priority - A number: 0 (Highest) - 8 (Lowest)
    """

    def __init__(self, dictionary):
        """ Creates a message object

        :param dictionary: A dictionary containing the requisite fields for a Message
                           (extra keys discarded, and missing keys given default values)
        """
        self._dictionary = dictionary
        for k in self._dictionary.keys():
            if k not in MESSAGE_KEYS:
                del self._dictionary[k]
        for k in MESSAGE_KEYS:
            self[k] = self[k] if k in self._dictionary.keys() else None

    @staticmethod
    def build_from_args(_type=None, _sender=None, _recip=None, _idkey="", _contents=None, _data=None, _priority=8):
        """ Allows a message to be built through arguments rather than a dictionary

        :param _type: CMD (issuing a command) or COM (networking)
        :param _sender: The name of the message's sender
        :param _recip: The name of the message's recipient (by System)
        :param _idkey: A string, "" auto-generates, "bypass" gives a special value
        :param _contents: (CMD) A string which specifies how Message[DATA] is to be interpreted,
        (COM) a network command string
        :param _data: A message, file, or object encoded according to _contents
        :param _priority: (See Information.txt; a number 0-8)
        :return: A Message object
        """
        self = Message(dict())
        self[TYPE] = _type
        self[SENDER] = _sender
        self[RECIP] = _recip
        self[IDKEY] = _idkey
        self[CONTENTS] = _contents
        self[DATA] = _data
        self[PRIORITY] = _priority
        return self

    def __setitem__(self, key, item):
        """ Protected [] mutation

        :param key: A string
        :param item: A value, which will be checked
        :return: None
        """
        if key not in MESSAGE_KEYS:
            raise KeyError(f"Message cannot have field {key}")
        if item is None:
            if key == PRIORITY:
                self._dictionary[key] = 8
            else:
                self._dictionary[key] = item
        else:
            if (key == TYPE) and (item not in MSG_TYPES):
                raise ValueError(f"Message type cannot be {key}, must be one of {MSG_TYPES}")
            if (key == RECIP) and (item not in ["AH", "LC", "SP", "MC", "_A", "_S"]):
                try:
                    item = rs_dir(item)
                except KeyError:
                    raise ValueError(f"Message recipient '{item}' must be on platform")
            if (key == CONTENTS) and not isinstance(item, str):
                raise ValueError(f"Message contents must be a string "
                                 f"(To pass objects, write to 'data' and specify data type in 'contents')")
            if key == PRIORITY:
                if not isinstance(item, int):
                    raise ValueError(f"Message priority must be [0,8], not {item}")
                item = 8 if item > 8 else 0 if item < 0 else item  #
            if key == IDKEY:
                if item == "bypass":
                    item = ID_KEY_BYPASS
                elif not isinstance(item, str):
                    item = self.gen_idemp_key(name=self[SENDER])
                elif bool(item):
                    # if the string has contents, no modification needed
                    # item = item
                    pass
                else:
                    item = self.gen_idemp_key(name=self[SENDER])
            self._dictionary[key] = item

    def __getitem__(self, key):
        """ Standard [] access

        Note [0] yields priority (implemented for priority queue)

        :param key: A string
        :return: The object stored in location [key]
        """
        if key == 0:
            return self._dictionary[PRIORITY]
        return self._dictionary[key]

    def get(self, key, default):
        try:
            return self._dictionary[key]
        except KeyError:
            return default

    # <Implemented for priority queue>

    def __lt__(self, other):
        return self[PRIORITY] < other[PRIORITY]

    def __le__(self, other):
        return self[PRIORITY] <= other[PRIORITY]

    def __gt__(self, other):
        return self[PRIORITY] > other[PRIORITY]

    def __ge__(self, other):
        return self[PRIORITY] >= other[PRIORITY]

    # </Implemented for priority queue>

    def encode(self):
        """ Converts a message into a bytes object for socket transmittal

        :return: A bytes object
        """
        msg_object = json.dumps(self._dictionary, default=lambda o: o.__dict__)
        binary_header = ""
        binary_header += self._dictionary.get(TYPE, "").rjust(3, "_")[0:3] + ">"
        binary_header += self._dictionary.get(RECIP, "").rjust(2, "_")[0:2] + ">"
        binary_header += self._dictionary.get(IDKEY, "").rjust(32, '0')[0:32] + ">"
        binary_header += f"{len(msg_object):023}" + ">"  # sys.getsizeof <--> len ?
        # print(f"DEBUG:: Msg {self[IDKEY]} given a size of {len(msg_object)} for '{msg_object}'")
        return (binary_header + msg_object).encode(FORMAT)

    def add_data(self, data=None, filename=None):
        """ Used to add the contents of the DATA field and handles loading files into a message

        Note: Providing a filename overrides the data argument

        :param data: String to be loaded into DATA
        :param filename: A file which will be converted into a string and loaded into DATA (takes priority)
        :return: None
        """
        if filename:
            with open(filename, 'rb') as file_handle:
                self[DATA] = base64.b64encode(file_handle.read()).decode(FORMAT)
        else:
            self[DATA] = data

    def read_data(self, filename=None):
        """ Used to fetch the contents of the DATA field and handles unloading files from a message

        :param filename: If given, will save the messaged file to 'filename'
        :return: If filename is given, returns True/False for success; otherwise, returns the string stored in DATA
        """
        if filename:
            try:
                with open(filename, 'wb+') as file_handle:
                    file_handle.write(base64.b64decode(self[DATA].encode(FORMAT)))
                return True
            except IOError as e:
                system_log.info(f"Failed to read contents of file {filename} due to: {repr(e)}")
                return False
        else:
            return self[DATA]

    def prints(self):
        """ A print-string for giving a short description of a Message

        :return: A string "type (sender>recipient) clipped-idemp-key: contents"
        """
        return f"{self[TYPE]} ({self[SENDER]}>{self[RECIP]}) {self[IDKEY][16:]}: '{self[CONTENTS]}' [{self[PRIORITY]}]"

    def print(self, sysio=True, trimmed=False):
        """ A print method for a Message

        Prints the contents for each key, but will slip DATA to the first 100 characters

        :param sysio: True - prints to console, otherwise prints to string
        :param trimmed: False - gives a tab indented list, True - removed tabs
        :return: A string
        """
        tab = "\n\t" if not trimmed else ", "
        msg_str = f"Message {self[IDKEY]}:"
        for k in self._dictionary.keys():
            data = self[k]
            if data is None:
                data = ""
            if (k == DATA) and len(data) > 256:
                msg_str += f"{tab}{k}: {data[:238]}...(clipped)"
            else:
                msg_str += f"{tab}{k}: {data}"
        msg_str = msg_str.replace(":, ", ": ", 1)
        if sysio:
            print(msg_str)
            return msg_str
        else:
            return msg_str

    def modified_copy(self, kw):
        """ Creates a copy of the message with some values changed

        Uses deepcopy to prevent cross-modification between the parent and copy Message

        :param kw: A dictionary of key-value pairs which shall overwrite the existing values
        :return:
        """
        msg_copy = deepcopy(self._dictionary)
        msg_copy = Message(msg_copy)
        for arg in kw:
            msg_copy[arg] = kw[arg]
        return msg_copy

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    @staticmethod
    def parse_header(binary_header):
        """ Decodes the bytes that lead a message through the socket into a list of pertinent data

        :param binary_header: A bytes object 64 bytes long
        :return: [Type, sender, idemp-key, length of message]
        """
        header_str = binary_header.decode(FORMAT)
        return header_str.split(">")[0:4]

    @staticmethod
    def load_message(binary_message):
        """ Decodes the bytes representing the message into a Message object

        :param binary_message: A bytes object of length specified by its header
        :return: A Message object
        """
        return Message(json.loads(binary_message.decode(FORMAT)))

    @staticmethod
    def gen_idemp_key(name=""):
        """ Generates idempotency keys

        :param name: Name of the system issuing the Message
        :return: A 32-character string of the form (2:sender name)(14:timestamp)(16:unique characters)
        """
        name = "" if name is None else name
        source = ID_KEY_SOURCE
        ret_string = name.rjust(2, "_")[0:2]
        ret_string += datetime.now().strftime("%Y%m%d%H%M%S")
        for _ in range(16):
            ret_string += source[random.randint(0, len(source)-1)]
        return ret_string

    @staticmethod
    def connect(name):
        """ Generates the message recognized by the server as a connect signal

        :param name: Name of the system connecting to the server
        :return: A Message object
        """
        cnt_msg = {TYPE: 'COM', SENDER: name, RECIP: "MC", IDKEY: "bypass", CONTENTS: "connect_signal", DATA: None,
                   PRIORITY: P_SERVER}
        return Message(cnt_msg)

    @staticmethod
    def disconnect(name):
        """ Generates the message recognized by the server as a disconnect signal

        Disambiguation with leave_server():  leave_server() is a System-level command which results in the System
        sending a signal to the server that it is leaving.  That signal is disconnect().  disconnect() if for System-
        to-Server and leave_server() is for System-to-System.

        :param name: Name of the system disconnecting from the server
        :return: A Message object
        """
        dcn_msg = {TYPE: 'COM', SENDER: name, RECIP: "MC", IDKEY: "bypass", CONTENTS: "disconnect_signal", DATA: None,
                   PRIORITY: P_SERVER}
        return Message(dcn_msg)

    @staticmethod
    def initialization_message(name):
        """ Creates a message which signals the initializaton of a system

        :param name: The name of the Agent
        :return: A Message
        """
        init_msg = {TYPE: "CMD", SENDER: name, RECIP: name, IDKEY: "",
                    CONTENTS: "ib.Read_Operation",
                    DATA: oprtn.Operation({oprtn.FUNC: "Initialize", oprtn.AGENT: name}).package(),
                    PRIORITY: 0}
        return Message(init_msg)

    @staticmethod
    def close_server():
        """ Generates the message recognized by the server as a signal to terminate

        :return: A Message object
        """
        end_msg = {TYPE: 'COM', SENDER: "MC", RECIP: "__", IDKEY: "bypass", CONTENTS: "close_server", DATA: None,
                   PRIORITY: P_SERVER}
        return Message(end_msg)

    @staticmethod
    def request_of_state(*, requester="", requestee=""):
        """ Generates a message the requestee will read and send its status back to requester

        :param requester: The system that wants the requestee's status
        :param requestee: The system whose status is being polled
        :return: A Message object
        """
        req_msg = {TYPE: "CMD", SENDER: requester, RECIP: requestee, IDKEY: "",
                   CONTENTS: "ib.Read_Operation",
                   DATA: oprtn.Operation({
                       oprtn.FUNC: "__.Send_State_to",
                       oprtn.KWARGS: {KWA_REQUESTER: requester},
                   }).package(),
                   PRIORITY: P_SERVER}
        return Message(req_msg)

    @staticmethod
    def status_report(status, *, sender="", recipient=""):
        """ Generates a message the recipient will read as a status update

        :param status: A str or json-made dict which represents the status update
        :param sender: The system whose status is in the status field
        :param recipient: The system who requested the status update
        :return: A Message object
        """
        if not isinstance(status, str):
            status = json.dumps(status, default=lambda o: o.__dict__())  # , default=lambda o: o.__dict__
        sts_msg = {TYPE: "CMD", SENDER: sender, RECIP: recipient, IDKEY: "", CONTENTS: "ib.Read_Operation",
                   DATA: oprtn.Operation({
                       oprtn.FUNC: "__.Read_Status_Report",
                       oprtn.KWARGS: {KWA_STS_REPORT: status}
                   }).package(),
                   PRIORITY: P_CMD_NR}
        return Message(sts_msg)

    @staticmethod
    def note_network_change(member_list: list):
        """ Generates a message the MC will recognize as an indication that a member has disconnected

        (Note: Used by the server upon a socket disconnect)

        :param member_list: Name of the system that has disconnected
        :return: A Message object
        """
        req_msg = {TYPE: "CMD", SENDER: '_S', RECIP: '_A', IDKEY: "",
                   CONTENTS: "ib.Read_Operation",
                   DATA: oprtn.Operation({
                       oprtn.FUNC: "__.Note_Network_Change",
                       oprtn.KWARGS: {KWA_NET_MEMBERS: member_list},
                   }).package(),
                   PRIORITY: P_SERVER}
        return Message(req_msg)

    @staticmethod
    def confirm_ro(name, recip, completed_message_id_key, *,
                   was_successful: bool, location=None, level: str = None, info: str = None):
        """ Message which is read as a confirmation of receipt and operation of a task

        :param name: System who has completed the task
        :param recip: System who issued the task (MC)
        :param completed_message_id_key: The unique ID for the issued task
        :param was_successful: True/False - Completed without issue
        :param location: If issue, the location of the problem
        :param level: If issue, the severity of the problem
        :param info: If issue, details
        :return: A Message object
        """
        sts_msg = {TYPE: "CMD", SENDER: name, RECIP: recip, IDKEY: "",
                   CONTENTS: "ib.Confirmation",
                   DATA: json.dumps({KWA_TASK_KEY: completed_message_id_key,
                                     KWA_TASK_RESULT: {'completion': was_successful,
                                                       'location': location,
                                                       'level': level,
                                                       'data': info}
                                     }),
                   PRIORITY: P_CMD_NR}
        return Message(sts_msg)

    @staticmethod
    def add_fault(sender, system_with_fault, fault):
        """ Message to add a fault

        :param sender: The reporter of the fault
        :param system_with_fault: The system to whom the fault belongs
        :param fault: the Fault
        :return: A Message
        """
        req_msg = {TYPE: "CMD", SENDER: sender, RECIP: system_with_fault, IDKEY: "",
                   CONTENTS: "ib.Read_Operation",
                   DATA: oprtn.Operation({
                       oprtn.FUNC: "__.Add_Fault",
                       oprtn.KWARGS: {KWA_FAULT_VAL: fault},
                   }).package(),
                   PRIORITY: P_CMD_NR}
        return Message(req_msg)

    @staticmethod
    def remove_fault(sender, system_removing_fault, fault):
        """ Message to remove a fault

        :param sender: The reporter of fault resolution.
        :param system_removing_fault: The system with the fault.
        :param fault: the Fault.
        :return: A Message
        """
        req_msg = {TYPE: "CMD", SENDER: sender, RECIP: system_removing_fault, IDKEY: "",
                   CONTENTS: "ib.Read_Operation",
                   DATA: oprtn.Operation({
                       oprtn.FUNC: "__.Resolve_Fault",
                       oprtn.KWARGS: {KWA_FAULT_VAL: fault},
                   }).package(),
                   PRIORITY: P_CMD_NR}
        return Message(req_msg)

    @staticmethod
    def resuscitation_message(sender, recipient, checkpoint_key, q_name, q_step, update_interval=None, timeout=None):
        """ Message to resume monitoring responses for an operation's checkpoint

        :param sender: The sender
        :param recipient: The system with the checkpoint
        :param checkpoint_key: ibid
        :param q_name: the name of the associated queue
        :param q_step: the step number of the associated operation
        :param update_interval: the interval with which to check for updates
        :param timeout: the new timeout field
        :return: A Message
        """
        req_msg = {TYPE: "CMD", SENDER: sender, RECIP: recipient, IDKEY: "",
                   CONTENTS: "ib.Read_Operation",
                   DATA: oprtn.Operation({
                       oprtn.FUNC: "ib.Resuscitate",
                       oprtn.KWARGS: {'checkpoint_key': checkpoint_key,
                                      'queue_name': q_name,
                                      'queue_step': q_step,
                                      'update_interval': update_interval,
                                      'timeout': timeout},
                   }).package(),
                   PRIORITY: P_CMD_NR}
        return Message(req_msg)

    @staticmethod
    def kill():
        """ A kill signal (still being worked out)

        :return: A Message object
        """
        kill_msg = {TYPE: 'CMD', SENDER: "__", RECIP: "_A", IDKEY: "bypass",
                    CONTENTS: None, DATA: None, PRIORITY: P_EMERGENCY}
        # script = Script("kill_signal")
        # script.add_operations({FUNC: '__.Exit'})
        # kill_msg[CONTENTS] = script
        return Message(kill_msg)

    @staticmethod
    def ping(recip="_S"):
        """ Creates a Ping for the server

        :param recip: The system being pinged
        :return: A Message
        """
        ping_msg = {TYPE: 'COM', SENDER: "__", RECIP: recip, IDKEY: "bypass",
                    CONTENTS: "ping", DATA: None, PRIORITY: P_SERVER}
        return Message(ping_msg)

    @staticmethod
    def leave_server(recip):
        """ Message signaling a member to disconnect (a command)

        Disambiguation: disconnect() is for messages to the server--it is a signal--, leave_server() is for the Agent--
        it is an operation.  leave_server() tells an Agent to send a disconnect() to the Server.

        :param recip: The Agent which should disconnect
        :return: A Message
        """
        req_msg = {TYPE: "CMD", SENDER: "_S", RECIP: recip, IDKEY: "bypass",
                   CONTENTS: "ib.Read_Operation",
                   DATA: oprtn.Operation({
                       oprtn.FUNC: "__.Disconnect",
                   }).package(),
                   PRIORITY: P_SERVER}
        return Message(req_msg)
