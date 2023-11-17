""" Module to allow a user to make their own Messages
@author: Ben C
"""

import tkinter
from tkinter import filedialog
from requests import exceptions as request_exceptions

import functionals
import json
import message as msgm
from constants import *
from database_constants import *
from gui_constants import *
from mcn_queues import MCNPriorityQueue
import custom_exceptions as cexc
import operations as oprtn
import database_interface as dbi
from mcn_logging_manager import system_log
from mcn_status import rs_dir, MCN_CFG, S_ALL_SYS

CONTACTS = functionals.list_subtract(MCN_CFG[S_ALL_SYS], ["_S", "_U", "__", "ui"])
ALLOWED_CONTENT_TYPES = ["Connect", "Disconnect",
                         "ib.Read_String", "ib.Read_File", "ib.Read_Operation",
                         "__.Run_From_DB", "ib.Run_Local"]


class MessageMaker:
    """ A GUI for Message Making """
    def __init__(self, master, data_pipes, cmd_list_keys, name="Not defined"):
        """ Creates a popup for making and submitting messages that adapts to the message being written and the System
        it is written from.

        :param master: a tkinter root
        :param data_pipes: For communication with the Parent System
        :param cmd_list_keys: To know the functionality of the System
        :param name: the name of the System
        """
        self.inbox = data_pipes[MSG_Q_IN]
        self.outbox = data_pipes[MSG_Q_OUT]
        self.name = name
        self.cmd_list_str = list(cmd_list_keys)
        self.core = master

        self.frame = tkinter.Frame(self.core)
        self.frame.winfo_toplevel().title("Message Maker")
        self.frame.pack(expand=True, side=tkinter.LEFT, fill=tkinter.BOTH, padx=10)
        self.labelframe_message = tkinter.LabelFrame(self.frame, text="Message Maker")
        self.labelframe_message.pack(**TB33)
        self.labelframe_data = tkinter.LabelFrame(self.frame, text="Data")
        self.labelframe_data.pack(expand=True, **TB33)
        self.submit_button = tkinter.Button(self.frame, text="Submit", command=self.submit, state=tkinter.DISABLED)
        self.submit_button.pack(**TB33)

        self.sender = tkinter.StringVar(self.labelframe_message)
        self.recipient = tkinter.StringVar(self.labelframe_message, value="MC")
        self.id_key = tkinter.StringVar(self.labelframe_message)
        self.contents = tkinter.StringVar(self.labelframe_message)
        self.priority = tkinter.IntVar(self.labelframe_message)

        self.selected_queue = tkinter.StringVar(self.labelframe_message)

        self.operation_func = tkinter.StringVar(self.labelframe_data)
        self.operation_agent = tkinter.StringVar(self.labelframe_data)

        self.msg_form = dict()
        self.data_form = dict()
        self.selected_file = ""

        self.build_gui()

    def build_gui(self):
        """ Constructs the message-making Form

        :return:
        """
        # row 1 - MSG sender
        self.msg_form[0] = dict()
        self.msg_form[0][0] = tkinter.Label(self.labelframe_message, text="Sender")

        self.msg_form[0][1] = tkinter.Radiobutton(self.labelframe_message, text='Self', indicatoron=0,  # noqa
                                                  variable=self.sender, value=self.name,
                                                  command=lambda: self.msg_sender_selector('self'))
        self.msg_form[0][1].select()
        self.msg_form[0][2] = tkinter.Entry(self.labelframe_message)
        self.msg_form[0][2].bind("<1>", lambda x: self.msg_sender_selector(x))

        # row 2 - MSG recipient
        self.msg_form[1] = dict()
        self.msg_form[1][0] = tkinter.Label(self.labelframe_message, text="Recipient")
        self.msg_form[1][1] = tkinter.OptionMenu(self.labelframe_message, self.recipient, *CONTACTS)

        # row 3 - IDKey
        self.msg_form[2] = dict()
        self.msg_form[2][0] = tkinter.Label(self.labelframe_message, text="ID Key")
        self.msg_form[2][1] = tkinter.Radiobutton(self.labelframe_message, text='Auto',
                                                  variable=self.id_key, value="auto")
        self.msg_form[2][1].select()
        self.msg_form[2][2] = tkinter.Radiobutton(self.labelframe_message, text='Bypass',
                                                  variable=self.id_key, value="bypass")

        # row 4 - Contents
        self.msg_form[3] = dict()

        self.msg_form[3][0] = tkinter.Label(self.labelframe_message, text="Contents")
        self.msg_form[3][1] = tkinter.OptionMenu(self.labelframe_message, self.contents,
                                                 *ALLOWED_CONTENT_TYPES, command=self.update_data_frame)

        # row 5 - Priority
        self.msg_form[4] = dict()
        self.msg_form[4][0] = tkinter.Label(self.labelframe_message, text="Priority")
        self.msg_form[4][1] = tkinter.Radiobutton(self.labelframe_message, text='Normal',
                                                  variable=self.priority, value=msgm.P_CMD_NR)
        self.msg_form[4][1].select()
        self.msg_form[4][2] = tkinter.Radiobutton(self.labelframe_message, text='High',
                                                  variable=self.priority, value=msgm.P_EMERGENCY)

        for r in self.msg_form.keys():
            for c in self.msg_form[r].keys():
                if r == 3 and c == 1:
                    self.msg_form[r][c].grid(column=c, row=r, columnspan=2, sticky=tkinter.W)
                else:
                    self.msg_form[r][c].grid(column=c, row=r, sticky=tkinter.W)

        self.data_form[0] = tkinter.Label(self.labelframe_data, text="    ")
        self.data_form[0].pack()

    def submit(self, *_):
        """ Converts the form into a Message object and sends it

        Closes on completion

        :return:
        """
        recipient = self.recipient.get()
        id_key = self.id_key.get()
        id_key = "" if id_key == 'auto' else id_key
        priority = self.priority.get()
        operation_func = self.operation_func.get()
        operation_agent = self.operation_agent.get()
        contents_type = self.contents.get()

        if contents_type == 'Connect':
            my_message = msgm.Message.connect(self.name)
            self.inbox.enqueue(my_message)
        elif contents_type == 'Disconnect':
            my_message = msgm.Message.disconnect(self.name)
            self.outbox.enqueue(my_message)
        elif contents_type == "ib.Read_String":
            my_message = msgm.Message.build_from_args(_type="CMD",
                                                      _sender=self.get_sender(),
                                                      _recip=recipient,
                                                      _idkey=id_key,
                                                      _contents=contents_type,
                                                      _data=self.data_form[1].get(),
                                                      _priority=priority)
            self.send(recipient, my_message)
        elif contents_type == "ib.Read_File":
            my_message = msgm.Message.build_from_args(_type="CMD",
                                                      _sender=self.get_sender(),
                                                      _recip=recipient,
                                                      _idkey=id_key,
                                                      _contents=contents_type,
                                                      _data=None,
                                                      _priority=priority)
            my_message.add_data(filename=self.selected_file)
            self.send(recipient, my_message)
        elif contents_type == "__.Run_From_DB":
            q_name = self.selected_queue.get()
            step_num = self.data_form[3].get()
            # Grab the queue
            try:
                decoded_content, _, _ = dbi.query_document(self.name, 'queue', 'queue_name', q_name)
            except request_exceptions.ConnectionError:
                system_log.exception("GUI - DB request failed - Queue lookup")
                return
            # Add supplementary materials
            try:
                print(decoded_content[Q_OPERATIONS_LIST][step_num])
                decoded_content[Q_OPERATIONS_LIST][step_num][Q_NAME] = q_name
                decoded_content[Q_OPERATIONS_LIST][step_num][DBQ_STEP_NUM] = step_num
                task = decoded_content[Q_OPERATIONS_LIST][step_num]
            except KeyError:
                system_log.exception(f"Encountered exception; step number likely invalid")
                return
            # Make message
            my_message = msgm.Message.build_from_args(
                "CMD", self.name, task[QOP_AGENT], id_key, "__.Run_From_DB",
                oprtn.Operation({oprtn.FUNC: task[QOP_OPERATION],
                                 oprtn.KWARGS: {Q_NAME: task[Q_NAME],
                                                oprtn.PID: task[QOP_CONTAINER],
                                                DBQ_STEP_NUM: task[DBQ_STEP_NUM]},
                                 oprtn.AGENT: task[QOP_AGENT]
                                 }).package(),
                msgm.P_CMD_NR)
            # Send message
            self.send(recipient, my_message)
        elif contents_type == "__.Run_Local":
            my_message = msgm.Message.build_from_args(_type="CMD",
                                                      _sender=self.get_sender(),
                                                      _recip=recipient,
                                                      _idkey=id_key,
                                                      _contents="ib.Read_Operation",
                                                      _data=None,
                                                      _priority=priority)
            operation_args = self.data_form[3].get()
            operation_args = tuple(json.loads(operation_args)) if operation_args else None
            operation_kwargs = self.data_form[5].get()
            operation_kwargs = json.loads(operation_kwargs) if operation_kwargs else None
            my_operation = oprtn.Operation.build_from_args(func=operation_func,
                                                           args=operation_args,
                                                           kwargs=operation_kwargs,
                                                           agent=operation_agent)
            my_message.add_data(data=my_operation.package())
            self.send(recipient, my_message)
        else:
            my_message = msgm.Message.build_from_args(_type="CMD",
                                                      _sender=self.get_sender(),
                                                      _recip=recipient,
                                                      _idkey=id_key,
                                                      _contents=contents_type,
                                                      _data=None,
                                                      _priority=priority)
            operation_args = self.data_form[3].get()
            operation_args = tuple(json.loads(operation_args)) if operation_args else None
            operation_kwargs = self.data_form[5].get()
            operation_kwargs = json.loads(operation_kwargs) if operation_kwargs else None
            my_operation = oprtn.Operation.build_from_args(func=operation_func,
                                                           args=operation_args,
                                                           kwargs=operation_kwargs,
                                                           agent=operation_agent)
            my_message.add_data(data=my_operation.package())
            self.send(recipient, my_message)
        self.core.destroy()

    def msg_sender_selector(self, *args):
        """ Handles the use of the 'self' as sender vs specify-other-as-sender feature

        :param args: Name of the sending system
        :return:
        """
        if args[0] == 'self':
            self.msg_form[0][2].delete(*INPUT_CONTENTS)
        else:
            self.msg_form[0][1].deselect()

    def select_file(self, *_):
        """ Opens the FTP dialog

        :return:
        """
        self.selected_file = filedialog.askopenfilename(title='Open a file', initialdir='/')

    def update_data_frame(self, *_):
        """ Updates the form for message making depending on the Message Contents specification

        :return:
        """
        self.submit_button["state"] = "active"
        for item in self.labelframe_data.winfo_children():
            item.destroy()
        contents_type = self.contents.get()

        if contents_type == '__.Run_From_DB':
            self.data_form[0] = tkinter.Label(self.labelframe_data, text="Queue Name:")
            self.data_form[0].grid(row=0, column=0, sticky=tkinter.W)
            try:
                self.data_form[1] = tkinter.OptionMenu(self.labelframe_data, self.selected_queue,
                                                       *dbi.get_available_queue_names())
            except cexc.DatabaseGeneralException:
                system_log.exception("Failed to pull a list of available queue names")
                self.data_form[1] = tkinter.OptionMenu(self.labelframe_data, self.selected_queue,
                                                       *["DB Not Connected!", ])
            self.data_form[1].grid(row=0, column=1, sticky=tkinter.W)
            self.data_form[2] = tkinter.Label(self.labelframe_data, text="Step #:")
            self.data_form[2].grid(row=1, column=0, sticky=tkinter.W)
            self.data_form[3] = tkinter.Entry(self.labelframe_data)
            self.data_form[3].grid(row=1, column=1, sticky=tkinter.W)
        elif contents_type == 'Connect':
            self.data_form[0] = tkinter.Label(self.labelframe_data, text="Do not use if already connected")
            self.data_form[0].pack()
        elif contents_type == 'Disconnect':
            self.data_form[0] = tkinter.Label(self.labelframe_data, text="Do not use if already disconnected")
            self.data_form[0].pack()
        elif contents_type == "ib.Read_String":
            self.data_form[0] = tkinter.Label(self.labelframe_data, text="Enter String")
            self.data_form[0].pack(**LB33)
            self.data_form[1] = tkinter.Entry(self.labelframe_data)
            self.data_form[1].pack(**LB33)
        elif contents_type == "ib.Read_File":
            self.data_form[0] = tkinter.Button(self.labelframe_data, text="Select File", command=self.select_file)
            self.data_form[0].pack(**TB33)
        else:
            self.labelframe_data.rowconfigure((1,2), weight=1)  # noqa: it (a) works and (b) is used in their docs...
            self.labelframe_data.columnconfigure(1, weight=1)
            self.data_form[0] = tkinter.Label(self.labelframe_data, text="Function")
            self.data_form[0].grid(row=0, column=0, sticky=tkinter.W)
            self.data_form[1] = tkinter.OptionMenu(self.labelframe_data, self.operation_func, *self.cmd_list_str)
            self.data_form[1].grid(row=0, column=1, sticky=tkinter.W)
            self.data_form[2] = tkinter.Label(self.labelframe_data, text="Args")
            self.data_form[2].grid(row=1, column=0, sticky=tkinter.W)
            self.data_form[3] = tkinter.Entry(self.labelframe_data, width=30)
            self.data_form[3].grid(row=1, column=1, sticky='nsew')
            self.data_form[4] = tkinter.Label(self.labelframe_data, text="Kwargs")
            self.data_form[4].grid(row=2, column=0, sticky=tkinter.W)
            self.data_form[5] = tkinter.Entry(self.labelframe_data, width=3)
            self.data_form[5].grid(row=2, column=1, sticky='nsew')
            self.data_form[6] = tkinter.Label(self.labelframe_data, text="Agent")
            self.data_form[6].grid(row=3, column=0, sticky=tkinter.W)
            self.operation_agent.set(self.recipient.get())
            self.data_form[7] = tkinter.OptionMenu(self.labelframe_data, self.operation_agent, *CONTACTS)
            self.data_form[7].grid(row=3, column=1, sticky=tkinter.W)

    def get_sender(self):
        """ Retrieves the sender

        :return:
        """
        sender = self.msg_form[0][2].get()
        if not sender:
            sender = self.sender.get()
        if not sender:
            print("Problem")
            sender = "__"
        return sender

    def send(self, recip, msg):
        """ Sends a message by putting it into the inbox/outbox

        :param recip: The recipient
        :param msg: the Message
        :return: None
        """
        if recip == "_A":
            self.outbox.enqueue(msg)
        elif rs_dir(recip) == self.name:
            self.inbox.enqueue(msg)
        else:
            self.outbox.enqueue(msg)


if __name__ == '__main__':
    import queue
    core = tkinter.Tk()

    command_list = dict()
    command_list.update({
        "__.Connect": None,
        "__.Disconnect": None,
        "__.Note_Disconnect": None,
        "__.Exit": None,
        "__.Resolve_Fault": None,
        "__.Add_Fault": None,
        "__.Send_State_to": None,
        "__.Read_Status_Report": None,
        "__.Run_From_DB": None,
        "ib.Confirmation": None,
        "__.Echo": None,
        "__.Except": None,
        "__.Wait": None,
        "ib.Read_String": None,
        "ib.Read_File": None,
        "ib.Read_Operation": None,
        "ib.Read_Script": None,
    })

    boxes = {MSG_Q_OUT: MCNPriorityQueue(), MSG_Q_IN: MCNPriorityQueue()}

    my_msg_maker = MessageMaker(tkinter.Toplevel(core),
                                boxes,
                                command_list.keys(),
                                "MC")

    core.mainloop()

    try:
        print(boxes[MSG_Q_OUT].get(block=False).print())
    except queue.Empty:
        pass
    try:
        print(boxes[MSG_Q_IN].get(block=False).print())
    except queue.Empty:
        pass
