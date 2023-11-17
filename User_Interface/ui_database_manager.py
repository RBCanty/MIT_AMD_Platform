""" Database Manager UI
@author: Ben C
"""

import tkinter as tk
from datetime import datetime
from io import StringIO
from pprint import pformat

import custom_exceptions as cexc
from custom_classes import Narg
import functionals
import database_interface as dbi
import yaml
import json
from database_constants import *
from gui_constants import *
from mcn_logging_manager import system_log
from ui_mea_culpa import Arbitrator

DEFAULT_OPERATION = {"-1": {QOP_AGENT: 'N.a.', QOP_OPERATION: 'N.a.', QOP_COMPLETED: 'N.a.'}}
QUICK_WIDTH = 64
TABULA_WIDTH = 140
ADV_MENU_WIDTH = 24


class DatabaseManager:
    """ A Simple UI for viewing Queue information """
    def __init__(self, master, name="Not defined"):
        """ Creates a UI for viewing Queue operation information

        :param master: a tkinter root
        :param name: Name of the System (for logon)
        """
        self.core = master
        self.name = name
        self.tk = dict()
        self.reg = dict()

        self.build_opr()
        self.build_details()

    def build_opr(self):
        """ Constructs and populates the Frames/Widgets for selecting queues and operations

        :return:
        """
        self.tk['opr'] = dict()
        self.tk['opr']['frame'] = tk.Frame(self.core)
        self.tk['opr']['frame'].pack(expand=True, **LBA3)
        self.tk['opr']['update'] = tk.Button(self.tk['opr']['frame'], text="Pull DB", command=self.update)
        self.tk['opr']['update'].pack(**TB33)  # .grid(row=0, column=0, **GRID_B33)

        self.tk['opr']['frame'].winfo_toplevel().title("Database Manager")

        self.reg['opr.q_menu'] = tk.StringVar(self.tk['opr']['frame'], value='Select Queue')
        try:
            self.reg['opr.q_list'] = ['Select Queue'] + dbi.get_available_queue_names()
        except cexc.DatabaseGeneralException:
            system_log.warning("Could not populate list of available queues", exc_info=True)
            self.reg['opr.q_list'] = ['Select Queue']
        self.tk['opr']['q_menu'] = tk.OptionMenu(self.tk['opr']['frame'], self.reg['opr.q_menu'],
                                                 *self.reg['opr.q_list'], command=self.q_select)
        self.tk['opr']['q_menu'].pack(**TB33)  # .grid(row=1, column=0, **GRID_B33)
        self.tk['opr']['disp'] = dict()
        self.tk['opr']['disp']['frame'] = tk.Frame(self.tk['opr']['frame'])
        self.tk['opr']['disp']['frame'].pack(expand=True, **TBA3)  # .grid(row=2, column=0, **GRID_B33)
        self.tk['opr']['disp']['scroll'] = tk.Scrollbar(self.tk['opr']['disp']['frame'])
        self.tk['opr']['disp']['scroll'].pack(side=tk.RIGHT, fill=tk.Y)
        self.tk['opr']['disp']['box'] = tk.Listbox(self.tk['opr']['disp']['frame'],
                                                   yscrollcommand=self.tk['opr']['disp']['scroll'].set,
                                                   width=32)
        self.tk['opr']['disp']['box'].bind('<<ListboxSelect>>', self.show_details)
        self.tk['opr']['disp']['box'].pack(expand=True, **TB33)
        self.tk['opr']['disp']['scroll'].config(command=self.tk['opr']['disp']['box'].yview)

    def build_details(self):
        """ Creates and populates the Frames/Widgets for viewing and (basic) editing of operation properties

        :return:
        """
        self.tk['det'] = dict()
        self.tk['det']['frame'] = tk.Frame(self.core)
        self.tk['det']['frame'].pack(**LBA3)
        self.tk['det']['info'] = tk.Text(self.tk['det']['frame'], width=QUICK_WIDTH, height=12)
        self.tk['det']['info'].insert(tk.END, "Select a queue and step for details")
        self.tk['det']['info'].pack(expand=True, **TB33)
        self.tk['det']['cmd'] = dict()
        self.tk['det']['cmd']['frame'] = tk.Frame(self.tk['det']['frame'])
        self.tk['det']['cmd']['frame'].pack(**TB33)
        self.tk['det']['cmd']['y'] = tk.Button(self.tk['det']['cmd']['frame'],
                                               text="Mark Complete",
                                               command=lambda: self.mark('yes'))
        self.tk['det']['cmd']['y'].pack(**LB33, padx=3, pady=6)
        self.tk['det']['cmd']['n'] = tk.Button(self.tk['det']['cmd']['frame'],
                                               text="Mark Incomplete",
                                               command=lambda: self.mark('no'))
        self.tk['det']['cmd']['n'].pack(**LB33, padx=3, pady=6)
        self.tk['det']['cmd']['a'] = tk.Button(self.tk['det']['cmd']['frame'],
                                               text="Advanced Editor",
                                               command=self.advanced)
        self.tk['det']['cmd']['a'].pack(**LB33, padx=3, pady=6)
        self.tk['det']['cmd']['b'] = tk.Button(self.tk['det']['cmd']['frame'],
                                               text="Fault Editor",
                                               command=self.fault_editor)
        self.tk['det']['cmd']['b'].pack(**LB33, padx=3, pady=6)

    def update(self, *_):
        """ Pulls the database for updates

        :return:
        """
        try:
            temp, _, _ = dbi.query_collection('MC', 'queue')
        except cexc.DatabaseRequestError:
            return
        if temp:
            self.reg['queue_collection'] = temp
        else:
            self.reg['queue_collection'] = dict()

        self.reg['opr.q_list'] = ['Select Queue'] + [v[Q_NAME]
                                                     for _, v in self.reg['queue_collection'].items()
                                                     if v.get(Q_NAME, False)]
        self.tk['opr']['q_menu']['menu'].delete(0, tk.END)
        for q in self.reg['opr.q_list']:
            self.tk['opr']['q_menu']['menu'].add_command(
                label=q,
                command=lambda v=q: [self.reg['opr.q_menu'].set(v), self.q_select()]
            )

    def get_sel_q(self):
        """ Retrieves the selected queue from the menu

        :return:
        """
        sel_q = self.reg['opr.q_menu'].get()

        if sel_q == 'Select Queue':
            return None

        return sel_q

    def q_select(self, *_):
        """ Updates operation list display when a queue is selected

        :return:
        """
        sel_q = self.get_sel_q()

        if sel_q is None:
            return

        self.tk['det']['info'].delete('1.0', tk.END)

        self.update()
        doc = functionals.search_dict_by_value(self.reg['queue_collection'], Q_NAME, sel_q)
        if doc is None:
            self.update()
            return

        # o: operations_list
        o = doc.get(Q_OPERATIONS_LIST, DEFAULT_OPERATION)
        # sk: sorted_keys
        sk = [str(i) for i in sorted([int(k) for k in o.keys()])]
        # pl: pretty_list
        try:
            pl = [f"{k} {o[k][QOP_AGENT]}:{o[k][QOP_OPERATION]} ({o[k][QOP_COMPLETED]})" for k in sk]
        except KeyError:
            system_log.exception(f"Error building list of operations for queue {sel_q}")
            return

        self.tk['opr']['disp']['box'].delete(0, tk.END)
        for item in pl:
            self.tk['opr']['disp']['box'].insert(tk.END, item)

    def show_details(self, *_):
        """ Updates operation details display when an operation is selected

        :return:
        """
        index = self.tk['opr']['disp']['box'].curselection()[:1]

        if not index:
            return

        sel = self.tk['opr']['disp']['box'].get(index)

        n, *_ = str(sel).split(" ")

        sel_q = self.get_sel_q()

        if sel_q is None:
            return

        doc = functionals.search_dict_by_value(self.reg['queue_collection'], Q_NAME, sel_q)

        self.tk['det']['info'].delete('1.0', tk.END)
        self.tk['det']['info'].insert(tk.END, f"{n}\n" + pformat(doc[Q_OPERATIONS_LIST][n], width=QUICK_WIDTH // 2))
        return pformat(doc[Q_OPERATIONS_LIST][n])

    def mark(self, val, *_):
        """ Calls updates to database to mark an operation as complete/incomplete

        :param val: yes/no - is the operation complete
        :return: None
        """
        index = self.tk['opr']['disp']['box'].curselection()[:1]

        if not index:
            return

        sel = self.tk['opr']['disp']['box'].get(index)

        n, *_ = str(sel).split(" ")

        sel_q = self.get_sel_q()

        if sel_q is None:
            return

        if val == 'yes':
            modus = {'stop': True}
        elif val == 'no':
            modus = {'reset': True}
        else:
            modus = None

        doc = functionals.search_dict_by_value(self.reg['queue_collection'], Q_NAME, sel_q)
        source = doc[Q_OPERATIONS_LIST][n][QOP_COMPLETED]

        upd8_request = {'request_type': 'update_database_field',
                        'document_name': sel_q,
                        'update_field': f"{Q_OPERATIONS_LIST}.{n}.{QOP_COMPLETED}",
                        'old_value': source,
                        'new_value': val
                        }

        # Note extra arg: this deactivates the raises and specifies a detailed return value
        try:
            dbi.database_request(self.name, upd8_request, True)
            if modus:
                dbi.mark_time(sel_q, n, **modus, logger=system_log)
        except cexc.DatabaseGeneralException:
            system_log.exception(f"Failed to make update request to Database")
        else:
            self.update()
            self.q_select()
            self.show_details()

    def advanced(self, *_):
        """ Opens the Advanced Database Manager GUI for viewing and editing DB documents

        :return:
        """
        AdvancedDatabaseManager(tk.Toplevel(self.core), self.name)

    def fault_editor(self, *_):
        """ Opens the Fault Editing GUI for viewing/adding/editing/removing Faults from queue documents

        :return:
        """
        Arbitrator(tk.Toplevel(self.core))


class AdvancedDatabaseManager:
    """ A more rigorous Database getter/setter """
    def __init__(self, root, name="Not defined"):
        """ Creates a GUI to interface with all aspects of the Database

        :param root: a tkinter root
        :param name: the name of the System (for logon)
        """
        self.core = root
        self.name = name
        self.tk = dict()
        self.reg = dict()
        self.var_args = dict()
        self.janus = True
        self.cursor = 0
        self.finds = []

        self.build_gui()

    def build_gui(self):
        """ Constructs the Widgets for making DB requests and viewing results

        :return:
        """
        self.tk['ctrl'] = dict()
        self.tk['ctrl']['frame'] = tk.Frame(self.core)
        self.tk['ctrl']['frame'].pack(**LB33)
        self.tk['ctrl']['frame'].winfo_toplevel().title("Advanced Database Manager")
        self.tk['ctrl']['req_type_label'] = tk.Label(self.tk['ctrl']['frame'], text="Request Type:", anchor=tk.W)
        self.tk['ctrl']['req_type_label'].pack(**TB33)
        self.reg['ctrl.req_type'] = tk.StringVar(self.tk['ctrl']['frame'], REQ_TYPES[0])
        self.tk['ctrl']['req_type_menu'] = tk.OptionMenu(self.tk['ctrl']['frame'],
                                                         self.reg['ctrl.req_type'],
                                                         *REQ_TYPES,
                                                         command=self.build_var_gui)
        self.tk['ctrl']['req_type_menu'].config(width=ADV_MENU_WIDTH)
        self.tk['ctrl']['req_type_menu'].pack(**TB33)
        self.tk['ctrl']['var_frame'] = tk.Frame(self.tk['ctrl']['frame'])
        self.tk['ctrl']['var_frame'].pack(**TB33)
        self.tk['ctrl']['submit'] = tk.Button(self.tk['ctrl']['frame'], text="Submit", command=self.submit)
        self.tk['ctrl']['submit'].pack(**TB33)
        self.tk['ctrl']['info'] = tk.Label(self.tk['ctrl']['frame'], text="--")
        self.tk['ctrl']['info'].pack(**TB33)

        tk.Button(self.tk['ctrl']['frame'], text="Reformat", command=self.reformat).pack(**TB33)
        tk.Label(self.tk['ctrl']['frame'],
                 text="When updating documents,\ndon't forget to trim the id").pack(**TB33)
        tk.Button(self.tk['ctrl']['frame'],
                  text="Reopen Document", command=self.reopen).pack(**TB33)
        tk.Button(self.tk['ctrl']['frame'],
                  text="Unload Saved Documents", command=self.clear_reg).pack(**TB33)

        tk.Label(self.tk['ctrl']['frame'],
                 text="Find").pack(**TB33)
        self.tk['ctrl']['find_field'] = tk.Entry(self.tk['ctrl']['frame'])
        self.tk['ctrl']['find_field'].bind("<1>", self._clear_tag)
        self.tk['ctrl']['find_field'].pack(**TB33)
        self.tk['ctrl']['find_button_1'] = tk.Button(self.tk['ctrl']['frame'],
                                                     text="Find Next", command=lambda: self.find(1))
        self.tk['ctrl']['find_button_1'].pack(**TB33)
        self.tk['ctrl']['find_button_2'] = tk.Button(self.tk['ctrl']['frame'],
                                                     text="Find Prev.", command=lambda: self.find(-1))
        self.tk['ctrl']['find_button_2'].pack(**TB33)

        self.build_var_gui()

        self.tk['io'] = dict()
        self.tk['io']['frame'] = tk.Frame(self.core)
        self.tk['io']['frame'].pack(expand=True, **LB33)
        self.tk['io']['scroll_y'] = tk.Scrollbar(self.tk['io']['frame'])
        self.tk['io']['scroll_y'].pack(side=tk.RIGHT, fill=tk.Y)
        self.tk['io']['scroll_x'] = tk.Scrollbar(self.tk['io']['frame'], orient='horizontal')
        self.tk['io']['scroll_x'].pack(side=tk.BOTTOM, fill=tk.X)
        self.tk['io']['tabula'] = tk.Text(self.tk['io']['frame'],
                                          yscrollcommand=self.tk['io']['scroll_y'].set,
                                          xscrollcommand=self.tk['io']['scroll_x'].set,
                                          width=TABULA_WIDTH,
                                          height=40,
                                          wrap=tk.NONE)
        self.tk['io']['tabula'].pack(expand=True, **TB33)
        self.tk['io']['scroll_y'].config(command=self.tk['io']['tabula'].yview)
        self.tk['io']['scroll_x'].config(command=self.tk['io']['tabula'].xview)
        self.reg['io.tabula'] = self.tk['io']['tabula']

    def build_var_gui(self, *_):
        """ Constructor for the DB requests panel---which mutates depending on the request being made

        :return:
        """
        # First Clear the Widget
        for item in self.tk['ctrl']['var_frame'].winfo_children():
            item.destroy()

        # Lookup request type (defines layout)
        req_type = self.reg['ctrl.req_type'].get()

        if not req_type:
            system_log.warning(f"Unknown request type somehow provided '{req_type}'.")
            return

        # Lookup the required inputs for the given request type
        args = REQ_ARGS.get(req_type, None)

        if args is None:
            system_log.warning(f"Unknown request type somehow provided '{req_type}'.")
            return

        # Manage the change in what tkinter variables will be required to manage the inputs
        if self.var_args:
            for k in list(self.var_args.keys()):
                if k in ['previous_document', 'new_document', 'incoming_dict']:
                    continue
                if hasattr(self.var_args[k], 'destroy'):
                    self.var_args[k].destroy()
                else:
                    del self.var_args[k]

        # Dynamically create UI
        vf = self.tk['ctrl']['var_frame']
        for arg in args:
            tk.Label(vf, text=f"{arg}:").pack(**TB33)
            field = VALID_VALUES.get(arg, None)
            if isinstance(field, list):
                self.var_args[arg] = tk.StringVar(vf)
                tk.OptionMenu(vf, self.var_args[arg], *VALID_VALUES[arg]).pack(**TB33)
            elif field == "Entry":
                self.var_args[arg] = tk.Entry(vf)
                self.var_args[arg].pack(**TB33)
            elif field == "Text":
                self.tk['ctrl'][f'{arg}_load'] = tk.Button(vf,
                                                           text="Load from Editor",
                                                           command=lambda x=arg: self.load(x))
                self.tk['ctrl'][f'{arg}_load'].pack(**TB33)
                if self.var_args.get(arg, None):
                    self.tk['ctrl'][f'{arg}_load'].config(fg='green3')
                    self.tk['ctrl'][f'{arg}_load'].config(text='Reload from Editor')
                else:
                    self.tk['ctrl'][f'{arg}_load'].config(fg='red3')
                    self.tk['ctrl'][f'{arg}_load'].config(text='Load from Editor')
            else:
                system_log.warning(f"Unknown field '{field}' for argument '{arg}'.")
                return

        # For 'query_document' add a helper that shows what is present in the collection
        if req_type == 'query_document':
            vf = self.tk['ctrl']['var_frame']
            tk.Button(vf, text="What's here?", command=self.list_options).pack(**TB33)

    def load(self, arg, *_):
        """ Saves the contents of the IO widget to a register

        The color of the button is updated to indicate if the register is empty (red) or occupied (green)

        :param arg: Name of the variable being occupied
        :return:
        """
        self.var_args[arg] = yaml.safe_load(self.reg['io.tabula'].get("1.0", tk.END))
        if self.var_args[arg]:
            self.tk['ctrl'][f'{arg}_load'].config(fg='green3')
            self.tk['ctrl'][f'{arg}_load'].config(text='Reload from Editor')
        else:
            self.tk['ctrl'][f'{arg}_load'].config(fg='red3')
            self.tk['ctrl'][f'{arg}_load'].config(text='Load from Editor')

    def clear_reg(self):
        """ Clears the registers storing information from the IO widget

        :return:
        """
        for doc in ["previous_document", "new_document", "incoming_dict"]:
            if doc in self.var_args:
                del self.var_args[doc]
                self.tk['ctrl'][f'{doc}_load'].config(fg='red3')
                self.tk['ctrl'][f'{doc}_load'].config(text='Load from Editor')

    def reformat(self, *_):
        """ Cleans the contents of the IO widget.

        Attempts to load it as a YAML (then as JSON if YAML fails) and update the IO widget with the clean YAML
        representation.

        :return:
        """
        try:
            text = yaml.safe_load(self.reg['io.tabula'].get("1.0", tk.END).replace("\t", "  "))
        except yaml.YAMLError as ye:
            try:
                text = json.loads(self.reg['io.tabula'].get("1.0", tk.END).replace("\t", "  "))
            except json.decoder.JSONDecodeError as je:
                system_log.warning("Could not parse as either a YAML or JSON")
                print(f"YAML:\n{repr(ye)}")
                print("==== = = = = ==== = = = =")
                print(f"JSON:\n{repr(je)}")
                self.tk['ctrl']['info']['text'] = f'Text could not be parsed'
                return
        self.tk['ctrl']['info']['text'] = "--"
        self.tk['io']['tabula'].delete('1.0', tk.END)
        stream = StringIO()
        yaml.safe_dump(text, stream)
        self.tk['io']['tabula'].insert(tk.END, stream.getvalue())
        del stream

    def submit(self, *_):
        """ Pools inputs to create and post a database request

        Results are printed to the IO widget

        :return:
        """
        req_type = self.reg['ctrl.req_type'].get()

        if not req_type:
            system_log.warning(f"Unknown request type somehow provided '{req_type}'.")
            return

        my_request = {'request_type': req_type}
        for arg in REQ_ARGS[req_type]:
            try:
                if isinstance(self.var_args[arg], dict):
                    value = self.var_args[arg]
                elif hasattr(self.var_args[arg], 'get'):
                    value = self.var_args[arg].get()
                else:
                    value = self.var_args[arg]
                if isinstance(value, str):
                    _s = value[0:1]
                    _e = value[-1:]
                    if _s == "[" and _e == "]":
                        try:
                            value = json.loads(value)
                        except json.decoder.JSONDecodeError:
                            self.tk['ctrl']['info']['text'] = f"Could not parse json based input"
                my_request.update({arg: value})
            except KeyError:
                self.tk['ctrl']['info']['text'] = f"Field '{arg}' not specified"
                return

        for k, v in my_request.items():
            if not v:
                self.tk['ctrl']['info']['text'] = f"Field '{k}' may be missing\n" \
                                                  f"Click submit again if you're sure"
                if self.janus:
                    self.janus = False
                    return

        self.janus = True

        self.tk['io']['tabula'].delete('1.0', tk.END)

        # Note extra arg to database_request: this deactivates the raises and specifies a detailed return value
        doc, err, code = dbi.database_request(self.name, my_request, True)
        if doc == "Connection Failed":
            self.tk['io']['tabula'].insert(tk.END, f"Failed to connect to DB, check connection")
        else:
            if code != 200:
                self.tk['io']['tabula'].insert(tk.END, f"DB Returned Error Code: {code}")
                self.tk['io']['tabula'].insert(tk.END, pformat(err, width=TABULA_WIDTH // 2))
            if doc == "Error":
                self.tk['io']['tabula'].insert(tk.END, f"DB Encountered an Error")
                self.tk['io']['tabula'].insert(tk.END, pformat(err, width=TABULA_WIDTH // 2))
            else:
                stream = StringIO()
                yaml.safe_dump(doc, stream)
                self.tk['io']['tabula'].insert(tk.END, stream.getvalue())
                del stream

        self.tk['ctrl']['info']['text'] = f'{datetime.now().strftime("%I:%M:%S %p")}: {code}'

    def reopen(self):
        """ Pulls contents of the register into the IO window

        :return:
        """
        text = None
        for doc in ["incoming_dict", "new_document", "previous_document"]:
            if doc in self.var_args:
                text = self.var_args[doc]
                break
        if text is None:
            return
        self.tk['ctrl']['info']['text'] = "--"
        self.tk['io']['tabula'].delete('1.0', tk.END)
        stream = StringIO()
        yaml.safe_dump(text, stream)
        self.tk['io']['tabula'].insert(tk.END, stream.getvalue())
        del stream

    def list_options(self):
        """ Performs a search on the Database and attempts to print results to the IO widget

        By default, simply lists the names of the items present in the collection (a what's in this collection)
        functionality.  But can be given a criterion by which to filter the results.

        The criterion is given in the 'search_field' entry and is formatted "key::criterion" where "key" is a flattened
        dictionary key (e.g. "toplevel.level1.level2") and "criterion" is a boolean expression.  The string specified
        in "criterion" will be evaluated using eval() [**NOTE**: Not safe for a public-facing server] and the value of
        the item keyed by "key" can be referenced as "_val_".

        :return:
        """
        # Figure out the collection
        sel_col = self.var_args.get('collection', None)
        meta_search = self.var_args.get('search_field', tk.Entry()).get()
        meta_value = Narg()
        if "::" in meta_search:
            meta_search, meta_value = meta_search.split("::")
        if sel_col is None:
            return
        sel_col = sel_col.get()
        # Pull al the named objects
        if sel_col:
            sel_col_doc, _, _ = dbi.query_collection('MC', sel_col)
            names = list()
            for _, v in sel_col_doc.items():
                name = v.get(DBG_CONTAINER_NAME, v.get(Q_NAME, None))
                if meta_search:
                    try:
                        _val_ = functionals.dictionary_direct_access(v, meta_search)
                        if not isinstance(meta_value, Narg):
                            try:
                                if not eval(f"{meta_value}"):
                                    continue
                            except Exception as e:
                                name = [repr(e), "Search criteria are specified by using '::' in the search_field "
                                                 "entry and expect a boolean-valued expression determined by eval() "
                                                 "and use _val_ to reference the value found by the keypath lookup"]
                    except KeyError:
                        continue
                if name is None:
                    continue
                names.append(name)
        else:
            names = ["Please select a collection from the dropdown menu to view a list of document names", ]
        if (not names) and meta_search:
            names = [f"Could not find any documents with a valid keypath of '{meta_search}' using '.' as "
                     f"a delimiter.  Clear 'search_field' to search for all documents or revise the key path.\n\n"
                     f"You may append '::' followed by a python expression to add a criterion to the search.  "
                     f"If referencing the value specified by the keypath, use '_val_'.  The expression should be "
                     f"a boolean statement that is parsable by eval()."]
        names.sort()
        names = "\n".join(names)
        # Printout
        self.tk['io']['tabula'].delete('1.0', tk.END)
        self.tk['io']['tabula'].insert(tk.END, names)

    def find(self, step, *_):
        """ A basic text-search function (wraps)

        Jumps the cursor of the IO widget to an instance of the searched term

        :param step: Indicates direction (find next [+1] vs find previous [-1])
        :return:
        """
        word = self.tk['ctrl']['find_field'].get()
        self.tk['io']['tabula'].tag_delete("search")
        index = "1.0"

        while True:
            if word:
                _start = self.tk['io']['tabula'].search(word, index, stopindex=tk.END)
                if not _start:
                    break
                self.finds.append(_start)
                _row, _col = _start.split('.')
                index = f"{_row}.{int(_col) + len(word)}"
            else:
                break
        if not self.finds:
            return

        if step < 0:
            self.cursor += step
        selection = self.finds[self.cursor % len(self.finds)]
        if step > 0:
            self.cursor += step

        _row, _col = selection.split('.')
        _end = f"{_row}.{int(_col) + len(word)}"
        self.tk['io']['tabula'].see(selection)
        self.tk['io']['tabula'].tag_add("search", selection, _end)
        self.tk['io']['tabula'].tag_configure("search", background="skyblue")

    def _clear_tag(self, *_):
        """ Clears the tags generated by the find() method and resets the cursor position

        :return:
        """
        self.cursor = 0
        self.finds = []
        self.tk['io']['tabula'].tag_delete("search")


if __name__ == '__main__':
    core = tk.Tk()

    # wellplate_coll, _, _ = dbi.query_collection('MC', COL_WELLPLATES)
    # for _, v in wellplate_coll.items():
    #     name = v.get(DBG_CONTAINER_NAME, None)
    #     contents = v.get(DBG_CONTENTS, {})
    #     if isinstance(contents, dict):
    #         size = len(contents.keys())
    #     elif isinstance(contents, str):
    #         size = 0
    #     else:
    #         size = -1
    #     if name:
    #         print(f"{name} ({size})")

    my_manager = DatabaseManager(core, "MC")

    core.mainloop()
