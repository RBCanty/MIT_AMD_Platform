""" A helper module to allow a tkinter GUI to be defined from a dictionary (collapsable in IDE) with the ability to
reference components by Name.

Prototype of the tkfire module on github, do not confuse the two as they operate differently.
Use this version for this project.

@author: Ben C
"""

import tkinter as tk
from tkinter import scrolledtext as st
from tkinter import ttk
from gui_constants import *


class Dispatcher:
    """ Defines the Method Resolution Order for TkFire: tk > st > ttk > custom """
    def __init__(self, custom_elements=None):
        """ Creates a Dispatcher which will be provided with keywords and then determine which
        module/object to call when constructing them.

        **Ordering**: tkinter.tk, tkinter.scrolledtext, tkinter.ttk, custom_elements --> raise ModuleNotFoundError

        The ordering can be changed by overwriting the __getattr__(self, item) method of the Dispatcher via a
        monkey-patch::

        | if item in dir(tk):
        |     return getattr(tk, item)
        | elif item in dir(st):
        |     return getattr(st, item)
        | elif item in dir(ttk):
        |     return getattr(ttk, item)
        | elif item in self.custom_elements:
        |     return self.custom_elements[item]
        | else:
        |     raise ModuleNotFoundError(item)

        References: :class:`TkFire`

        :param custom_elements: (Optional) a dict of user-made classes (keyed by class name as
          referenced in TkFire's mother)
        """
        if custom_elements is None:
            custom_elements = dict()
        self.custom_elements = custom_elements

    def __getattr__(self, item):
        if item in dir(tk):
            return getattr(tk, item)
        elif item in dir(st):
            return getattr(st, item)
        elif item in dir(ttk):
            return getattr(ttk, item)
        elif item in self.custom_elements:
            return self.custom_elements[item]
        else:
            raise ModuleNotFoundError(item)


LAYOUT = 'LAYOUT'
TYPE = 'TYPE'
CHILDREN = 'CHILDREN'
# __gui_constants__________________________________________________
# BOTH33 = {'fill': 'both', 'ipadx': 3, 'ipady': 3}
# TB33 = {'side': 'top', 'fill': 'both', 'ipadx': 3, 'ipady': 3}
# LB33 = {'side': 'left', 'fill': 'both', 'ipadx': 3, 'ipady': 3}
# T33 = {'side': 'top', 'ipadx': 3, 'ipady': 3}
# L33 = {'side': 'left', 'ipadx': 3, 'ipady': 3}
# TX33 = {'side': 'top', 'fill': 'x', 'ipadx': 3, 'ipady': 3}
# LX33 = {'side': 'left', 'fill': 'x', 'ipadx': 3, 'ipady': 3}
# TY33 = {'side': 'top', 'fill': 'y', 'ipadx': 3, 'ipady': 3}
# LY33 = {'side': 'left', 'fill': 'y', 'ipadx': 3, 'ipady': 3}
# P33 = {'ipadx': 3, 'ipady': 3}
# PACK_SCROLL = {'side': 'right', 'fill': 'y'}
# ST_CONTENTS = ["1.0", 'end']
# INPUT_CONTENTS = [0, "end"]
# _________________________________________________________________


def grid_arg(row, col, kwargs=None):
    """ Used to abridge the specification of the tkinter grid packer

    :param row: the row number (0-indexed)
    :param col: the column number (0-indexed)
    :param kwargs: And additional kwargs for pack (like columnspan)
    :return: A dict of keyword arguments for grid()
    """
    if kwargs is None:
        kwargs = {}
    args = {'row': row, 'column': col}
    args.update(kwargs)
    return args


class TkFire:
    """ A Wrapper for tkinter for a nicer experience in the IDE

    Elements of the GUI are of the form {name: {layout: [], type: [], children: {}}}

    - name will define how the element is references in the self.gui dictionary
    - layout keys to a list whose first element is the layout method (e.g. pack, grid) and whose second element is a
      dictionary of the kwargs which are passed to the layout method.
    - type keys to a list whose first element is the name of the Widget (e.g. LabelFrame, Button) and whose second
      element is a dictionary of the kwargs given to the widget's constructor
    - children is an optional element which keys to a dictionary with more TkFire Elements.

    This nestable definition of elements is passed into the constructor as the 'mother'.  Elements of the 'mother' can
    reference the elements of a companion dictionary 'memory' which can store tkinter variables and functions to bind
    to things like buttons.

    Technically, any element keyed by a string which maps to a list will be interpreted as [constructor, {arguments}],
    and so tkinter variables can be specified in the constructor of other Elements (see Opt1 in the example below).

    Elements outside TkFire can be accessed via the self.gui dictionary where the path is bang-separated.

    Example::

    | # Defining the memory
    | memory = {
    |     'cmd_continue': lambda: print("Continuing..."),
    |     'var_opt_1': tk.IntVar,
    | }
    |
    | # Defining the mother
    | mother = {
    |     'main_panel': {
    |         layout: ['pack', {'side': 'left', 'fill': 'both', 'ipadx': 3, 'ipady': 3}],
    |         type: ['LabelFrame', {'text': "Left Side", 'width': 16}],
    |         children: {
    |             'Button_01': {
    |                 type: ['Button', {'text': "Play Again", 'command': 'cmd_continue'}],
    |                 layout: ['grid', {'row': 0, 'column': 0}],
    |             },
    |             'Button_02': {
    |                 type: ['Button', {'text': "Quit"}],
    |                 layout: ['grid', {'row': 0, 'column': 0}],
    |             },
    |             'Opt1': {
    |                 TYPE: ['OptionMenu', {'variable': ['var_opt_1', {'value': 0}], 'values': 'options'}],
    |                 LAYOUT: ['pack', TB33],
    |             },
    |         },
    |     ...
    |     }
    |
    | my_core = tk.Tk()
    |
    | # Construct the TkFire object
    | my_ui = TkFire(my_core, mother, memory)
    |
    | # Bind a command to a button after construction
    | my_ui.gui['main_panel!Button_02']['command'] = lambda *_: print("Exiting...")
    | ...
    | my_core.mainloop()

    """
    def __init__(self, core, memory: dict, mother: dict, dispatcher=Dispatcher()):
        """ Creates a TkFire object which allows a tkinter GUI to be created and modified using dictionary syntax.

        References: :class:`Dispatcher`

        :param core: a tkinter root
        :param memory: a dictionary of variables (keys may be referenced in mother)
        :param mother: a dictionary specifying the GUI
        :param dispatcher: Determines how the objects references in mother will be constructed
        """
        self.core = core
        self.memory = memory
        self.mother = mother
        self._dispatcher = dispatcher
        self.gui = dict()
        self._build("", self.mother)

    def _build(self, parent, structure):
        """ Constructs the dictionary which contains all the widgets (recursively)

        commits to self.gui

        :param parent: The parent widget
        :param structure: a (sub) dictionary of the mother
        :return: None
        """
        if parent:
            root = self.gui[parent]
            path = parent + "!"
        else:
            root = self.core
            path = ""
        for k, child in structure.items():
            # index = parent+k+"!"+child[TYPE][0]
            index = path + k
            # print("Parent: ", parent)
            # print("Child: ", child)
            # print("Child Type: ", child[TYPE])
            kwargs = child[TYPE][1]
            args = list()
            has_scrolly = False
            has_scrollx = False
            if 'variable' in kwargs:
                self.memory[kwargs['variable'][0]] = self.memory[kwargs['variable'][0]](root, **kwargs['variable'][1])
                args.append(self.memory[kwargs['variable'][0]])
                del child[TYPE][1]['variable']
            if 'command' in kwargs:
                child[TYPE][1]['command'] = self.memory[kwargs['command']]
            if 'values' in kwargs:
                args += self.memory[kwargs['values']]
                del child[TYPE][1]['values']
            if 'scrolly' in kwargs:
                has_scrolly = True
                del child[TYPE][1]['scrolly']
            if 'scrollx' in kwargs:
                has_scrollx = True
                del child[TYPE][1]['scrollx']
            self.gui[index] = getattr(self._dispatcher, child[TYPE][0])(root, *args, **child[TYPE][1])
            getattr(self.gui[index], child[LAYOUT][0])(**child[LAYOUT][1])
            if has_scrolly:
                self.gui[index + '!Scrolly'] = tk.Scrollbar(root)
                self.gui[index]['yscrollcommand'] = self.gui[index + '!Scrolly'].set
                self.gui[index + '!Scrolly'].config(command=self.gui[index].yview)
                self.gui[index + '!Scrolly'].pack(**PACK_SCROLL)
            if has_scrollx:
                self.gui[index + '!Scrollx'] = tk.Scrollbar(root, orient=tk.HORIZONTAL)
                self.gui[index]['xscrollcommand'] = self.gui[index + '!Scrollx'].set
                self.gui[index + '!Scrollx'].config(command=self.gui[index].xview)
                self.gui[index + '!Scrollx'].pack(**{'side': tk.BOTTOM, 'fill': tk.X})
            if CHILDREN in child.keys():
                self._build(index, child[CHILDREN])

    def destroy(self):
        """ Calls destroy() on the tkinter core

        :return: None
        """
        self.core.destroy()


if __name__ == '__main__':
    from pprint import pprint

    my_core = tk.Tk()
    my_memory = dict()
    my_memory['Opt1'] = tk.IntVar
    my_memory['Opt2'] = tk.IntVar
    my_memory['Opt3'] = tk.IntVar
    my_memory['options'] = [1, 2, 3, 4, 5]
    my_memory['btn_cmd'] = lambda: print("Hello")
    my_mother = {
        'left_panel': {
            LAYOUT: ['pack', LB33],
            TYPE: ["LabelFrame", {'text': "Left Side", 'width': 16}],
            CHILDREN: {
                'Button1': {
                    TYPE: ['Button', {'text': "Button 1", 'command': 'btn_cmd'}],
                    LAYOUT: ['grid', {'row': 0, 'column': 0}],
                },
                'Button2': {
                    TYPE: ['Button', {'text': "Button 2", 'command': 'btn_cmd'}],
                    LAYOUT: ['grid', {'row': 1, 'column': 0}],
                },
                'Button3': {
                    TYPE: ['Button', {'text': "Button 3", 'command': 'btn_cmd'}],
                    LAYOUT: ['grid', {'row': 2, 'column': 0}],
                },
                'scrolled_text': {
                    TYPE: ["ScrolledText", {'wrap': tk.WORD, 'width': 16, 'height': 12}],
                    LAYOUT: ['grid', {'row': 3, 'column': 1}],
                    CHILDREN: {
                        'my_entry': {
                            TYPE: ['Entry', {}],
                            LAYOUT: ['pack', LB33]
                        }
                    }
                }
            },
        },
        'right_panel': {
            'LAYOUT': ['pack', LB33],
            'TYPE': ["LabelFrame", {'text': "Right Side", 'width': 16}],
            'CHILDREN': {
                'Opt1': {
                    TYPE: ['OptionMenu', {'variable': ['Opt1', {'value': 0}], 'values': 'options'}],
                    LAYOUT: ['pack', TB33],
                },
                'Opt2': {
                    TYPE: ['OptionMenu', {'variable': ['Opt2', {'value': 0}], 'values': 'options'}],
                    LAYOUT: ['pack', TB33],
                },
                'Opt3': {
                    TYPE: ['OptionMenu', {'variable': ['Opt3', {'value': 0}], 'values': 'options'}],
                    LAYOUT: ['pack', TB33],
                },
            },
        },
    }

    my_gui = TkFire(my_core, my_memory, my_mother)
    my_core.mainloop()

    pprint(my_gui.gui)
