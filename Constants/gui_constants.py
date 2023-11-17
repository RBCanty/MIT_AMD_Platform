""" Constants for use with tkinter

Arguments are designed to be used via unpacking (e.g. tk.object.pack(**BOTH33))

**Packing**: General order (Side, Fill, Internal Pad X, Internal Pad Y); exceptions for GRID_B33 and PACK_SCROLL

**Indexing**: Tkinter indexes by character (#) and by line-offset ("{row}.{col}").
"""

BOTH33 = {'fill': 'both', 'ipadx': 3, 'ipady': 3}
TB33 = {'side': 'top', 'fill': 'both', 'ipadx': 3, 'ipady': 3}
TB33E = {'side': 'top', 'fill': 'both', 'ipadx': 3, 'ipady': 3, 'expand': True}
LB33 = {'side': 'left', 'fill': 'both', 'ipadx': 3, 'ipady': 3}
LB33E = {'side': 'left', 'fill': 'both', 'ipadx': 3, 'ipady': 3, 'expand': True}
LBA3 = {'side': 'left', 'fill': 'both', 'ipadx': 3, 'ipady': 3, 'pady': 3, 'padx': 3}
TBA3 = {'side': 'top', 'fill': 'both', 'ipadx': 3, 'ipady': 3, 'pady': 3, 'padx': 3}
T33 = {'side': 'top', 'ipadx': 3, 'ipady': 3}
L33 = {'side': 'left', 'ipadx': 3, 'ipady': 3}
TX33 = {'side': 'top', 'fill': 'x', 'ipadx': 3, 'ipady': 3}
LX33 = {'side': 'left', 'fill': 'x', 'ipadx': 3, 'ipady': 3}
TY33 = {'side': 'top', 'fill': 'y', 'ipadx': 3, 'ipady': 3}
LY33 = {'side': 'left', 'fill': 'y', 'ipadx': 3, 'ipady': 3}
P33 = {'ipadx': 3, 'ipady': 3}
GRID_B33 = {'ipadx': 3, 'ipady': 3, 'sticky': 'nsew'}
PADDING_3 = {'ipadx': 3, 'ipady': 3, 'pady': 3, 'padx': 3}
PACK_SCROLL = {'side': 'right', 'fill': 'y'}

ST_CONTENTS = ["1.0", 'end']
INPUT_CONTENTS = [0, "end"]
