""" Collection of constants used in the MCN

* Names of the data pipes
* MM/DD/YYYY and timeout constants
* Wellplate well names
* Library of nicknames for the systems

@author: Ben C
"""

# Constants for Initialization
# # Data Pipes
MSG_Q_IN = "msg_queue_in"
MSG_Q_OUT = "msg_queue_out"
ERROR_Q = "error_queue"
CHILD_COM = "child-grandchild_communication_pipe"
INTERNAL = 'internal_pipe'
FOREMAN = 'analytics_manager'

# Constants for internal operations
TIMEOUT = 30  # TODO: Create a library of standard wait times (0.1, 1.0, 5.0, 30.0, 60*60, 8*60*60)
DATE_FORMAT = "%m/%d/%Y"

# Database Queue keywords for the operation (and those which Checklist class adds upon loading them into the MCN)
# DBQ_PLATE_ID = 'plate_id'
# DBQ_STEP_NUM = 'step'

# Wellplates
WELLPLATE_COLUMNS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
WELLPLATE_ROWS = ["A", "B", "C", "D", "E", "F", "G", "H"]
WELL_IDS = [f"{L}{N}" for L in WELLPLATE_ROWS for N in WELLPLATE_COLUMNS]
WELL_IDS_FLIPPED = [f"{L}{N}" for L in WELLPLATE_ROWS.__reversed__() for N in WELLPLATE_COLUMNS.__reversed__()]


def get_flipped_well_id(well_id):
    try:
        return WELL_IDS_FLIPPED[WELL_IDS.index(well_id)]
    except ValueError:
        return None


# Names/Roles
NN = {"MC": "Master Controller",
      "LC": "Liquid Chromatography",
      "AH": "Automation Handler",
      "SP": "Special Processes",
      "NM": "NMR Machine",
      "Ra": "Robotic Arm",
      "Lh": "Liquid Handler",
      "Pr": "Plate Reader",
      "Sw": "Shimadzu Wrapper",
      "as": "Autosamper",
      "fc": "Fraction Collector",
      "Fs": "FTIR Shell",
      "Ss": "Storage Shell",
      "Th": "Thermoreactor",
      }


def banner(name):
    """ Purely aesthetic and for fun """
    try:
        ret_val = f"+{'-'*58}+\n" \
                  f"| Greetings, this is {NN[name].ljust(37, ' ')} |\n" \
                  f"+{'-'*58}+\n"
    except:  # noqa
        ret_val = ""
    return ret_val
