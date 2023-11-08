"""
Default keyword arguments for message passing

@author: Ben C
"""

# keyword arguments carry the original message
ORIGINAL_MSG = 'original_message'  # only used in Control/functionals.py

# Keywords for operation's kwargs
KWA_FAULT_KEY = "fault_keys"             # Not used
KWA_FAULT_VAL = "fault_value"            # Used in Communication/message.py & Control/functionals.py
KWA_TASK_KEY = "task_key"                # Used in Communication/message.py & Control/functionals.py
KWA_TASK_RESULT = "result"               # Used in Communication/message.py & Control/functionals.py
KWA_NET_MEMBERS = "updated_member_list"  # Used in Communication/message.py & Control/functionals.py
KWA_STS_REPORT = 'status_report'         # Used in Communication/message.py & Control/functionals.py
KWA_FILENAME = 'filename'                # Not used
KWA_REQUESTER = 'requester'              # Used in Communication/message.py & Control/functionals.py
