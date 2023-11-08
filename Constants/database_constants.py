"""Database Queue keywords for the operation (and those which Checklist class adds upon loading them into the MCN)

Contains many constants (prefix legend):

* DB - Lists of keywords for a given DB context
* DBG - Database General (keywords used in multiple contexts)
* COL - Collections
* Q - Queue
* DBQ - List of keywords for a given DB.Queue context
* QOP - Queue Operations
* QDET - Queue Operation Detail
* W - Wellplate
* DBW - List of keywords for a given DB.Wellplate context
* R - Reactants
* DBR - List of keywords for a given DB.Reagent context
* CON - Consumables
* S - Solvent
* DBS - List of keywords for a given DB.Solvent context
* DBQS - Database Queue Status

Contains some mapping for request generation:

* REQ_TYPES
* REQ_ARGS
* VALID_VALUES
* OPERATIONS_DETAIL_REQUIREMENTS

"""

DBG_ID = '_id'
DBG_CONTAINER_NAME = 'container_name'
DBG_CONTAINER_PLATETYPE = 'plate_type'
DBG_CONTENTS = 'contents'
DBG_CONFIRM = 'confirmation'
DBG_REAGENTS = 'reagents'
DBG_SOLVENT = 'solvent'
DBG_TARGET = 'target_product'
DBG_TOTAL_VOLUME = 'total_volume'
DBG_DATE_CREATED = 'date_created'
DBG_DATE_UPDATED = 'date_updated'
DBG_LABWARE_TYPE = 'labware_type'
DBG_LOCATION = 'location'
DBG_PLATE_WELL = 'plate_well'
DBG_CHEM_NAME = 'chemical_name'
DBG_CHEM_SMILES = 'chemical_smiles'
DBG_INT_STD = 'internal_standard'
DBG_INT_STD_CONC = 'internal_standard_concentration'
DBG_VOLUME_UL = 'volume_ul'

COL_REAGENTS = 'reagents'
COL_WELLPLATES = 'wellplates'
COL_QUEUE = 'queue'
COL_CONSUMABLES = 'consumables'
COL_SOLVENTS = 'solvents'
COL_HISTORICAL = 'historical'
DB_COLLECTIONS = [COL_REAGENTS, COL_WELLPLATES, COL_QUEUE, COL_CONSUMABLES, COL_SOLVENTS, COL_HISTORICAL]

Q_CONTAINERS = 'containers'
Q_OPERATIONS_LIST = 'operations'
Q_NAME = 'queue_name'
Q_STATUS = 'status'
Q_DEPENDENCY = 'dependency'
Q_RECORD = 'fault_record'
DB_QUEUE = [DBG_ID, Q_CONTAINERS, DBG_DATE_CREATED, Q_OPERATIONS_LIST, Q_NAME, Q_STATUS, Q_DEPENDENCY]

DBQ_CONTAINERS = [DBG_CONTAINER_NAME, DBG_CONTAINER_PLATETYPE, DBG_CONTENTS]
DBQ_CONTENTS = [DBG_CONFIRM, DBG_PLATE_WELL, DBG_REAGENTS, DBG_SOLVENT, DBG_TARGET, DBG_TOTAL_VOLUME]

QOP_AGENT = 'agent'
QOP_COMPLETED = 'completed'
QOP_CONTAINER = 'container'
QOP_DETAILS = 'details'
QOP_OPERATION = 'operation'
QOP_TIME = 'time_est'
QOP_START = 'start_time'
QOP_END = 'end_time'
DBQ_OPERATIONS = [QOP_AGENT, QOP_COMPLETED, QOP_CONTAINER, QOP_DETAILS, QOP_OPERATION, QOP_TIME, QOP_START, QOP_END]

QDET_PAIRED = 'is_paired'
QDET_TIME = 'schedule_time'

W_CONT_CATEGORY = 'container_category'
W_BARCODE = 'barcode'
DB_WELLPLATE = [DBG_ID, W_CONT_CATEGORY, DBG_CONTAINER_NAME, DBG_CONTENTS, DBG_DATE_CREATED, DBG_DATE_UPDATED, DBG_LABWARE_TYPE, DBG_LOCATION, W_BARCODE]

WCONT_POTENTIAL_VOLUME = 'potential_volume'
WCONT_INT_STD_VOL = 'internal_standard_volume'
DBW_CONTENTS = [DBG_CONFIRM, DBG_PLATE_WELL, WCONT_POTENTIAL_VOLUME, DBG_REAGENTS, DBG_SOLVENT, DBG_TARGET, DBG_TOTAL_VOLUME, DBG_INT_STD, DBG_INT_STD_CONC, WCONT_INT_STD_VOL]

DB_REAGENTS = [DBG_ID, DBG_CONTAINER_NAME, DBG_CONTENTS, DBG_DATE_CREATED, DBG_DATE_UPDATED, DBG_LABWARE_TYPE, DBG_LOCATION]

R_CONCENTRATION = 'concentration_molar'
R_PLATE_ID = 'well_plate_id'
DBR_CONTENTS = [DBG_CHEM_NAME, DBG_CHEM_SMILES, R_CONCENTRATION, DBG_PLATE_WELL, DBG_VOLUME_UL, R_PLATE_ID]

CON_NAME = 'consumable_name'
CON_NUM = 'consumable_number'
CON_TIP_LOCATIONS = 'tip_locations'
CON_STATUS = 'consumable_status'
DB_CONSUMABLES = [DBG_ID, CON_NAME, CON_NUM, DBG_CONTAINER_NAME, DBG_DATE_UPDATED, DBG_LABWARE_TYPE, DBG_LOCATION, CON_TIP_LOCATIONS, CON_STATUS]

DB_SOLVENTS = [DBG_ID, DBG_CONTAINER_NAME, DBG_CONTENTS, DBG_DATE_CREATED, DBG_DATE_UPDATED, DBG_LABWARE_TYPE, DBG_LOCATION]

S_XX = 'XX'
DBS_CONTENTS = [S_XX]

S_DATE = 'chemical_date'
S_CONTAINER_ID = 'container_id'
DBS_CONTENTS_XX = [S_DATE, DBG_CHEM_NAME, DBG_CHEM_SMILES, S_CONTAINER_ID, DBG_INT_STD, DBG_INT_STD_CONC, DBG_VOLUME_UL]

REQ_TYPES = ["query_collection", "query_document", "query_platform_location",
             "update_database_document", "update_database_field",
             "add_to_database",
             "move_database_document"]
REQ_ARGS = {"query_collection": ['collection'],
            "query_document": ['collection', 'search_field', 'search_term'],
            "query_platform_location": ['location_tag', 'sublocation'],
            "update_database_document": ['collection', 'previous_document', 'new_document'],
            "update_database_field": ['document_name', 'update_field', 'old_value', 'new_value'],
            "add_to_database": ['collection', 'incoming_dict'],
            "move_database_document": ['doc_id', 'current_collection', 'new_collection']}
VALID_VALUES = {"collection": DB_COLLECTIONS,
                "current_collection": DB_COLLECTIONS,
                "new_collection": DB_COLLECTIONS,
                "search_term": "Entry",
                "search_field": "Entry",
                "location_tag": "Entry",
                "sublocation": "Entry",
                "doc_id": "Entry",
                "previous_document": "Text",
                "new_document": "Text",
                "incoming_dict": "Text",
                "document_name": "Entry",
                "update_field": "Entry",
                "old_value": "Entry",
                "new_value": "Entry"
                }

# Native
NO_DEPENDENCIES = 'none'
DBQ_STATUS = 'status'
DBQS_IDLE = 'idle'
DBQS_WORK = 'in_progress'
DBQS_FAIL = 'error'
DBQS_DONE = 'complete'
DBQ_DATE_FORMAT = "%m/%d/%Y"  # e.g. "09/24/2021"
DB_YES = 'yes'
DB_NO = 'no'

# Added
DBQ_STEP_NUM = 'step'  # Database Queue.containers keywords
DBQ_PLATE_ID = 'plate_id'

# Operations Details
"""
operation_name: {
    'required': [list of each required detail keyword],
    'dependant': {'detail keyword': [list of detail keywords required if its key is included]},
    'optional': [list of each optional detail keyword -- Not implemented]
}
"""
OPERATIONS_DETAIL_REQUIREMENTS = {
    None: {},
    'db_print': {'required': ['field', 'term'],
                 'dependant': {},
                 'optional': []},
    'run_semiprep_batch': {'required': ['target_container', ],
                           'dependant': {},
                           'optional': []},
    'transfer_wellplate': {'required': ['target_destination'],
                           'dependant': {},
                           'optional': ['disassemble_inert', 'temperature']},
    'start_stop_heater_shaker': {'required': ['power'],
                                 'dependant': {'power': ['rpms', 'temperature']},
                                 'optional': ['rpms', 'temperature']},
    'run_spark': {'required': ['method_type', ],
                  'dependant': {},
                  'optional': ['file_path', ]},
    'run_spark_method': {'required': ['method_path', ],
                         'dependant': {},
                         'optional': ['file_path', ]},
    'move_wellplate': {'required': ['target_destination', ],
                       'dependant': {},
                       'optional': ['scan_barcode', ]},
    'Ss.request': {'required': ['plate_type', ],
                   'dependant': {},
                   'optional': ['plate_name', ]},
    'Ss.stow': {'required': ['plate_name', ],
                'dependant': {},
                'optional': ['plate_type', ]},
    'Th.run': {'required': ['heating_profile', ],
                'dependant': {},
                'optional': []},
}


if __name__ == '__main__':
    #print(dbq)
    pass
