"""
Created on Wed Aug  4 16:52:57 2021

@author: mattmcdonald

# NOTE: Elements of this code had to be redacted before public distribution, items noted with "Redacted:"
"""

import os
import sys
import time
import queue
import copy
import mcn_status
import requests
import serial
from database_constants import *
from mcn_status import V_FATAL as FATAL_LEVEL
from constants import WELL_IDS as VIAL_NAMES
import tkinter as tk
from threading import Thread
from disposable_wrapper import disposable
from threading import Lock
from ui_hplc_volumes import VolumeSet
import synchronization as sync
from mcn_logging_manager import system_log
import custom_exceptions as cexc
import database_interface as dbi
from datetime import datetime
from pprint import pprint
from LC_data_extract import lc_data_analysis, run_optimization, log_e_chemprop_push, log_e_chemprop_pull, calc_logP
from LC_data_DB_interface import add_reaction, add_extinction_coef_prediction, add_HPLC_logD_measurement
import aceso

update_queue = True
collect_all_unique_products = False
print("IN REAL MODE, connecting to instrument...")
# Redacted: Adding the control interface files to the path

# enumerate READY, BUSY (try again soon), PROBLEM (quick human fix), FATAL (needs troubleshooting)
# Redacted: Which status codes map two which Fault level...
READY = list()
BUSY = list()
PROBLEM = list()
FATAL = list()  # ...end would-be mapping

COLUMN_OPTIONS = ["analytical", "semiprep", "UPLC"]

TEMPLATES_FOLDER = r"C:\Shimadzu\Data\Automated\Templates"
DATA_FOLDER = r"C:\Shimadzu\Data\Automated"
CV_PORT = 'COM3'

# First pass at getting a working shell together... error handling is next on list
class ShimadzuController:
    def __init__(self, pipes):
        self._pipes = pipes
        self.lock = Lock()
        # Redacted: Create instances of the controller and...
        self.Shimadzu_obj = None  # ...bind to "Shimadzu object"
        print("Connected to instrument")
        # info required to keep track of individual batches/plates
        self.samples = queue.Queue()
        self.currentSample = ("LS_initialized", -1, "A0", "SMILES")  # (filename, mol weight, well, SMILES)
        self.running = (False, "Initialized")
        # valid values: "Initialized", "Batch running", "Batch Stopped", "Batch Ended",
        #   "Instrument Error", "Offline", and "Resuming"
        self.folderName = ""
        # self.FRC_next = [0,0]  # this is the next spot the FRC will collect into: [plate, well]
        # copies of the queue, specific operation, wellplate, and fraction plate
        self.queue = {}
        self.wellplate = {}
        self.running = (False, "Initialized")
        self.operation = {}
        self.opn_num = ""
        self.FRC_plate = []
        self.Error_status = ""
        self.chemprop_complete = False
        self.lc_column = "none"

        # add functions to the listener for events fired by Shimadzu... e.g. to commence data analysis
        # Redacted: Binding of our methods to the Shimadzu event listeners

        # track info that could require user intervention
        self.solvent_levels = {'Aqueous': 0.0, 'Organic': 0.0,
                               'MS_makeup': 0.0, 'last_update_time': -1.0,
                               'last_sample': "A0"}

        # Stable UI to handle volumes specified by user
        volumes = VolumeSet(tk.Tk(), title="Solvent Startup Check", dialog="Enter solvent volumes in mL",
                            categories=['Aqueous', 'Organic', 'MS_makeup'], defaults=[1000, 1000, 1000]).run()
        self.solvent_levels['Aqueous'] = volumes[0]
        self.solvent_levels['Organic'] = volumes[1]
        self.solvent_levels['MS_makeup'] = volumes[2]

        system_log.info(f"Initialized solvent volumes to {self.solvent_levels}")

        self.valve = serial.Serial(CV_PORT, baudrate=9600, bytesize=8, parity='N', stopbits=1, timeout=1)
        self.valve.write(b'CP\r\n')
        valve_position = str(self.valve.readline())
        print(valve_position)
        if 'CP' not in valve_position:
            print(f"Problem connecting to column switching valve, check that it is connected to port {CV_PORT}")
            sys.exit()
        if "A" not in valve_position:
            init_column = "Semiprep"
        elif "B" not in valve_position:
            init_column = "UPLC"
        else:
            init_column = "ERROR"
        system_log.info(f"Initial column is {init_column}")
        self.lc_mode = init_column

    def __del__(self):
        self.batch_stop()
        self.disconnect()

    # Needs to notify the MCN that the batch has started successfully (? or does it ?)
    def batch_start(self):
        self.running = (True, "Batch running")
        self.update_solvent_levels(pumps_off=False)
        system_log.info("Batch Started")

    # Notify the MCN that the batch has ended and whether or not it was successful
    def batch_end(self):
        system_log.info("Batch Ended")
        if not self.running[0]:
            return True, "Received 2 or more batch end signals"
        else:
            self.running = (False, "Batch Ended")
            self.update_solvent_levels(pumps_off=True)
            # at the end of a semiprep run need to populate the FRC plate
            if self.lc_column == "semiprep":
                self.frc_plate_populate()
            elif self.lc_column == "analytical":
                if 'reagent' not in self.operation['container']:
                    pass
                    # self.estimate_yield()
            elif self.lc_column == "optimization":
                self.optimization_queue_populate()
            elif self.lc_column == "uplc":
                system_log.info("UPLC Batch has ended")
            # reset key variables
            self.queue = {}
            self.wellplate = {}
            self.operation = {}
            self.opn_num = ""
            self.FRC_plate = []
            self.folderName = ""
            self.chemprop_complete = False

            if self.samples.empty():
                return True, f"{self.lc_column} batch completed successfully."
            else:
                return False, f"F_{self.lc_column} batch failed on sample: {self.currentSample[0]}"

    # Needs to notify the MCN that a batch has been stopped (by issuing a stop command)
    def batch_stop(self):
        # Call a shutdown batch immediately to ensure actual stoppage of the pumps
        result_bool, msg = self.lc_stop()
        if result_bool:
            self.running = (False, "Batch Stopped")
            self.update_solvent_levels(pumps_off=True)
            system_log.info("Batch Stopped")
        else:
            system_log.info("Batch stopped before completion. Problem stopping pumps: " + msg)

    # Keep track of the current sample (removes sample from sampleQueue at start of run)
    def run_start(self):
        system_log.info("Run Started")
        self.currentSample = self.samples.get()
        self.update_solvent_levels(pumps_off=False)

    # Starts analysis of a sample file as soon as it is complete
    def run_end(self):
        self.update_solvent_levels(pumps_off=False)
        if self.lc_mode == 'photochem':
            print("No postrun analysis during photochem runs")
            #save_txt_data(self.currentSample[0], 'photochem', self.currentSample[2], self.folderName)
        elif (self.currentSample[1] != -1) and (self.lc_mode == 'analytical'):
            target_sample = self.currentSample
            system_log.info(f"Beginning analysis of {target_sample[0]}")
            if not self.chemprop_complete:
                stat, msg = self.chemprop_pull()
                if stat:
                    self.chemprop_complete = True
                else:
                    print(msg)
            analysis_thread = Thread(target=self.data_analysis, args=(target_sample,), daemon=True)
            analysis_thread.start()
        else:
            system_log.info("Run Ended")

    # honestly not sure what to do with this
    # (does this only appear when a stop command is issued? Different from Stopped?)
    def data_acq_stop(self):
        self.running = (False, "Batch Stopped")
        system_log.info("Data acquisition stopped")

    # Need to parse error codes, determine if error is problem or fatal, and pass message to MCN
    #  errors seem to be sent with a code that the current implementation doesn't handle!!!
    def instrument_error(self, error_code):
        pressure_error_code = None  # Redacted: Code literal
        d2_lamp_error_code = None  # Redacted: Code literal
        if str(error_code) == pressure_error_code:
            system_log.warning("Pressure exceeded maximum pressure!")
        elif str(error_code) == d2_lamp_error_code:
            system_log.warning("D2 Lamp failed to start")
            #TODO: Some logic to try restarting
        else:
            system_log.warning(f"Unknown Instrument error, error code {error_code}")
        self.running = (False, "Instrument error")
        self.update_solvent_levels(pumps_off=True)
        # reload the current sample into the sample queue, will need to create a new batch on recovery from error
        self.samples.put(self.currentSample)

    def update_time_est(self):
        new_est = 380 * (1 + self.samples.qsize())
        try:
            dbi.set_time_est(self.queue['queue_name'], self.opn_num, new_est)
            print(f"Operation time updated from {self.operation['time_est']} to {new_est}")
        except cexc.DatabaseRequestError:
            print("Failed to update the time estimate")
            return False, f"F_Failed to update time_est for {self.operation['operation']}"

    def chemprop_pull(self):
        # return True, "Not doing this today"
        wellplate_logE = log_e_chemprop_pull(self.wellplate, self.folderName)
        for well, content in wellplate_logE['contents'].items():
            target = content['target_product'][0][0]
            campaign_name = self.queue['queue_name'].rsplit('_', 2)[0]
            add_extinction_coef_prediction(target, content['properties'][target]['log_e'], campaign_name)
        try:
            dbi.update_document("LC", 'wellplates', self.wellplate, wellplate_logE)
        except cexc.DatabaseRequestError:
            return (False,
                    "Failed to update wellplates collection with wellplate containing ChemProp predicted logE values")

        if update_queue:
            updated_queue = copy.deepcopy(self.queue)
            updated_queue['containers'][self.operation['container']] = wellplate_logE
            try:
                dbi.update_document("LC", 'queue', self.queue, updated_queue)
            except cexc.DatabaseRequestError:
                return False, "Failed to update queue with wellplate containing ChemProp predicted logE values"
        return True, "Success updating logE values"

    # Interpret MS and PDA data for multiple runs at the end of a batch
    #  NOT CURRENTLY USED, THIS COULD BE USEFUL WHEN REAGENTS NEED TO BE RUN TO LOOK AT A REACTION
    def batch_data_analysis(self, samples):
        success = True
        msg = "Post batch sample analysis completed successfully"
        msgs = []
        up_to_date_plate, *_ = dbi.query_document('LC', 'wellplate', 'container_name', self.wellplate['container_name'])
        for sample in samples:
            smiles = sample.split('__')[0]
            well = sample.split('__')[1]

            target_sample = ("no MS file needed", "placeholder_mw", sample.split('__')[0], sample.split('__')[1])
            status, ret_msg = self.data_analysis(target_sample)
            if not status:
                print(f"Problem analyzing {sample[3]}: {ret_msg}")
                msgs.append(ret_msg)
                success = False
                msg = "Post batch analysis failed on sample(s): {msgs}"
        return success, msg

    # Interpret MS and PDA data at the end of a run
    def data_analysis(self, target_sample):
        reaction_info = lc_data_analysis(self.folderName, self.wellplate, target_sample)

        # going to make some assumptions...
        #  1. the first solvent in the solvents list will always be the sample solvent
        #  2. all of the volatile solvent is gone by the time it gets to the hplc
        # these assumptions can be fixed if we switch to using DMSO with internal standard in it and label as such
        if 'fraction_transferred' in self.wellplate['contents'][target_sample[2]].keys():
            fraction_transferred = self.wellplate['contents'][target_sample[2]]['fraction_transferred']
            solvents = self.wellplate['contents'][target_sample[2]]['solvents']
            remaining_volume = 0
            for solv in solvents:
                if solv[0] in ('dmso', 'dmf', 'water'):
                    remaining_volume += solv[1]
            # TODO! Check this math... it still seems a little weird
            mols_originally_in_reaction = reaction_info['conc'] * remaining_volume / fraction_transferred
            mols_remaining_in_filtrate = mols_originally_in_reaction * (1 - fraction_transferred)
        else:
            mols_originally_in_reaction = reaction_info['conc'] * self.wellplate['contents'][target_sample[2]]['total_volume']
            mols_remaining_in_filtrate = mols_originally_in_reaction

        reagents = self.wellplate['contents'][target_sample[2]]['reagents']
        pred_reagents = self.wellplate['contents'][target_sample[2]]['predicted_reagents']
        if not pred_reagents:
            pred_reagents = []
        pred_catalyst = self.wellplate['contents'][target_sample[2]]['predicted_catalyst']
        if not pred_catalyst:
            pred_catalyst = []
        reactants = []
        limiting_reactant = ['', 1.0, '']
        for reagent in reagents:
            try:
                if (reagent[2] not in pred_reagents) and (reagent[2] not in pred_catalyst):
                    reactants.append(reagent)
                    if reagent[1] < limiting_reactant[1]:
                        limiting_reactant = reagent
            except TypeError:
                print(f"not sure why this is still an issue... reagent = \n\t {reagent}")

        print(f"DEBUGGING YIELD ESTIMATES:\n\t product conc = {reaction_info['conc']}")
        print(f"\t limiting reactant = {limiting_reactant}")
        print(f"\t remaining volume = {remaining_volume}")
        print(f"\t mols remaining in filtrate = {mols_remaining_in_filtrate}")
        rxn_yield = min(1.0, mols_originally_in_reaction / limiting_reactant[1])
        reaction_info['yield'] = rxn_yield
        add_reaction(self.wellplate['contents'][target_sample[2]], reaction_info, self.queue['operations'])
        # TODO: ^ Add argument = self.pipes[INTERNAL].get([mcs.S_AUXILIS, "LC", "RxnDB"], None)
        #  to allow IP:PORT changes w/o restart

        ret_time = reaction_info['retention_time']
        campaign_name = self.queue['queue_name'].rsplit('_', 2)[0]
        add_HPLC_logD_measurement(target_sample[3][0][0], calc_logP(ret_time), campaign_name)

        best_mz = reaction_info['best_m/z']
        updated_queue = self.queue

        print_info = "No product observed."
        for plate in self.operation['details']['update_plates']:
            original_plate, _, _ = dbi.query_document('LC', 'wellplates', 'container_name',
                                                      self.queue['containers'][plate]['container_name'])
            updated_plate = copy.deepcopy(original_plate)
            u_plate_cont = updated_plate.get(DBG_CONTENTS, "Empty")
            if isinstance(u_plate_cont, str):
                u_plate_cont = dict()
            for well, content in u_plate_cont.items():
                # this needs to be updated checked during reaction optimization
                if content['target_product'] == target_sample[3]:
                    if ret_time > 0:
                        print_info = f"Sample observed with MW = {str(best_mz)} at {str(ret_time)} minutes."
                        updated_plate[DBG_CONTENTS][well]['confirmation'] = str(ret_time)
                        updated_plate[DBG_CONTENTS][well]['product_remaining'] = [target_sample[3],
                                                                                  mols_remaining_in_filtrate]
                        # updated_plate[DBG_CONTENTS][well]['target_product'][0][1] = mols_originally_in_reaction
                        if update_queue:
                            updated_queue[Q_CONTAINERS][plate][DBG_CONTENTS][well]['confirmation'] = str(ret_time)
                            updated_queue[Q_CONTAINERS][plate][DBG_CONTENTS][well]['product_remaining'] = [
                                target_sample[3],
                                mols_remaining_in_filtrate]
                            # updated_queue[Q_CONTAINERS][plate][DBG_CONTENTS][well]['target_product'][0][1] = mols_originally_in_reaction
                    else:
                        updated_plate[DBG_CONTENTS][well]['confirmation'] = 'no mass hit'
                        updated_plate[DBG_CONTENTS][well]['product_remaining'] = [target_sample[3], 0.0]
                        if update_queue:
                            updated_queue[Q_CONTAINERS][plate][DBG_CONTENTS][well]['confirmation'] = 'no mass hit'
                            updated_queue[Q_CONTAINERS][plate][DBG_CONTENTS][well]['product_remaining'] = [
                                target_sample[3], 0.0]

            try:
                dbi.update_document("LC", 'wellplates', original_plate, updated_plate)
            except cexc.DatabaseRequestError:
                return False, f"F_Failed to update retention time in {self.wellplate['container_name']}."

        if update_queue:
            try:
                dbi.update_document("LC", 'queue', self.queue, updated_queue)
            except cexc.DatabaseRequestError:
                return False, f"F_Failed to update retention time in {self.queue['queue_name']}."

        if not reaction_info['yield']:
            system_log.info(print_info + " Yield estimate not yet calculated.")
            return (False,
                    "Retention time updated but yield not yet estimated, may be waiting on new reagents to be added")
        system_log.info(print_info + f" Yield estimated to be {reaction_info['yield']}")
        return True, "Retention time updated in wellplates and queue"

    def frc_plate_populate(self):
        # read txt file to populate contents
        frc_report_file = os.path.join(self.folderName, "FRC_Data.txt")
        with open(frc_report_file) as f:
            report_lines = f.readlines()
        if type(self.FRC_plate[0]['contents']) is str:
            current_vial = 0
        else:
            current_vial = len(self.FRC_plate[0]['contents'])
        current_sample = ''
        collected_frac = False
        idx = 0
        vials = [[]]
        vial_contents = [[]]
        for line in report_lines:
            if line[0:8] == '[Header]':
                collected_frac = False
            if line[0:9] == "Data File":
                current_sample = line[15:]
            if collected_frac:
                vials[idx].append(VIAL_NAMES[current_vial])
                vial_contents[idx].append(current_sample)
                current_vial += 1
                if current_vial > 95:
                    idx += 1
                    vials.append([])
                    vial_contents.append([])
                    current_vial = 0
            if line[0:9] == "Fraction#":
                collected_frac = True
        min_acn = 1000.0
        updated_queue = copy.deepcopy(self.queue)
        for j in range(0, len(vials)):
            updated_contents = self.FRC_plate[j]['contents']
            if updated_contents == 'Empty':
                updated_contents = {}
            i = 0
            for sample in vials[j]:
                sample_name = vial_contents[j][i]
                sample_name = sample_name[-8:-5]
                if sample_name[0] == '\\':
                    sample_name = sample_name[1:]
                if sample_name == '-1':
                    # sanity check
                    if sample != 'H12':
                        print(f"FRACTION COLLECTION MATERIAL ASSIGNED TO WRONG WELLS!! {sample} != H12")
                    updated_contents['H12'] = {'confirmation': 'reference_well',
                                               'final_products': [],
                                               'plate_well': 'H12',
                                               'predicted_catalyst': [],
                                               'predicted_reagents': [],
                                               'product_remaining': 0.0,
                                               'reaction_smiles': '',
                                               'reagents': [],
                                               'target_product': [],
                                               'templates': []}
                    acn_conc = min_acn
                else:
                    updated_contents[sample] = self.wellplate['contents'][sample_name]
                    ret_time = float(updated_contents[sample]['confirmation'])
                    acn_conc = max([0.05, min([0.2 * (ret_time - 0.5), 1])]) * 1000
                    min_acn = min(min_acn, acn_conc)
                updated_contents[sample]['solvents'] = [['acetonitrile', acn_conc, 'CC#N'],
                                                        ['water', 1000 - acn_conc, 'O']]
                updated_contents[sample]['date_updated'] = datetime.now().strftime('%m/%d/%Y %H:%M:%S')
                updated_contents[sample]['total_volume'] = 1000
                updated_contents[sample]['plate_well'] = sample
                i += 1
            if 'H12' in updated_contents.keys():
                updated_contents['H12']['solvents'] = [['acetonitrile', min_acn, 'CC#N'],
                                                       ['water', 1000 - min_acn, 'O']]
                updated_contents['H12']['date_updated'] = datetime.now().strftime('%m/%d/%Y %H:%M:%S')
                updated_contents['H12']['total_volume'] = 1000
                updated_contents['H12']['plate_well'] = 'H12'

            updated_plate = copy.deepcopy(self.FRC_plate[j])
            updated_plate['contents'] = updated_contents
            try:
                dbi.update_document('LC', 'wellplates', self.FRC_plate[j], updated_plate)
            except cexc.DatabaseRequestError:
                return False, f"P_failed to update {self.FRC_plate[j]['container_name']} (fraction plate) in DB"
            for q_container_name, q_container_details in self.queue['containers'].items():
                if q_container_details['container_name'] == self.FRC_plate[j]['container_name']:
                    updated_queue['containers'][q_container_name]['contents'] = updated_contents
                    try:
                        updated_queue['containers'][q_container_name]['plate_type'] = updated_plate['plate_type']
                    except KeyError:
                        updated_queue['containers'][q_container_name]['plate_type'] = updated_plate['labware_type']
        try:
            dbi.update_document('LC', 'queue', self.queue, updated_queue)
        except cexc.DatabaseRequestError:
            return False, f"P_failed to update {self.queue['queue_name']} (queue) in DB"

        return True, "Fraction collection plate successfully updated in wellplates and queues"

    # check to see if there are new reagents that will need to be added to the chromatogram library, and if so, add them
    def check_new_reagents(self):
        return False, "Not running reagent standards yet"
        # get most recent version of the queue (just in case someone else has updated it??)
        # up_to_date_queue, _, _ = dbi.query_document('LC', 'queue', 'queue_name', self.queue['queue_name'])
        # new_containers_queue, new_ops = add_reagent_prep_to_queue(up_to_date_queue, self.wellplate)
        #
        # if new_containers_queue:
        #     dbi.update_document('LC', 'queue', up_to_date_queue, new_containers_queue)
        #     dbi.insert_queue_steps(self.queue['queue_name'], self.opn_num, new_ops)
        #     return True, "New reagent standards required. Added operations to generate new standards"
        # return False, "No new reagents needed"

    # Analyze and update an optimization batch and queue
    def optimization_queue_populate(self):
        # get most recent version of the queue (just in case someone else has updated it??)
        up_to_date_queue, _, _ = dbi.query_document('LC', 'queue', 'queue_name', self.queue['queue_name'])
        contents = self.wellplate['contents']
        new_containers_queue, new_ops = run_optimization(up_to_date_queue, contents)

        dbi.update_document('LC', 'queue', up_to_date_queue, new_containers_queue)
        dbi.insert_queue_steps(self.queue['queue_name'], -1, new_ops)

    # TODO! recalculate how solvent levels are tracked... or better yet get them from the instrument
    def update_solvent_levels(self, pumps_off: bool):
        # Update the solvent levels, return true if ok for next sample, return 
        #  false, stop batch, open user dialog, and send fault to MC if solvent low
        if self.Error_status == "low solvent":
            return True, "B_In the process of stopping"

        threshold = 100

        if self.lc_column == 'analytical':
            flowrate = 1.0
            single_run_time = 6
        elif self.lc_column == 'semiprep':
            flowrate = 3.0
            single_run_time = 6
        elif self.lc_column == 'uplc':
            flowrate = 0.5
            single_run_time = 6
        else:
            return False, "F_LC mode not analytical, semiprep, or uplc"

        if float(self.solvent_levels['last_update_time']) < 0:
            # indicates the pumps are not currently running
            time_passed = 0.0
        else:
            time_passed = (time.time() - self.solvent_levels['last_update_time']) / 60

        # if the current sample and the last sample match, no gradient was run
        if self.currentSample[2] == self.solvent_levels['last_sample']:
            aqueous_used = flowrate * (0.95 * time_passed)
            organic_used = flowrate * (0.05 * time_passed)
            ms_makeup_used = 0.2 * time_passed
        else:  # math that turns standard gradient conditions into solvent amounts
            aqueous_used = flowrate * (0.95 * (time_passed - single_run_time) + 2.85)
            organic_used = flowrate * (0.05 * (time_passed - single_run_time) + 3.15)
            ms_makeup_used = time_passed * 0.2
        self.solvent_levels['Aqueous'] -= aqueous_used
        self.solvent_levels['Organic'] -= organic_used
        self.solvent_levels['MS_makeup'] -= ms_makeup_used

        if pumps_off:
            self.solvent_levels['last_update_time'] = -1. * time.time()
        else:
            self.solvent_levels['last_update_time'] = time.time()

        is_low = (self.solvent_levels['Aqueous'] < threshold) or \
                 (self.solvent_levels['Organic'] < threshold) or \
                 (self.solvent_levels['MS_makeup'] < threshold)

        if not is_low:
            return True, f"Greater than {threshold} mL of each solvent remaining"
        else:
            self.Error_status = "low_solvent"
            # report fault and prompt user to input new volume amounts
            sync.synchronized_fault_add(self._pipes, mcn_status.Fault("LC",
                                                                      FATAL_LEVEL,
                                                                      "Need to refill LC solvent",
                                                                      queue=self.queue.get('queue_name', None)))
            _ = self.lc_stop()

            volumes = VolumeSet(tk.Tk(), title="Solvent Startup Check", dialog="Enter solvent volumes in mL",
                                categories=['Aqueous', 'Organic', 'MS_makeup'], defaults=[0, 0, 0]).run()
            self.solvent_levels['Aqueous'] = volumes[0]
            self.solvent_levels['Organic'] = volumes[1]
            self.solvent_levels['MS_makeup'] = volumes[2]
            sync.synchronized_remove_fault(self._pipes,
                                           mcn_status.Fault("LC", FATAL_LEVEL, "Need to refill LC solvent"))

            return True, f"Greater than {threshold} mL of each solvent remaining"

    @staticmethod
    def disconnect():
        # Redacted: Disconnect from the HPLC controller
        return True

    def status_ping(self):
        # Redacted: Pull and translate status of the HPLC...
        status, details = None, "None"  # ...map to "status" (code) and "details" (semantic)
        if status in READY:
            return "Ready", details
        elif status in BUSY:
            return "Busy", details
        elif status in PROBLEM:
            return "Problem", details
        else:
            return "Fatal", details

    def lc_stop(self):
        self.running = (False, "Batch Stopped")
        # Redacted: Halt HPLC run and save result
        #   If the result is not successful return (False, f"F_{additional details}")

        # Start a new method with flow rate 0, no pda, and oven set < RT, etc.
        # Redacted: Attempt to load the batch file from os.path.join(TEMPLATES_FOLDER, "blank_batch_interrupt.lcb")...
        batch = None  # ...bind result to "batch"
        if batch is None:
            return False, "F_Failed to load interrupt template batch file"

        # Redacted: Populate details of "batch", bind the "shutdown.lcm" method to it and the output file "E_STOP.lcd"
        #   Save the results of the interrupted batch to "interrupt.lcb".
        #   If the save failed, return (False, "F_Failed to save batch file as: interrupt.lcb")

        # Load and execute shutdown batch
        time.sleep(10)
        # Redacted: load the batch results, and record return code
        #   If the result is not a successful one, return (False, f"F_{additional details}")
        time.sleep(5)
        # this will cause issues when batch_start, batch_end, run_start and run_end events are fired
        # Redacted: start the batch, and record return code
        #   If the result is not a successful one, return (False, f"F_{additional details}")

        return True, "Batch stopped"

    def update_from_queue(self, queue_id, operation_id):
        if self.Shimadzu_obj is None:
            return False, "F_Failed to Initialize Shimadzu"
        if not self.samples.empty():
            return False, "F_LC sample queue not empty, unfinished samples from previous run"
        stat1, stat2 = self.status_ping()
        # Check that the instrument isn't busy or in an error state
        if stat1 != "Ready":
            return False, "F_" + stat1 + ": " + stat2
        # from library request the queue name that contains the current operation
        query = {'request_type': 'query_document', 'collection': 'queue',
                 'search_field': 'queue_name', 'search_term': queue_id}
        try:
            query_return, _, _ = dbi.database_request('LC', query)
        except cexc.DatabaseRequestError:
            return False, "F_Failed to load " + queue_id + " from library"

        self.queue = query_return
        self.operation = query_return['operations'][operation_id]
        self.opn_num = operation_id
        if self.operation['agent'] != 'LC':
            return False, "F_Operation agent does not match assigned instrument (LC)"
        if self.operation['completed'] == 'yes':
            return False, "P_Operation already marked as completed"
        plate_name = self.operation['container']
        plate_id = query_return['containers'][plate_name]['container_name']
        if plate_id is None:
            return False, "F_Sample plate not defined (or assigned by previous operation)"

        query = {'request_type': 'query_document', 'collection': 'wellplates',
                 'search_field': 'container_name', 'search_term': plate_id}
        try:
            query_return, _, _ = dbi.database_request('LC', query)
        except cexc.DatabaseRequestError:
            return False, "F_Failed to load " + plate_id + " from library"
        self.wellplate = query_return
        if self.wellplate['location'][0] != 'autosampler':
            return False, "P_Sample plate (" + plate_name + ") not in the autosampler"
        if self.operation['operation'] == 'run_analytical_batch':
            self.lc_column = "analytical"
            if "hplc_processing_method" in self.operation['details'].keys():
                if self.operation['details']['hplc_processing_method'] == 'kinetics':
                    self.lc_column = "uplc"
                elif self.operation['details']['hplc_processing_method'] == 'photochemistry':
                    self.lc_column = "uplc"
        elif self.operation['operation'] == 'run_semiprep_batch':
            self.lc_column = "semiprep"
        elif self.operation['operation'] == 'run_optimization_batch':
            self.lc_column = "uplc"
        elif self.operation['operation'] == 'run_uplc_batch':
            self.lc_column = "uplc"
        else:
            self.lc_column = "none"

        # create a new folder (name) that will house all the files associated with the batch 
        #  requires that optimization take place in a separate Q-doc
        self.folderName = os.path.join(DATA_FOLDER, queue_id, 
                                       self.wellplate['container_name'],
                                       self.lc_column + "_data")
        try:
            os.makedirs(self.folderName)
        except FileExistsError:
            self.folderName = self.folderName + "_" + datetime.now().strftime('%d-%m-%Y_%H-%M-%S')
            os.makedirs(self.folderName)
            print('File already exists, making new file, appended with current time')
        batch_file_name = os.path.join(self.folderName, self.lc_column + "_batch.lcb")
        if os.path.isfile(batch_file_name):
            return False, "P_There is already LC data associated with queue: " + queue_id

        # check for new reagents
        self.check_new_reagents()

        # calculate log_e where not already catalogued
        if self.lc_column != 'uplc':
            log_e_chemprop_push(self.wellplate, self.folderName)

        return True, batch_file_name

    def run_batch(self, batch_dict, **kwargs):
        queue_id = self.queue['queue_name']
        # Redacted: Load the batch file from os.path.join(TEMPLATES_FOLDER, "blank_batch_" + self.lc_column + ".lcb")...
        batch = None  # ...save to "batch"

        if batch is None:
            print(f"batch name is: blank_batch_{self.lc_column}.lcb")
            return False, "F_Failed to load template batch: try restarting LS with different batch file open"

        with disposable(batch):
            if self.lc_column == 'semiprep':
                if 'contents' not in self.FRC_plate[0].keys():
                    print("Designated fraction collection plate missing 'contents' field")
                    return False, "F_Designated fraction collection plate missing 'contents' field"
                elif self.FRC_plate[0]['contents'] == 'Empty':
                    start_well = 0
                else:
                    start_well = len(self.FRC_plate[0]['contents'])
                inj_vol = 50
                # if this is the first time this batch is running then reset the fraction collector
                # Redacted: Load details into the batch file, including a startup method
                #   from f"Semiprep_column_startup_{kwargs['plate_num']}_{VIAL_NAMES[start_well]}.lcm"
                self.samples.put((os.path.join(self.folderName, "X_startup.lcd"), -1, "none", ""))
                # Redacted: Progress to next item in the batch
            elif self.lc_column == 'analytical':
                inj_vol = 5
            elif self.lc_column == 'uplc':
                if self.lc_mode == 'photochem':
                    inj_vol = 2.5
                else:
                    inj_vol = 1
                # if this is the first time this batch is running then reset the fraction collector
                # Redacted: Load details into the batch file, including a startup method from "UPLC_column_startup.lcm"
                self.samples.put((os.path.join(self.folderName, "X_startup.lcd"), -1, "none", ""))
                # Redacted: Progress to next item in the batch
            else:
                return False, f"F_Unknown LC mode: {self.lc_column}"

            # Run through wells and populate batch file lines based on the contents of each well            
            for well in batch_dict['method_file'].keys():
                # Redacted: Load details into the batch file
                self.samples.put((os.path.join(self.folderName, well + ".lcd"),
                                  batch_dict['mol_mw'][well],
                                  well.split('_')[0],
                                  batch_dict['target_SMILES'][well]))
                # Redacted: Progress to next item in the batch

            # Redacted: Load details into the batch file, including a column storage method
            #   from os.path.join(TEMPLATES_FOLDER, "Analytical_column_storage.lcm")
            # Redacted: Progress to next item in the batch
            self.samples.put((os.path.join(self.folderName, "X_cleanup.lcd"), -1, "none", ""))
            # Redacted: Load details into the batch file, including a column storage method
            self.samples.put((os.path.join(self.folderName, "X_shutdown.lcd"), -1, "none", ""))

            # Overwrite template batch file
            # Redacted: Save the batch file as batch_dict['batch_file_name']
            #   If the result was not a successful one,
            #   return (False, "F_Failed to save batch file as: " + batch_dict['batch_file_name'])

        # Load and execute batch
        # Redacted: Load batch from batch_dict['batch_file_name'], recording return code
        #   If the code is not a successful one, return (False, f"F_{additional details}")
        self.update_time_est()
        time.sleep(5)
        # Redacted: Start batch, recording return code
        #   If the code is not a successful one, return (False, f"F_{additional details}")

        return True, "Run " + self.lc_column + " Batch Command Received Successfully"

    def run_analytic_batch(self, queue_id, operation_id):
        # check that you're not running the same reaction multiple times
        # check that Shimadzu has been initialized and in the expected state
        result, update_msg = self.update_from_queue(queue_id, operation_id)
        if result:
            batch_file_name = update_msg
        else:
            return False, update_msg

        self.valve.write(b'GOB\r\n')
        self.valve.write(b'CP\r\n')
        if 'CPA' not in str(self.valve.readline()):
            print(f"Problem switching to Analytical Column, check Column Switching Valve at {CV_PORT}")
            return False, 'F_Failed to switch to Analytical Column'

        batch_dict = {'method_file': {},
                      'mol_mw': {},
                      'target_SMILES': {},
                      'batch_file_name': batch_file_name}
        unique_wells = []
        final_products = False


        for well, info in self.wellplate['contents'].items():
            # should reactions with slightly different concentrations be considered unique?
            #rxn_identifier = [info['reaction_smiles'], info['reagents']]
            if info['confirmation'] != 'none':
                continue
            if info['final_product'][1] == 'yes':
                final_products = True

            if 'hplc_processing_method' in self.operation['details'].keys():
                if self.operation['details']['hplc_processing_method'] == 'photochemistry':
                    self.lc_mode = 'photochem'
                    batch_dict['method_file'][well] = os.path.join(TEMPLATES_FOLDER, "photochem_UPLC.lcm")
                    batch_dict['mol_mw'][well] = -1
                    batch_dict['target_SMILES'][well] = info['reaction_smiles']
                    continue
                # kinetics does not care if the reactions all have the same product
                elif self.operation['details']['hplc_processing_method'] == 'kinetics':
                    self.lc_mode = 'kinetics'
                    batch_dict['method_file'][well] = os.path.join(TEMPLATES_FOLDER, "binary_grad_45min.lcm")
                    batch_dict['mol_mw'][well] = -1
                    # Need a unique identifier for each aliquot so that the HPLC recognizes them as different
                    batch_dict['target_SMILES'][well] = well
                    continue
                else:
                    return False, f"{self.operation['details']['hplc_processing_method']} is not valid HPLC method"
            else:
                mw = float(requests.post('http://localhost:40888',
                                         json={"exact_mass": info['target_product'][0][0]}).content.decode('utf-8'))

                rxn_identifier = info['reaction_smiles']
                if rxn_identifier in unique_wells:
                    continue
                unique_wells.append(rxn_identifier)

                print(f"{well}: {info['target_product'][0][0]} = {mw} g/mol")
                if mw < 350:
                    method_file = os.path.join(TEMPLATES_FOLDER, "Analytical_100_350MWcutoff.lcm")
                elif mw < 600:
                    method_file = os.path.join(TEMPLATES_FOLDER, "Analytical_350_600MWcutoff.lcm")
                elif mw < 850:
                    method_file = os.path.join(TEMPLATES_FOLDER, "Analytical_600_850MWcutoff.lcm")
                else:
                    return False, f"F_Not yet configured for products with a MW > 850. Target MW = {mw}"
                batch_dict['method_file'][well] = method_file
                batch_dict['mol_mw'][well] = mw
                batch_dict['target_SMILES'][well] = info['target_product']

        stat, msg = self.check_new_reagents()
        if stat:
            num_extra_ops = 7
        else:
            num_extra_ops = 0

        # Add the operations for the semiprep batch and then move that container to the liquid handler
        if final_products and (
                self.queue['operations'][str(int(operation_id) + 1)]['operation'] != 'run_semiprep_batch'):
            print("Adding a semiprep batch to the queue")
            new_container_k = self.add_plate_to_queue('fraction_plate')
            new_op = {'agent': 'LC',
                      'completed': 'no',
                      'container': self.operation['container'],
                      'end_time': None,
                      'start_time': None,
                      'operation': 'run_semiprep_batch',
                      'time_est': 3000,
                      'details': {'target_container': new_container_k,
                                  'queue_number': self.operation['details']['queue_number'],
                                  'is_paired': 'yes'}}
            dbi.insert_queue_steps(self.queue['queue_name'], int(operation_id) + num_extra_ops, [new_op])
            self.queue, _, _ = dbi.query_document('LC', 'queue', 'queue_name', queue_id)
        return self.run_batch(batch_dict)

    def run_semiprep_batch(self, queue_id, operation_id):
        # check that Shimadzu has been initialized and in the expected state
        result, update_msg = self.update_from_queue(queue_id, operation_id)
        if result:
            batch_file_name = update_msg
        else:
            return False, update_msg

        self.valve.write(b'GOB\r\n')
        self.valve.write(b'CP\r\n')
        if 'CPB' not in str(self.valve.readline()):
            print(f"Problem switching to Semiprep Column, check Column Switching Valve at {CV_PORT}")
            return False, 'F_Failed to switch to Semiprep Column'

        batch_dict = {'method_file': {},
                      'mol_mw': {},
                      'target_SMILES': {},
                      'batch_file_name': batch_file_name}
        unique_wells = []
        wettest_method = 1
        try:
            collect_wells = self.operation['details']['specified_wells']
        except KeyError:
            collect_wells = []
        if collect_wells:
            for well in collect_wells:
                confirmation = ''
                smiles = ''
                try:
                    confirmation = self.wellplate['contents'][well]['confirmation']
                    smiles = self.wellplate['contents'][well]['target_product'][0][0]
                except KeyError:
                    print(f'specified sample well ({well}) not sample plate ({self.wellplate["container_name"]})')
                    continue
                if confirmation == 'none':
                    return False, "Analytical batch not yet run for " + queue_id
                elif confirmation == 'no mass hit':
                    continue
                counter = 1
                while well in batch_dict['method_file'].keys():
                    well = well.split('_')[0] + f'_{counter}'
                    counter += 1

                ret_time = float(confirmation) * 0.9 + 0.55
                method_file_num = int(36 - 6 * ret_time)
                if method_file_num == 0:
                    method_file_num = 1
                wettest_method = max(wettest_method, method_file_num)
                lcm = "Isolate_" + str(method_file_num) + ".lcm"
                method_file = os.path.join(TEMPLATES_FOLDER, lcm)
                batch_dict['method_file'][well] = method_file
                mw = float(requests.post('http://localhost:40888',
                                         json={"exact_mass": smiles}).content.decode('utf-8'))
                print(f"{well}: {smiles} = {mw} g/mol")
                batch_dict['mol_mw'][well] = mw
                batch_dict['target_SMILES'][well] = self.wellplate['contents'][well]['target_product']
        else:
            for well, info in self.wellplate['contents'].items():
                if info['confirmation'] == 'none':
                    return False, "Analytical batch not yet run for " + queue_id
                elif info['confirmation'] == 'no mass hit':
                    continue
                elif (info['final_product'][1] == 'yes') or collect_all_unique_products:
                    if info['target_product'][0][0] in unique_wells:
                        continue
                    unique_wells.append(info['target_product'][0][0])
                    ret_time = float(info['confirmation'])
                    method_file_num = int(36 - 6 * ret_time)
                    if method_file_num == 0:
                        method_file_num = 1
                    wettest_method = max(wettest_method, method_file_num)
                    lcm = "Isolate_"+str(method_file_num)+".lcm"
                    method_file = os.path.join(TEMPLATES_FOLDER, lcm)
                    batch_dict['method_file'][well] = method_file
                    mw = float(requests.post('http://localhost:40888',
                                             json={"exact_mass": info['target_product'][0][0]}).content.decode('utf-8'))
                    print(f"{well}: {info['target_product'][0][0]} = {mw} g/mol")
                    batch_dict['mol_mw'][well] = mw
                    batch_dict['target_SMILES'][well] = info['target_product']
        num_products = len(batch_dict['method_file'])
        if num_products == 0:
            return True, "No product hits for " + queue_id

        # find any plates already in the fraction collector
        tot_num_wells = 95
        query_location, _, _ = dbi.query_location('LC', 'fraction_collector')
        frc_locations = {}
        for match in query_location['partial_matches'].values():
            if match['labware_type'] != '96 Well DeepWell':
                return False, "F_Disallowed plate type in FRC. Need: '96 Well DeepWell' not '{match['labware_type']}'"
            sublocation = str(match['location'][1][0])
            num_empty_wells = tot_num_wells if match['contents'] == 'Empty' else tot_num_wells - len(match['contents'])
            frc_locations[sublocation] = {'container_name': match['container_name'],
                                          'num_empty_wells': num_empty_wells}

        # Manual override to use a specific plate
        op_details = self.operation[QOP_DETAILS]
        if 'user_specified_plate' in op_details.keys() and op_details['user_specified_plate'] == 'yes':
            frc_container_name = self.queue['containers'][op_details['target_container']]['container_name']
            if frc_container_name:
                try:
                    frc_container, _, _ = dbi.query_document('LC', 'wellplates', 'container_name', frc_container_name)
                except cexc.DatabaseRequestError:
                    return (False,
                            f"F_FRC container defined ({frc_container_name}), but does not appear in database")
                if (tot_num_wells - len(frc_container['contents'])) < num_products:
                    return (False,
                            f"F_The specified fraction plate ({frc_container_name}) does not have enough empty wells")
                if frc_container['location'][0] != 'fraction_collector':
                    if frc_container['location'][0] == 'lpx':
                        new_ops = []
                        if len(frc_locations) == 3:
                            new_ops = self.remove_plate_from_frc(frc_locations['1']['container_name'])
                        new_ops = new_ops + self.fetch_plate_from_lpx(frc_container_name)
                        dbi.insert_queue_steps(queue_id, int(operation_id), new_ops)
                        return (False,
                                f"B_Need to fetch the specified fraction plate ({frc_container_name}) from the lpx")
                    else:
                        return (False,
                                f"F_The fraction plate ({frc_container_name}) is not accessible (in FRC or lpx)")
                else:
                    self.FRC_plate = [frc_container]
                    return self.run_batch(batch_dict, plate_num=frc_container['location'][1][0])

        # find any fraction plates associated with the campaign
        campaign_name = self.queue['queue_name'].rsplit('_', 2)[0]
        campaign_name = campaign_name
        campaign_dict, _, _ = dbi.query_document('LC', 'wellplates', 'campaign_name', campaign_name)
        plate_dict = {}
        pprint(f"Plates currently occupying the fraction collector:\n {frc_locations}")
        for campaign_plate_db_name, campaign_plate in campaign_dict.items():
            if campaign_plate.get('container_category', None) == 'fraction_plate':
                plate_dict[campaign_plate_db_name] = campaign_plate
        assigned_plate = []
        if len(plate_dict) == 0:
            # check if plates already on FRC are available for this campaign
            lowest_available_plate = 4
            for idx, info in frc_locations.items():
                if (info['num_empty_wells'] == tot_num_wells) and (int(idx) < lowest_available_plate):
                    lowest_available_plate = int(idx)
                    if assigned_plate:
                        assigned_plate[0] = info['container_name']
                    else:
                        assigned_plate.append(info['container_name'])
            # if no plate on FRC available, fetch a new plate from the LPX
            if not assigned_plate:
                aceso.resource_request(queue_name=queue_id,
                                       step=operation_id,
                                       target_destination='fraction_collector',
                                       is_new=False,
                                       name_basis=self.operation['details']['target_container'],
                                       scan_barcode='no')
                if len(frc_locations) == 3:
                    dbi.insert_queue_steps(queue_id, int(operation_id),
                                           self.remove_plate_from_frc(frc_locations['1']['container_name']))
                return False, "B_Fraction collector needs new plate"

        elif len(plate_dict) > 0:
            # logic tree for getting the right plates on the FRC with minimum number of new queue operations
            for plate_details in plate_dict.values():
                num_empty_wells = tot_num_wells - len(plate_details['contents'])
                location = plate_details['location']
                if num_empty_wells == 0:
                    continue
                elif num_empty_wells > num_products:
                    if location[0] == 'fraction_collector':
                        assigned_plate.append(plate_details['container_name'])
                        break
                    elif location[0] == 'lpx':
                        # if the fraction collector is full remove the fraction plate in position 1
                        new_ops = []
                        if len(frc_locations) == 3:
                            new_ops = self.remove_plate_from_frc(frc_locations['1']['container_name'])
                        new_ops = new_ops + self.fetch_plate_from_lpx(plate_details['container_name'])
                        dbi.insert_queue_steps(queue_id, int(operation_id), new_ops)
                        return False, f"B_Fraction plate ({plate_details['container_name']}) needs to be fetched"
                    else:
                        # plate is not in an accessible location, report problem, maybe it becomes available
                        return (False,
                                f"P_Fraction plate ({plate_details['container_name']}) location is {location[0]}, "
                                f"must be lpx or fraction_collector")
                else:  # this is where collection will spill over into a second DWP
                    if location[0] == 'fraction_collector':
                        new_ops = []
                        request_new_plate = False
                        if len(frc_locations) == 3:
                            if int(location[1][0]) == 1:
                                if frc_locations['2']['num_empty_wells'] != tot_num_wells:
                                    new_ops = self.remove_plate_from_frc(frc_locations['2']['container_name'])
                                    request_new_plate = True
                                else:
                                    assigned_plate = [frc_locations['1']['container_name'],
                                                      frc_locations['2']['container_name']]
                                    break
                            elif int(location[1][0]) == 2:
                                if frc_locations['3']['num_empty_wells'] != tot_num_wells:
                                    new_ops = self.remove_plate_from_frc(frc_locations['3']['container_name'])
                                    request_new_plate = True
                                else:
                                    assigned_plate = [frc_locations['2']['container_name'],
                                                      frc_locations['3']['container_name']]
                                    break
                            elif int(location[1][0]) == 3:
                                new_ops = self.remove_plate_from_frc(frc_locations['2']['container_name'])
                                new_ops = new_ops + self.remove_plate_from_frc(frc_locations['3']['container_name'],
                                                                               '', 'fraction_collector')
                                request_new_plate = True
                        elif len(frc_locations) == 2:
                            if int(location[1][0]) == 1:
                                if '2' not in frc_locations.keys():
                                    request_new_plate = True
                                elif frc_locations['2']['num_empty_wells'] != tot_num_wells:
                                    new_ops = self.remove_plate_from_frc(frc_locations['2']['container_name'])
                                    request_new_plate = True
                                else:
                                    assigned_plate = [frc_locations['1']['container_name'],
                                                      frc_locations['2']['container_name']]
                                    break
                            elif int(location[1][0]) == 2:
                                if '3' not in frc_locations.keys():
                                    request_new_plate = True
                                elif frc_locations['3']['num_empty_wells'] == tot_num_wells:
                                    assigned_plate = [frc_locations['2']['container_name'],
                                                      frc_locations['3']['container_name']]
                                    break
                                else:
                                    new_ops = self.remove_plate_from_frc(frc_locations['2']['container_name'],
                                                                         '', 'fraction_collector')
                                    request_new_plate = True
                            elif int(location[1][0]) == 3:
                                if '1' not in frc_locations.keys():
                                    if frc_locations['2']['num_empty_wells'] != tot_num_wells:
                                        new_ops = self.remove_plate_from_frc(frc_locations['2']['container_name'],
                                                                             '', 'fraction_collector')
                                        request_new_plate = True
                                else:
                                    request_new_plate = True
                                new_ops = new_ops + self.remove_plate_from_frc(frc_locations['3']['container_name'],
                                                                               '', 'fraction_collector')
                        elif len(frc_locations) == 1:
                            if int(location[1][0]) in [2, 3]:
                                new_ops = self.remove_plate_from_frc(frc_locations[location[1][0]]['container_name'],
                                                                     '', 'fraction_collector')
                            request_new_plate = True
                        else:
                            print("This should not print, already on FRC")
                        if request_new_plate:
                            aceso.resource_request(queue_name=queue_id,
                                                   step=operation_id,
                                                   target_destination='fraction_collector',
                                                   is_new=False,
                                                   name_basis=self.operation['details']['target_container'],
                                                   scan_barcode='no')
                        dbi.insert_queue_steps(queue_id, int(operation_id), new_ops)
                        return False, "B_Need to rearrange fraction collector and possibly fetch new plate"
                    elif location[0] == 'lpx':
                        new_ops = []
                        request_new_plate = False
                        if len(frc_locations) == 3:
                            if frc_locations['2']['num_empty_wells'] == tot_num_wells:
                                new_ops = self.remove_plate_from_frc(frc_locations['1']['container_name'])
                            elif frc_locations['3']['num_empty_wells'] == tot_num_wells:
                                new_ops = self.remove_plate_from_frc(frc_locations['2']['container_name'])
                            else:
                                new_ops = self.remove_plate_from_frc(frc_locations['1']['container_name'])
                                new_ops = new_ops + self.remove_plate_from_frc(frc_locations['2']['container_name'])
                                request_new_plate = True
                        elif len(frc_locations) == 2:
                            if '3' not in frc_locations.keys():
                                if frc_locations['2']['num_empty_wells'] == tot_num_wells:
                                    new_ops = self.remove_plate_from_frc(frc_locations['1']['container_name'])
                                else:
                                    new_ops = self.remove_plate_from_frc(frc_locations['2']['container_name'])
                                    request_new_plate = True
                            elif '2' not in frc_locations.keys():
                                if frc_locations['2']['num_empty_wells'] != tot_num_wells:
                                    new_ops = self.remove_plate_from_frc(frc_locations['1']['container_name'])
                                    request_new_plate = True
                            elif '1' not in frc_locations.keys():
                                if frc_locations['2']['num_empty_wells'] != tot_num_wells:
                                    new_ops = self.remove_plate_from_frc(frc_locations['2']['container_name'])
                                    request_new_plate = True
                        elif len(frc_locations) == 1:
                            if '2' in frc_locations.keys():
                                if frc_locations['2']['num_empty_wells'] != tot_num_wells:
                                    new_ops = new_ops + self.remove_plate_from_frc(frc_locations['2']['container_name'],
                                                                                   '', 'fraction_collector')
                                    request_new_plate = True
                            else:
                                request_new_plate = True
                        else:
                            print('this should not print, coming from lpx')
                        new_ops = new_ops + self.fetch_plate_from_lpx(plate_details['container_name'])
                        if request_new_plate:
                            aceso.resource_request(queue_name=queue_id,
                                                   step=operation_id,
                                                   target_destination='fraction_collector',
                                                   is_new=False,
                                                   name_basis=self.operation['details']['target_container'],
                                                   scan_barcode='no')
                        dbi.insert_queue_steps(queue_id, int(operation_id), new_ops)
                        return False, "B_Need to fetch plate(s) and possibly rearrange fraction collector"
            else:
                print("This should not print, exiting for loop without break or return")

        # Update the plate in the FRC to be the fraction plate
        end_ops = []
        self.FRC_plate = []
        plate_position = 3
        for frc_plate in assigned_plate:
            update_q_containers = True
            frc_container, _, _ = dbi.query_document('LC', 'wellplates', 'container_name', frc_plate)
            old_frc_container = copy.deepcopy(frc_container)
            frc_container['campaign_name'] = campaign_name
            frc_container['container_category'] = 'fraction_plate'
            dbi.update_document('LC', 'wellplates', old_frc_container, frc_container)
            for q_container in self.queue['containers'].values():
                if q_container['container_name'] == frc_container['container_name']:
                    # check that all the information lines up??
                    update_q_containers = False
                    break
            if update_q_containers:
                self.add_plate_to_queue('fraction_plate', frc_plate)
                print(f"Is this where extra fraction plates are coming from? {frc_container['container_name']}")
            self.FRC_plate.append(frc_container)
            plate_position = min(plate_position, frc_container['location'][1][0])
            # If last queue in campaign, send to heater shaker
            if 'queue_number' in op_details.keys():
                if int(op_details['queue_number'][0]) == int(op_details['queue_number'][1]):
                    end_ops = end_ops + self.remove_plate_from_frc(frc_plate, '', 'heater_shaker')

        # If a plate is going to fill up, send it to a heater shaker to evaporate with H12 dummy well to monitor evap
        H12_well_num = 0
        if len(self.FRC_plate) > 1:
            H12_well_num = 96 - len(self.FRC_plate[0]['contents'])
            # need to reorder the batch dict....
            temp_batch_dict = {'method_file': {},
                               'mol_mw': {},
                               'target_SMILES': {},
                               'batch_file_name': batch_file_name}
            for i, well in enumerate(batch_dict['method_file'].keys()):
                temp_batch_dict['method_file'][well] = batch_dict['method_file'][well]
                temp_batch_dict['mol_mw'][well] = batch_dict['mol_mw'][well]
                temp_batch_dict['target_SMILES'][well] = batch_dict['target_SMILES'][well]
                if i == H12_well_num:
                    temp_batch_dict['method_file']['-1'] = os.path.join(TEMPLATES_FOLDER,
                                                                        "Isolate_" + str(wettest_method) + ".lcm")
                    temp_batch_dict['mol_mw']['-1'] = -1
                    temp_batch_dict['target_SMILES']['-1'] = 'reference_well'
            batch_dict = temp_batch_dict
            end_ops = end_ops + self.remove_plate_from_frc(assigned_plate[0], '', 'heater_shaker')
            if int(self.FRC_plate[0]['location'][1][0]) != (int(self.FRC_plate[1]['location'][1][0]) + 1):
                print("This should not print, FRC plates must be in consecutive FRC positions!")
                return False, "F_fraction plate locations not consecutive! This should not print"

        dbi.insert_queue_steps(queue_id, int(operation_id) + 2, end_ops)
        self.queue, _, _ = dbi.query_document('LC', 'queue', 'queue_name', queue_id)

        return self.run_batch(batch_dict, plate_num=plate_position)

    def run_optimization_batch(self, queue_id, operation_id):
        self.run_analytic_batch(queue_id, operation_id)

    def resume_batch(self, queue_id, operation_id):
        status = self.status_ping()
        if status == "Busy":
            return False, "B_Cannot resume batch while instrument is running"
        elif status == "Problem":
            return False, "P_Instrument in problem state"
        elif status == "Fatal":
            return False, "F_Critical instrument error, cannot resume batch"
        
        # If resuming after stop with continuous LC_module operation (i.e. not a total crash)
        #  rebuild batch_dict from self.samples
        if not self.samples.empty():
            batch_dict = {'method_file': {},
                          'mol_mw': {},
                          'target_SMILES': {},
                          'batch_file_name': os.path.join(self.folderName, "resume.lcb")}
            while not self.samples.empty():
                sample = self.samples.get()
                well = sample[2]
                mw = float(sample[1])
                if mw < 350:
                    method_file = os.path.join(TEMPLATES_FOLDER, "Analytical_100_350MWcutoff.lcm")
                elif mw < 600:
                    method_file = os.path.join(TEMPLATES_FOLDER, "Analytical_350_600MWcutoff.lcm")
                elif mw < 850:
                    method_file = os.path.join(TEMPLATES_FOLDER, "Analytical_600_850MWcutoff.lcm")
                else:
                    return False, f"F_Not yet configured for products with a MW > 850. Target MW = {mw}"
                batch_dict['method_file'][well] = method_file
                batch_dict['mol_mw'][well] = mw
                batch_dict['target_SMILES'][well] = sample[3]

        # If resuming after a module restart: need to rebuild sample queue from database
        else:
            result, update_msg = self.update_from_queue(queue_id, operation_id)
            if result:
                batch_file_name = update_msg.split('.')[0] + '_resume.lcb'
            else:
                return False, update_msg

            resume_folder = os.path.join(self.folderName, "resume_" + datetime.now().strftime('%d-%m-%Y_%H-%M-%S'))
            self.folderName = resume_folder
            batch_dict = {'method_file': {},
                          'mol_mw': {},
                          'target_SMILES': {},
                          'batch_file_name': batch_file_name}
            
            for well, info in self.wellplate['contents'].items():
                if info['confirmation'] != 'none':
                    continue
                mw = float(requests.post('http://localhost:40888',
                                         json={"exact_mass": info['target_product'][0][0]}).content.decode('utf-8'))
                print(f"{well}: {info['target_product'][0][0]} = {mw} g/mol")
                batch_dict['mol_mw'][well] = mw
                if mw < 350:
                    method_file = os.path.join(TEMPLATES_FOLDER, "Analytical_100_350MWcutoff.lcm")
                elif mw < 600:
                    method_file = os.path.join(TEMPLATES_FOLDER, "Analytical_350_600MWcutoff.lcm")
                elif mw < 850:
                    method_file = os.path.join(TEMPLATES_FOLDER, "Analytical_600_850MWcutoff.lcm")
                else:
                    return False, "P_Cannot find method file for mw > 850"
                batch_dict['method_file'][well] = method_file
                batch_dict['mol_mw'][well] = mw
                batch_dict['target_SMILES'][well] = info['target_product']
        
        self.running = (False, "Resuming")
        return self.run_batch(batch_dict)

    def add_plate_to_queue(self, name_prefix, database_container_name=None,
                           plate_type='96 Well DeepWell', category='fraction_plate'):
        for q_container_name, container in self.queue['containers'].items():
            if database_container_name == container['container_name']:
                return q_container_name
        new_container_k = name_prefix
        count = 1
        while new_container_k in self.queue['containers']:
            if self.queue['containers'][new_container_k]['contents'] == 'Empty':
                break
            new_container_k = name_prefix + '_' + str(count)
            count += 1

        if database_container_name:
            try:
                new_container_complete_v, _, _ = dbi.query_document('LC', 'wellplates', 'container_name',
                                                                    database_container_name)
                new_container_v = {'container_name': new_container_complete_v['container_name'],
                                   'plate_type': new_container_complete_v['labware_type'],
                                   'contents': new_container_complete_v['contents'],
                                   'container_category': category}
            except cexc.DatabaseRequestError:
                print(f"Database query for {database_container_name} failed... likely not in DB")

        else:
            new_container_v = {'container_name': None,
                               'plate_type': plate_type,
                               'contents': 'Empty',
                               'container_category': category}
        old_containers_queue, _, _ = dbi.query_document('LC', 'queue', 'queue_name', self.queue['queue_name'])
        new_containers_queue = copy.deepcopy(old_containers_queue)
        containers = new_containers_queue['containers']
        containers[new_container_k] = new_container_v
        new_containers_queue['containers'] = containers
        dbi.update_document('LC', 'queue', old_containers_queue, new_containers_queue)
        self.queue = new_containers_queue
        return new_container_k

    def fetch_plate_from_lpx(self, db_container_name, q_container_name='', destination='fraction_collector'):
        new_ops = []
        if not db_container_name:
            db_container_name = self.queue['containers'][q_container_name]['container_name']
        elif not q_container_name:
            q_container_name = self.add_plate_to_queue('fraction_plate', db_container_name)
        else:
            print("Need to specify exactly one of either db_container_name or q_container_name")
            return new_ops
        try:
            check_plate, _, _ = dbi.query_document('LC', 'wellplates', 'container_name', db_container_name)
        except cexc.DatabaseRequestError:
            print(f"Unable to find {db_container_name} in wellplates collection")
            return new_ops
        if check_plate['location'][0] != 'lpx':
            print(f"Cannot fetch {db_container_name} from lpx because its location is {check_plate['location'][0]}")
            return new_ops
        new_op1 = {'agent': 'Ss',
                   'completed': 'no',
                   'container': q_container_name,
                   'end_time': None,
                   'start_time': None,
                   'operation': 'Ss.fetch',
                   'time_est': 30,
                   'details': {'is_paired': 'yes',
                               'schedule_time': 0}}
        new_ops.append(new_op1)
        new_op2 = {'agent': 'Ra',
                   'completed': 'no',
                   'container': q_container_name,
                   'end_time': None,
                   'start_time': None,
                   'operation': 'move_wellplate',
                   'time_est': 30,
                   'details': {'target_destination': 'fraction_collector',
                               'is_paired': 'yes',
                               'schedule_time': 0}}
        new_ops.append(new_op2)
        return new_ops

    def remove_plate_from_frc(self, db_container_name, q_container_name='', destination='lpx'):
        # can ignore db_container_name if q_container_name is specified
        if destination not in ['lpx', 'storage_hotel', 'heater_shaker', 'fraction_collector']:
            print(f"Can not transfer plate from FRC to {destination}")
            return []
        if not q_container_name:
            q_container_name = self.add_plate_to_queue('fraction_plate', db_container_name)
        # Create operations for move and stow
        if destination in ['storage_hotel', 'heater_shaker']:
            ra_dest = 'liquid_handler'
            is_paired = 'yes'
        else:
            ra_dest = destination
            is_paired = 'yes'
        if destination == 'lpx':
            new_op1 = {'agent': 'Ra',
                       'completed': 'no',
                       'container': q_container_name,
                       'end_time': None,
                       'start_time': None,
                       'operation': 'move_wellplate',
                       'time_est': 30,
                       'details': {'target_destination': ra_dest,
                                   'is_paired': is_paired,
                                   'schedule_time': 0}}
            new_op2 = {'agent': 'Ss',
                       'completed': 'no',
                       'container': q_container_name,
                       'end_time': None,
                       'start_time': None,
                       'operation': 'Ss.stow',
                       'time_est': 30,
                       'details': {'is_paired': is_paired,
                                   'schedule_time': 0}}
            return [new_op1, new_op2]
        elif destination == 'storage_hotel':
            new_op0 = {'agent': 'Lh',
                       'completed': 'no',
                       'container': q_container_name,
                       'end_time': None,
                       'start_time': None,
                       'operation': 'transfer_wellplate',
                       'time_est': 30,
                       'details': {'target_destination': 'transfer_prep',
                                   'is_paired': is_paired,
                                   'schedule_time': 0}}
            new_op1 = {'agent': 'Ra',
                       'completed': 'no',
                       'container': q_container_name,
                       'end_time': None,
                       'start_time': None,
                       'operation': 'move_wellplate',
                       'time_est': 30,
                       'details': {'target_destination': ra_dest,
                                   'is_paired': is_paired,
                                   'schedule_time': 0}}
            new_op2 = {'agent': 'Lh',
                       'completed': 'no',
                       'container': q_container_name,
                       'end_time': None,
                       'start_time': None,
                       'operation': 'transfer_wellplate',
                       'time_est': 30,
                       'details': {'target_destination': 'transfer_cleanup',
                                   'is_paired': is_paired,
                                   'schedule_time': 0}}
            new_op3 = {'agent': 'Lh',
                       'completed': 'no',
                       'container': q_container_name,
                       'end_time': None,
                       'start_time': None,
                       'operation': 'transfer_wellplate',
                       'time_est': 30,
                       'details': {'target_destination': destination}}
            return [new_op0, new_op1, new_op2, new_op3]
        elif destination == 'heater_shaker':
            new_op0 = {'agent': 'Lh',
                       'completed': 'no',
                       'container': q_container_name,
                       'end_time': None,
                       'start_time': None,
                       'operation': 'transfer_wellplate',
                       'time_est': 30,
                       'details': {'target_destination': 'transfer_prep',
                                   'is_paired': is_paired,
                                   'schedule_time': 0}}
            new_op1 = {'agent': 'Ra',
                       'completed': 'no',
                       'container': q_container_name,
                       'end_time': None,
                       'start_time': None,
                       'operation': 'move_wellplate',
                       'time_est': 30,
                       'details': {'target_destination': ra_dest,
                                   'is_paired': is_paired,
                                   'schedule_time': 0}}
            new_op2 = {'agent': 'Lh',
                       'completed': 'no',
                       'container': q_container_name,
                       'end_time': None,
                       'start_time': None,
                       'operation': 'transfer_wellplate',
                       'time_est': 30,
                       'details': {'target_destination': 'transfer_cleanup',
                                   'is_paired': is_paired,
                                   'schedule_time': 0}}
            new_op3 = {'agent': 'Lh',
                       'completed': 'no',
                       'container': q_container_name,
                       'end_time': None,
                       'start_time': None,
                       'operation': 'transfer_wellplate',
                       'time_est': 30,
                       'details': {'target_destination': destination,
                                   'temperature': 79}}
            new_op4 = {'agent': 'Lh',
                       'completed': 'no',
                       'container': q_container_name,
                       'end_time': None,
                       'start_time': None,
                       'operation': 'start_stop_heater_shaker',
                       'time_est': 30,
                       'details': {'power': 'on',
                                   'rpms': 600,
                                   'temperature': 79}}
            new_op5 = {'agent': 'Lh',
                       'completed': 'no',
                       'container': q_container_name,
                       'end_time': None,
                       'start_time': None,
                       'operation': 'detect_liquid_level',
                       'time_est': 75,
                       'details': {'detection_well': 'H12'}}
            new_op6 = {'agent': 'Lh',
                       'completed': 'no',
                       'container': q_container_name,
                       'end_time': None,
                       'start_time': None,
                       'operation': 'start_stop_heater_shaker',
                       'time_est': 30,
                       'details': {'power': 'off',
                                   'rpms': 200,
                                   'temperature': 25}}
            new_op7 = {'agent': 'Lh',
                       'completed': 'no',
                       'container': q_container_name,
                       'end_time': None,
                       'start_time': None,
                       'operation': 'transfer_wellplate',
                       'time_est': 30,
                       'details': {'target_destination': 'storage_hotel'}}
            return [new_op0, new_op1, new_op2, new_op3, new_op4, new_op5, new_op6, new_op7]
        else:
            new_op1 = {'agent': 'Ra',
                       'completed': 'no',
                       'container': q_container_name,
                       'end_time': None,
                       'start_time': None,
                       'operation': 'move_wellplate',
                       'time_est': 30,
                       'details': {'target_destination': ra_dest,
                                   'is_paired': is_paired,
                                   'schedule_time': 0}}
            return [new_op1]
