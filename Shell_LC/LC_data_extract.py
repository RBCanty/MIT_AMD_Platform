#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 11:53:26 2021

@author: mattmcdonald
"""

import pandas as pd
import numpy as np
import os
import json
import shutil
import glob
import requests
import os.path
import time
import scipy.signal as sp
from pyteomics import mzxml
from datetime import datetime
from copy import deepcopy
from constants import WELL_IDS as vial_names
from ui_quick_dialog import tk, QuickUI


compound_data_root = r"C:\Shimadzu\Data\Automated\Compound_data"
compound_dictionary = r"C:\Shimadzu\Data\Automated\Compound_data\Compound_dictionary.json"

log_e_to_calc_file = 'log_e_to_calc.csv'
log_e_preds_file = 'log_e_predictions.csv'

def chemprop_server_ping():
    try:
        resp = requests.post('http://localhost:40888', json={"ping": "ping"})
        if "success" in resp.content.decode('utf-8'):
            return True
        else:
            return False
    except:
        return False

def log_e_chemprop_push(wellplate: dict, path: str):
    #return True
    wellplate['path'] = path
    try:
        resp = requests.post('http://localhost:40888', json=wellplate)
        #dict.pop(wellplate['path'])
        return resp.content.decode('utf-8')
    except:
        print("Error during chemprop prediction")

def log_e_chemprop_pull(wellplate: dict, path: str):
    """
    A function to pull log_e values that chemprop has already calculated, push them to data DB

    Parameters
    ----------
    wellplate : dict
        a standard wellplate dictionary, should be the one that log_e values will be pushed to.
    path : str
        the path to the chemprop data (use Shimadzu_API self.FolderName).

    Returns
    -------
    updated wellplate to be pushed into the plate database and used to update data DB

    """
    if "labware_type" in wellplate.keys():
        wellplate['plate_type'] = wellplate['labware_type']
        wellplate.pop('labware_type')
    log_e_vals = {}
    predictions_file = os.path.join(path, log_e_preds_file)
    if not os.path.exists(predictions_file):
        print("Chemprop predictions not available!")
        return wellplate
    with open(predictions_file, 'r') as preds_in:
        for line in preds_in.readlines():
            if "Results pending..." in line:
                print("Chemprop results not ready, please wait... sleeping 30s")
                time.sleep(30)
                return log_e_chemprop_pull(wellplate, path)
            elif "SMILES" in line:
                continue
            info = line.split(',')
            smiles = info[0].strip()
            solvent = info[1].strip()
            log_e = info[2].strip()
            log_e_unc = info[4].strip()
            if smiles in log_e_vals.keys():
                log_e_vals[smiles][solvent] = [float(log_e), float(log_e_unc)]
            else:
                log_e_vals[smiles] = {solvent: [float(log_e), float(log_e_unc)]}
        
    for key, well in wellplate['contents'].items():
        target = well['target_product'][0][0]
        if 'properties' in well.keys():
            wellplate['contents'][key]['properties'][target] = {'log_e': log_e_vals[target]}
        else:
            wellplate['contents'][key]['properties'] = {target: {'log_e': log_e_vals[target]}}
    # remove or archive result files??
    return wellplate

def calc_logP(ret_time: float):
    """
    Quick function to calculate logP from HPLC retention time

    Parameters
    ----------
    ret_time : float
        DESCRIPTION.

    Returns
    -------
    None.

    """
    # Fitting parameters
    slope = 1.667
    slope_unc = 0.1036
    intercept = -0.6912
    intercept_unc = 0.2916
    mean_fit_ret_time = 2.235
    sum_x2 = 88.14

    # Calculate LogP and uncertainty
    logP = slope*ret_time + intercept
    base_uncertainty = (slope_unc*ret_time)**2 + (intercept_unc**2)
    correlation_uncertainty = (sum_x2 - 2*ret_time*mean_fit_ret_time)*(slope_unc**2)
    logP_unc = (base_uncertainty + correlation_uncertainty)**0.5

    return [logP, logP_unc]

def add_reagent_prep_to_queue(up_to_date_queue: dict, wellplate: dict):
    """
    Returns two dictionaries, one for the queue steps and one for the new plate. Feed into aceso?

    Parameters
    ----------
    up_to_date_queue : dict
        A copy of the up to date queue, used to name the reagent plate
    wellplate : dict
        standard wellplate dictionary.

    Returns
    -------
    None.

    """
    new_wellplate = {'container_name': None, 'contents': {}, 'plate_type': '96 Well Microplate'}
    
    with open(compound_dictionary, 'r') as comps_in:
        standard_reagents = json.load(comps_in)
    
    new_reagents = []
    waiting_reactions = []
    new_well = 0
    for key, well in wellplate['contents'].items():
        for reagent in well['reagents']:
            smiles = reagent[2]
            if ((smiles not in standard_reagents.keys()) or 
                (well['solvents'][0][0] not in standard_reagents[smiles]['solvent'])):
                waiting_reactions.append(smiles + '__' + key)
                if smiles not in new_reagents:                
                    new_wellplate['contents'][vial_names[new_well]] = {'confirmation': 'none',
                                                                       'plate_well': vial_names[new_well],
                                                                       'reagents': [reagent],
                                                                       'solvents': well['solvents'],
                                                                       'target_product': [smiles, reagent[1]],
                                                                       'total_volume': well['total_volume']}
                    new_reagents.append(reagent[2])
                    new_well +=1
        
    if not new_reagents:
        return {}, {}
    
    plate_num = 1
    old_container_name = ''
    for container_name in up_to_date_queue['containers'].keys():
        if container_name['container_name'] == wellplate['container_name']:
            old_container_name = container_name
        if 'reagent_plate' in container_name:
            plate_num = max(int(container_name.split('_')[2]) + 1, plate_num)        
    new_container_name = f'reagent_plate_{plate_num}'
    
    new_ops = []
    # prepare reagent plate
    new_ops.append({'agent': 'Lh', 'completed': 'no', 'container': new_container_name,
                    'details': {'empty_wells': True, 'tip_prep_rinse': True, 'washing_frequency': None},
                    'end_time': None, 'operation': 'prepare_wellplate', 'start_time': None, 'time_est': 600})
    # move plate currently in autosampler to liquid_handler hotel
    new_ops.append({'agent': 'Ra', 'completed': 'no', 'container': old_container_name,
                    'details': {'target_destination': 'liquid_handler'}, 'end_time': None,
                    'operation': 'move_wellplate', 'start_time': None, 'time_est': 20})
    # move the reagent plate to the transfer hotel
    new_ops.append({'agent': 'Lh', 'completed': 'no', 'container': new_container_name,
                    'details': {'target_destination': 'transfer_hotel'}, 'end_time': None,
                    'operation': 'transfer_wellplate', 'start_time': None, 'time_est': 30})
    # move the reagent plate to the autosampler
    new_ops.append({'agent': 'Ra', 'completed': 'no', 'container': new_container_name,
                    'details': {'target_destination': 'autosampler'}, 'end_time': None, 
                    'operation': 'move_wellplate', 'start_time': None, 'time_est': 20})
    # run the analytical batch on the reagent plate
    new_ops.append({'agent': 'LC', 'completed': 'no', 'container': new_container_name,
                    'details': {'waiting_for_reagents': waiting_reactions}, 
                    'end_time': None, 'operation': 'run_analytical_batch',
                    'start_time': None, 'time_est': 360*(1+len(new_wellplate['contents']))})
    # move the reagent plate to (???, long term storage?, trash?) liquid handler
    new_ops.append({'agent': 'Ra', 'completed': 'no', 'container': new_container_name,
                    'details': {'target_destination': 'liquid_handler'}, 'end_time': None, 
                    'operation': 'move_wellplate', 'start_time': None, 'time_est': 20})
    # move the reaction plate back to the autosampler to begin semiprep batch
    new_ops.append({'agent': 'Ra', 'completed': 'no', 'container': old_container_name,
                    'details': {'target_destination': 'autosampler'}, 'end_time': None, 
                    'operation': 'move_wellplate', 'start_time': None, 'time_est': 20})
    
    new_container_queue = deepcopy(up_to_date_queue)
    new_container_queue['containers'][new_container_name] = new_wellplate
    return new_container_queue, new_ops

def lc_data_analysis(path: str, wellplate: dict, target_sample: tuple):
    """

    Parameters
    ----------
    path : str
        data save folder path
    wellplate : dict
        as defined in a standard queue document
    target_sample : tuple
        (.lcd file name, molecular weight, well, smiles)
        
    Returns
    -------
    reaction_info : dict
        information about the reaction.
    """
    ret_time, best_mz = ms_data_analysis(target_sample)
    yield_info = pda_data_analysis(path, wellplate, target_sample, ret_time)
    
    reaction_info = {'conc': yield_info['conc'],
                     'retention_time': ret_time,
                     'best_m/z': best_mz,
                     'PDA_Data': yield_info['chromatogram_report'],
                     'MS_Data': target_sample[0][:-4] + ".mzXML",
                     'reaction_class': 'placeholder, query ASKCOS for reaction classification',
                     'templates': wellplate['contents'][target_sample[2]]['templates']}
    # there will need to be more in reaction_info eventually, not sure what that is yet
    return reaction_info

def ms_data_analysis(target_sample):
    # peak identification parameters: number of consecutive points m std devs above the baseline noise
    n_con_pts = 10
    m_sigma = 4
    def baseline_helper(data, n_segs, n_toss):
        # helper function to find standard deviation of the baseline
        # randomly remove points to resize data into n_segs equal sized segments
        data = np.delete(data, np.random.permutation(len(data))[0:(len(data) % n_segs)])
        data = data.reshape((-1, n_segs), order='F')
        mus = np.mean(data, axis=0)
        idx = np.argsort(mus)
        data = np.delete(data, idx[-n_toss:], 1)
        return np.mean(data), np.std(data)

    mw = target_sample[1]
    # look for adducts [+proton, -proton, +ACN+H, -formate]... can add more (+Na, +K, etc.) but increases false pos.
    mws = np.array([mw + 1, 1 - mw, mw + 42])

    file_name = target_sample[0][:-4] + ".mzXML"
    experiment = mzxml.read(file_name)
    spec = experiment.next()
    # rows ought to be the total length of the experiment divided by 2... which is 360 but this is janky
    rows = int(360)
    cols = int(2500)
    upsample_pts = np.linspace(spec['startMz'], spec['endMz'], cols)

    # even scans are positive polarity, odd scans are negative polarity
    flipper = 0
    rt_pos = np.ndarray(shape=(rows, 1), dtype=np.float32)
    rt_neg = np.ndarray(shape=(rows, 1), dtype=np.float32)
    data_pos = np.ndarray(shape=(cols, rows), dtype=np.float32)
    data_neg = np.ndarray(shape=(cols, rows), dtype=np.float32)
    i = 0
    rt_pos[i] = spec['retentionTime']
    for spec in experiment:
        mzs = spec['m/z array']
        ints = spec['intensity array']
        rt = spec['retentionTime']
        if flipper:
            data_pos[:, i] = np.interp(upsample_pts, mzs, ints)
            flipper -= 1
            rt_pos[i] = rt
        else:
            data_neg[:, i] = np.interp(upsample_pts, mzs, ints)
            flipper += 1
            rt_neg[i] = rt
            i += 1

    mzs = upsample_pts
    ret_time = -1.0
    max_intensity = 0
    intensities = [] # list of [ret_time, summed_intensity] pairings
    best_mz = 0
    # indentify peaks based on N consecutive points > M*sigma baseline
    for mz in mws:
        if mz > 0:
            #print(mz)
            mic_low = np.where(mzs > (mz - 0.5))[0][0]
            mic_high = np.where(mzs < (mz + 0.5))[0][-1]
            mic_pos = np.sum(data_pos[mic_low:mic_high, :], 0)
            mic_pos = mic_pos / max(mic_pos)
            #plt.plot(rt_pos, mic_pos)
            mu, sigma = baseline_helper(mic_pos, 10, 2)
            #print(mu, sigma)
            idx = np.array(np.where(mic_pos > (mu + m_sigma * sigma)))
            for i in range(0, n_con_pts - 1):
                idx = idx[np.diff(idx, append=-1) == 1]
            j = 0 
            for i in idx:
                if i == j+1:
                    intensities[-1][0] += mic_pos[i]
                else:
                    intensities.append([0, float(rt_pos[i])])
                j = i
            
            if intensities:
                ret_time = max(intensities)[1]
                if max(intensities)[0] > max_intensity:
                    max_intensity = max(mic_pos)
                    best_mz = mz
        else:
            mic_low = np.where(mzs > (abs(mz) - 0.5))[0][0]
            mic_high = np.where(mzs < (abs(mz) + 0.5))[0][-1]
            mic_neg = np.sum(data_neg[mic_low:mic_high, :], 0)
            mic_neg = mic_neg / max(mic_neg)
            mu, sigma = baseline_helper(mic_neg, 10, 2)
            idx = np.array(np.where(mic_neg > (mu + m_sigma * sigma)))
            for i in range(0, n_con_pts - 1):
                idx = idx[np.diff(idx, append=-1) == 1]
            if idx.size:
                if max(mic_neg) > max_intensity:
                    ret_time = rt_neg[idx[0]]
                    max_intensity = max(mic_neg)
                    best_mz = mz
            j = 0 
            for i in idx:
                if i == j+1:
                    intensities[-1][0] += mic_pos[i]
                else:
                    intensities.append([0, float(rt_pos[i])])
                j = i
            
            if intensities:
                ret_time = max(intensities)[1]
                if max(intensities)[0] > max_intensity:
                    max_intensity = max(mic_pos)
                    best_mz = mz

    if (0.6 > ret_time > 0) or (ret_time > 5.7):
        print(f"Peak detected too early/late to isolate ({ret_time})... possibly unrealistic ")
        ret_time = -1.0

    return ret_time, best_mz

def pda_data_analysis(path: str, wellplate: dict, target_sample: tuple, ret_time: float):
    """

    Parameters
    ----------
    path : str
        data save folder path
    wellplate : dict
        as defined in a standard queue document
    target_sample : tuple
        (.lcd file name, molecular weight, well, [[smiles, amount]])
    ret_time : float
        retention time as calculated by MS

    Returns
    -------
    conc_info : dict
        information about the reaction.

    """
    
    if ret_time < 0:
        yield_est = 0.0
        conc_info = {'conc': yield_est,
                     'chromatogram_csv': 'no mass hit',
                     'chromatogram_report': 'temp report file'}
    else:
        smiles = target_sample[3][0][0]
        compound_file = save_txt_data(path, smiles, target_sample[2])
        preds_file = os.path.join(path, log_e_preds_file)
        # do analysis with the compound_file using MOCCA
        
        conc = calc_amount(compound_file, preds_file, smiles, ret_time)
        
        
        conc_info = {'conc': conc,
                     'chromatogram_csv': compound_file,
                     'chromatogram_report': 'temp report file'}
    return conc_info

def save_txt_data(data_path, smiles, well, save_path=r"C:\Shimadzu\Data\Automated\Compound_data"):
    """
    Reads the 3D data exported by the Shimadzu software. Saves data to local folder
    Parameters
    ----------
    path : str
        The directory, in which the experimental data are stored.

    Returns
    -------
    save_file : str
        path to txt of data
    """

    if data_path[-4:] == '.lcd':
        data_path = os.path.split(data_path)[0]
    shimadzu_file = os.path.join(data_path, 'PDA_Data.txt')

    # build a local database of compound PDA data and save new compounds by smile string and date recorded
    save_file_root = os.path.join(save_path, datetime.now().strftime('%m%d%Y'))
    if not os.path.exists(save_file_root):
        os.mkdir(save_file_root)
    filename = smiles_to_filename(smiles)
    save_file = os.path.join(save_file_root, filename + "__" + well + "_0.txt")
    filenum = 0
    while os.path.exists(save_file):
        filenum += 1
        save_file = os.path.join(save_file_root, filename + "__" + well + "_" + str(filenum) + ".txt")
    # eventually need a way to record which sample is the one to be used for the next reaction
    shutil.move(shimadzu_file, save_file)
    return save_file

def smiles_to_filename(in_string, reverse=False):
    """
    Takes a smiles and generates a filename corresponding to that smiles without disallowed characters

    Parameters
    ----------
    in_string : str
        smiles string for compound, should be in rdkit standard form.
    reverse : bool
        set to True to go from filename to smiles

    Returns
    -------
    filename : str
        filename to be used in local storage as .csv

    """
    disallowed_chars = ['.', '%', '\\', '/', ':']
    replacement_chars = ['(dot)', '(pct)', '(bs)', '(fs)', '(cln)']
    
    out_string = in_string
    
    if reverse:
        for i, char in enumerate(replacement_chars):
            out_string = out_string.replace(char, disallowed_chars[i])
    else:
        for i, char in enumerate(disallowed_chars):
            out_string = out_string.replace(char, replacement_chars[i])
    
    return out_string

def smiles_file_lookup(smiles: str, well=""):
    """
    Lookup a SMILES string to see if there is data associated with that SMILES (for reagents)

    Parameters
    ----------
    smiles : str
        SMILES string for the reagent of interest.
    well : str
        Specify a given well, for referencing previous iterations of optimization. If not given will 
        return the most recent version of that chemical. If specified SMILES-well comnbination does not 
        exist print warning

    Returns
    -------
    hit_file : str
        The name of the text file for the compound (empty string if file doesn't exist')
    """
    file_smiles = smiles_to_filename(smiles)
    all_compound_files = glob.iglob(compound_data_root + r'\*\*.txt')
    hit_file = ""
    if well:
        file_smiles = file_smiles + '__' + well
        hit_file = "x"
    most_recent = 0
    for file in all_compound_files:
        if file_smiles in file:
            filenum = file.split('_')[-1]
            filenum = filenum.split('.')[0]
            if int(filenum) > most_recent:
                hit_file = file
                most_recent = int(filenum)
    if hit_file == "x":
        print(f'Warning: No data HPLC for {smiles} in well {well}!')
        hit_file = ""
    return hit_file
    
def run_optimization(queue, contents):
    """
    

    Parameters
    ----------
    contents : TYPE
        DESCRIPTION.

    Returns
    -------
    new_containers_queue : dict
        The same queue but with the containers updated to reflect the new containers that 
        are needed to do the optimization
    new_queue_ops : dict
        New opererations that need to be inserted into the queue

    """
    # need to clone the containers section of the queue
    old_containers = queue['containers']
    new_containers = {}

    new_ops = []
    all_ops = queue['operations']
    num_optimization_iters = 0
    for i in range(len(all_ops)):
        op = all_ops[str(i+1)]
        if op['operation'] == 'complete_queue':
            break 
        
        new_container = op['container'] + "_" + num_optimization_iters
        new_container_dict = old_containers[op['container']]
        
        # remember to skip over reagent plate preps
        if (op['operation'] == 'prepare_wellplate') and ('reagent' not in op['container']):
            new_ops = []
            new_containers = {}
            num_optimization_iters += 1
        else:
            # TODO! create new_contents generator
            new_container_dict['contents'] = {}

        op_clone = deepcopy(op)
        op_clone['completed'] = 'no'
        op_clone['container'] = new_container
        new_ops.append(op_clone)
        
        if new_container not in new_containers.keys():
            new_containers[new_container] = new_container_dict
    
    for k, container in old_containers.items():
        if k not in new_containers.keys():
            new_containers[k] = container
    
    """
    need a different function
    """
    
    # need to edit new_ops[0] to do next iteration of optimization
    new_containers_queue = deepcopy(queue)
    new_containers_queue['containers'] = new_containers
        
    if num_optimization_iters > 3:  # arbitrarily limit to 3 iters, need to rethink constraint
        new_ops = []
    return new_containers_queue, new_ops

def initialize_optimization_plate(reactions, original_plate):
    """
    Create initial plates for optimization (i.e. no prior knowledge)

    Parameters
    ----------
    reactions : TYPE
        DESCRIPTION.
    original_plate : TYPE
        DESCRIPTION.

    Returns
    -------
    A list of new plates to run... in the event that all the optimization doesn't fit on one plate
    """

    num_rxns = 0
    num_wells = 0
    next_plate = {}
    new_plates = []
    for well, contents in original_plate['contents'].values():
        for reaction in reactions:
            if reaction in contents['reaction_smiles']:
                new_contexts = generate_new_contexts(contents)
                num_wells += len(new_contexts)
                if num_wells > 96:
                    new_plates.append(next_plate)
                    next_plate = {}
                    num_wells = len(new_contexts)
                    num_rxns = 0
                for new_context in new_contexts:
                    next_plate[vial_names[num_rxns]] = new_context
                    num_rxns += 1
    new_plates.append(next_plate)
    return new_plates

def generate_new_contexts(contents):
    template_contents = deepcopy(contents)
    template_contents['confirmation'] = 'none'
    new_contexts = []
    for i, reagent in enumerate(template_contents['reagents']):
        template_contents[i] = [reagent[0], reagent[1]/2]
        new_contexts.append(template_contents)
        template_contents[i] = [reagent[0], reagent[1]*2]
        new_contexts.append(template_contents)
    return new_contexts

def reaction_yield(data_file_path, ret_time, smiles):
    pass

def isolation_yield(data_file_path, ret_time, smiles):
    pass

def read_shimadzu_pda_data(data_file: str):
    """
    Translate a data .txt file into a pandas dataframe

    Parameters
    ----------
    data_file : str
        path to .txt file.

    Returns
    -------
    dataframe with data
    """
    with open(data_file, 'r') as data_in:
        lines = data_in.readlines()
    for i, line in enumerate(lines):
        if len(line) > 200:
            break
    data = pd.read_csv(data_file, header=i-2, index_col=0)
    return data
            
def calc_amount(data_file: str, preds_file: str, smiles: str, ret_time: float):
    """
    calc_amount uses the raw data output from Shimadzu and a prediction of the UVvis absorbance and extinction
    coefficient to estimate the amount of material present in a peak, indicated by the retention time.
    Currently hard coded for 10 uL injection volume

    Parameters
    ----------
    data_file : str
        path to the data file (.txt).
    preds_file : str
        file with the chemprop predictions
    smiles : str
        smiles for the molecule of interest.
    ret_time : float
        retention time for the molecule of interest.

    Returns
    -------
    concentration of the product of interest in the sample

    """
    flow_rate = 1e-3/60
    inj_vol = 10e-6
    path_length = 1
    
    data = read_shimadzu_pda_data(data_file)
    wl = list(data.columns)
    wl = [float(i)/100 for i in wl]
    
    with open(preds_file, 'r') as preds_in:
        lines = preds_in.readlines()
    wl = {'O': 0.0, 'CC#N': 0.0}
    logE = {'O': 0.0, 'CC#N': 0.0}
    for line in lines:
        vals = line.split(',')
        if smiles in vals:
            wl[vals[1]] = float(vals[3])
            logE[vals[1]] = float(vals[2])
    # assume the extinction coef is the weighted average of the 'O' and 'CC#N' values based on solv composition
    y = min(max(95/4*(ret_time-0.5),5),100)/100
    ext_coef = 10**(y*logE['CC#N'] + (1-y)*logE['O'])
    
    
    # just out of curiosity, compare the predicted UVvis to the measured one at the retention time
    # get max wavelength with prominence of XXXX at MS retention time
    ms_rt_idx = np.where(data.index<ret_time)[0][-1]   # get the index of the MS retention time
    wl_peaks, _ = sp.find_peaks(data.iloc[ms_rt_idx], prominence=5e4)   # get the max wavelengths at that ret time
    if not wl_peaks.any():
        wl_peaks = np.array([1])
    pda_peaks, _ = sp.find_peaks(data[data.columns[wl_peaks[-1]]], prominence=1e4)   # get pda ret time from largest wl
    # get peak closest to the ms peak
    matched_peak = pda_peaks[np.argmin(abs(pda_peaks-ms_rt_idx))]
    # get peak bounds... used to integrate peak area
    peak_width = sp.peak_widths(data[data.columns[wl_peaks[-1]]], [matched_peak], rel_height=1.0)
    signal_to_integrate = data[data.columns[wl_peaks[-1]]].iloc[int(peak_width[2]):int(peak_width[3])]
    area = np.trapz(signal_to_integrate.values, signal_to_integrate.index)
    
    abs_units = area*1e-6
    
    conc = abs_units * flow_rate / (ext_coef * path_length * inj_vol)
    return conc


if not chemprop_server_ping():
    root = tk.Tk()
    prompt = "Open Anaconda powershell \ncd Desktop \nconda activate chemprop \npython chemprop_handler.py"
    QuickUI(root, title="Start Chemprop Server", dialog=prompt, buttons={}).run()

