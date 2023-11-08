#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 11:01:00 2021

@author: mattmcdonald

Shimadzu_API helper functions.
Functions here do not require access to the Shimadzu_obj.
"""

import numpy as np
from pyteomics import mzxml
import database_interface as dbi

def data_analysis_helper(target_sample):
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
    # look for common adducts [+proton, -proton, +ACN+H, +Na, -formate]... 
    # can add more (+K, e.g.) but increases false pos.
    mws = np.array([mw + 1, 1 - mw, mw + 42, mw + 23, 1 - (46 + mw)])

    # need to sleep to allow Shimadzu to populate the mzXML file
    file_name = target_sample[0][:-4] + ".mzXML"
    experiment = mzxml.read(file_name)
    spec = experiment.next()
    # hard coded to a 6 minute method with 1 Hz scan rate
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
    best_mz = 0
    # indentify peaks based on N consecutive points > M*sigma baseline
    n_con_pts = 5
    m_sigma = 3
    for mz in mws:
        if mz > 0:
            mic_low = np.where(mzs > (mz-0.5))[0][0]
            mic_high = np.where(mzs < (mz+0.5))[0][-1]
            mic_pos = np.sum(data_pos[mic_low:mic_high, :], 0)
            mic_pos = mic_pos/max(mic_pos)
            mu, sigma = baseline_helper(mic_pos, 10, 2)
            idx = np.array(np.where(mic_pos > (mu+m_sigma*sigma)))
            for i in range(0, n_con_pts-1):
                idx = idx[np.diff(idx, append=-1) == 1]
            if idx.size:
                if max(mic_pos) > max_intensity:
                    ret_time = rt_pos[idx[0]]
                    max_intensity = max(mic_pos)
                    best_mz = mz
        else:
            mic_low = np.where(mzs > (abs(mz)-0.5))[0][0]
            mic_high = np.where(mzs < (abs(mz)+0.5))[0][-1]
            mic_neg = np.sum(data_neg[mic_low:mic_high, :], 0)
            mic_neg = mic_neg/max(mic_neg)
            mu, sigma = baseline_helper(mic_neg, 10, 2)
            idx = np.array(np.where(mic_neg > (mu+m_sigma*sigma)))
            for i in range(0, n_con_pts-1):
                idx = idx[np.diff(idx, append=-1) == 1]
            if idx.size:
                if max(mic_neg) > max_intensity:
                    ret_time = rt_neg[idx[0]]
                    max_intensity = max(mic_neg)
                    best_mz = mz
                    
    return ret_time, best_mz

def FRC_plate_request_helper(container, queue_id, operation_id):
    print("use Aceso.py version please")
    lpx_op = {'operation': 'Ss.request',
              'agent': 'Ss',
              'completed': 'no',
              'container': container,
              'end_time': None,
              'start_time': None,
              'time_est': 20,
              'details': {'plate_type': '96 Well DeepWell'}}
    ra_op = {'operation': 'move_wellplate',
             'agent': 'Ra',
             'completed': 'no',
             'container': container,
             'end_time': None,
             'start_time': None,
             'time_est': 20,
             'details': {'target_destination': 'fraction_collector',
                         'scan_barcode': 'no'}}  # 'yes' if using barcode scanner
    dbi.insert_queue_steps(queue_id, int(operation_id)-1, [lpx_op, ra_op])