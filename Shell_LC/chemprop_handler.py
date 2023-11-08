#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 12:53:39 2022

Testing HTTP server to handle chemprop predictions run remotely

@author: mattmcdonald
"""

import json
import csv
import os
import subprocess
import threading
from http.server import BaseHTTPRequestHandler, HTTPServer

model_dir = r'C:\Shimadzu\Data\Automated\Chemprop_models\model_03222022\fold_2'

class handler(BaseHTTPRequestHandler):
    def do_GET(self):
        self.send_response(200)
        self.send_header('Content-type','application/json')
        self.end_headers()
        message = "Server is running"
        self.wfile.write(message.encode(encoding='utf_8'))

    def do_POST(self):
        content_len = int(self.headers.get('Content-Length'))
        post_body = self.rfile.read(content_len)
        
        info_in = json.loads(post_body.decode(encoding='utf-8'))
        if info_in == "ping":
            info_out = "chemprop server ping success"
        else:
            info_out = log_e_chemprop_predict(info_in)
        json_out = json.dumps(info_out)
        self.send_response(200)
        self.send_header('Content-type','application/json')
        self.end_headers()
        self.wfile.write(json_out.encode(encoding='utf-8'))
        
            #self.send_response(400)
            #self.send_header('Content-type','application/json')
            #self.end_headers()
              
# update these to match LC computer
log_e_to_calc_file = 'log_e_to_calc.csv'
log_e_preds_file = 'log_e_predictions.csv'

def log_e_chemprop_predict(wellplate: dict, solvents=None):
    """
    Run chemprop log_e predictions on a wellplate. Dump results into file_name. Can also give solvents, acceptable
    solvents are (currently) water ('O'), acetonitrile ('CC#N'), and dichloromethane ('ClCCl'). Also checks to see
    that all the reagents are catalogued for MOCCA

    Parameters
    ----------
    smiles : dict
        contents dictionary
    file_name : str
        file_name to be used to generate and retrieve predictions
    solvents : list, optional
        list of solvents to predict in. The default is ['O','CC#N'].

    Returns
    -------
    The name of the file that the chemprop predictions will be available in
    """

    if solvents is None:
        solvents = ['O', 'CC#N']

    def chemprop_thread_func(command):
        process = subprocess.run(command, shell = True, capture_output = True)    
        formatted_stderr = []
        for line in process.stderr.decode().split('\r'):
            if 'WARNING' in line:
                continue
            clean_line = line.replace('\r', '').strip('\n')
            if clean_line.strip() == '':
                continue
            if '#' not in clean_line:
                continue
            formatted_stderr.append(clean_line)
        print(process.returncode)
        print('Finishing!: %s' % command)
    
    # generate csv file to be fed to chemprop
    to_calc_file = os.path.join(wellplate['path'], log_e_to_calc_file)
    preds_file = os.path.join(wellplate['path'], log_e_preds_file)
    mol_solv_list = []
    for well in wellplate['contents'].values():
        for solv in solvents:
            mol_solv = [well['target_product'][0][0], solv]
            if mol_solv not in mol_solv_list:
                mol_solv_list.append(mol_solv)
    header = ['SMILES', 'solvent']
    if not os.path.exists(wellplate['path']):
        os.makedirs(wellplate['path'])
    with open(to_calc_file, 'w', newline='', encoding='utf-8') as new_csv:
        write = csv.writer(new_csv)
        write.writerow(header)
        write.writerows(mol_solv_list)
        
    # generate chemprop command string
    command_l1 = f"chemprop_predict --test_path '{to_calc_file}' --preds_path '{preds_file}' --checkpoint_dir '{model_dir}'"
    command_l2 = " --ensemble_variance --number_of_molecules 2"
    command = command_l1 + command_l2
    chemprop_thread = threading.Thread(target=chemprop_thread_func, args=(command,), daemon=True)
    chemprop_thread.start()
    
    #with open(preds_file, 'w', newline='', encoding='utf-8') as new_csv:
    #    write = csv.writer(new_csv)
    #    write.writerow(['Results pending...'])
        
    return preds_file
            
print("started server")
with HTTPServer(('', 40888), handler) as server:
    server.serve_forever()