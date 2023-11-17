# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 17:24:08 2020

@author: Brent

This script takes the Carrier.cfg file and converts it to a JSON dictionary to 
make it easier to deal with in scripting. Because of this, the script needs to
update that carrier_labware definition file in order to properly work with 
EVOWare scripting.

But this script is infrequently executed, only when there is a problem with the
carrier_labware definitions or when new things are added to the CFG database

Alot of the information about the carriers define the default positions of all
of the sites on the carrier, important to EVOWare training but not to us really

The information about the labware tracks the allowed carriers and some 
dimensional information that EVO uses to prepopulate default values, we keep
the allowed carriers for (future) error checking purposes

"""

import json, sys

filepath = r'.\Carrier.cfg'
final_dict = {'carriers':{}, 'labware':{}}

number_allowed = 0
with open(filepath, 'r') as infile:
    for line in infile:
        if line.startswith('13'):
            temp_split = line.strip('\n').split(';')
            key_number = int(temp_split[2].split('/')[0])
            final_dict['carriers'][key_number] = {'name': temp_split[1], 'positions': int(temp_split[-3]), 'line': line.strip('\n')}
        elif line.startswith('15'):
            temp_split = line.strip('\n').split(';')
            name = temp_split[1]
            number_allowed = int(temp_split[-2])
            dimension_split = line.split(';')[3].split('/')
            #print(dimension_split)
            final_dict['labware'][name] = {'name': name, 'allowed_carriers': [], 'line': line.strip('\n'),
                                           'nrows': int(dimension_split[1]), 'ncols': int(dimension_split[0])}
            continue
        if number_allowed > 0:
            temp_split = line.strip('\n').split(';')
            number_allowed -= 1
            if temp_split[1] == '0/0/0/0/0':
                continue
            final_dict['labware'][name]['allowed_carriers'].append(int(temp_split[1]))
            
with open(r'.\carriers_and_labware.json', 'w') as jsonfile:
    json.dump(final_dict, jsonfile)
    
