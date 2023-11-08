"""
Code to generate a YAML Well-Position Configuration File for the Bruker HTS-XT
Original implementation was for a 96-well silicon plate based off of Bruker Plates
Version 0.0.1
Updated on 12/7/2019
@author: brent
"""
import os
import yaml
import string
import sys

letters= list(string.ascii_uppercase)

os.chdir('C:\\AMD_FTIR_v2\\')
microns_per_step = 7.2
reference_position = [int(93240+7000//7.2),int(19123+1500//7.2)]
microns_between_wells = 12000
distance = microns_between_wells

array = [[[] for i in range(9)] for j in range(6)]
well_positions = {}
for i in range(0,6):
    x = reference_position[0] - i*distance
    for j in range(0,9):
        y = reference_position[1] + j*distance
        array[i][j] = [x,y]
        well_positions[letters[i]+str(j+1)] = [x,y]
with open(os.getcwd() + '\\si_54_well_positions.yaml','w') as yaml_file:
    yaml.dump(well_positions, yaml_file)



sys.exit(-1)
microns_per_step = 7.2
#reference_position = [93240,19123] #Pre-Kick position
reference_position = [int(93240+8000//7.2),19123]
microns_between_wells = 9000
distance = microns_between_wells

array = [[[] for i in range(12)] for j in range(8)]
well_positions = {}
for i in range(0,8):
    x = reference_position[0] - i*distance
    for j in range(0,12):
        y = reference_position[1] + j*distance
        array[i][j] = [x,y]
        well_positions[letters[i]+str(j+1)] = [x,y]
with open(os.getcwd() + '\\si_well_positions.yaml','w') as yaml_file:
    yaml.dump(well_positions, yaml_file)

