# -*- coding: utf-8 -*-
"""
Created on 10/17/2022

Code that creates standard groups from csv files or pd DataBases and submits these to queues

@author: Jakob Dahl
"""
import pandas as pd
import numpy as np
import glob
import datetime
import os
from queue_workflows import queue_workflows
from groups_to_queues import generate_queue_document
import yaml
import warnings
warnings.simplefilter("default")


class ReactionSet():
    
    def __init__(self):
        self.database = pd.DataFrame()
        self.group = {}
        self.queue = ''

    def load(self, filename):
        try:
            self.database = pd.read_csv(filename, index_col=0)
        except Exception as E:
            print(E)
            warnings.warn('Could not read file: ' + filename)
        return self.database

    def set_db(self, dataframe):
        self.database = dataframe
        return self.database


    def generate_group(self, time_points, name, use_UV=True, Temperature=30, wellplate_type = '96 Well MediumWell', aliquot_solvent = 'dmso',
                       aliquot_volume=10, aliquot_target=200, total_volume=250, solvent='COc1ccccc1', intermediate_purification=False, separate_hplc_plate=False,
                       hplc_analysis=False, order=[3, 5, 1, 4, 2, 1], reagent_group_indices={'reactants': [0, 1], 'catalysts': [2, 3], 'reagents': [4, 5]},
                       additional_prep_details={'single_transfer': True}):
        # time points: List of numbers denoting time at which aliquot samples are taken, in seconds
        # name: name of campaign
        # volumes are in microliters (aliquot volume, aliquot target, total_volume)
        # concentrations are in mols/liter or mmol / mL
        
        if self.database.empty:
            raise ValueError('Database may not be empty!')
        
        if use_UV is True:
            optional_characterization_flag = True
            optional_characterization = {'characterization_type': 'uvvis','frequency': 'every_plate'}
        elif isinstance(use_UV, dict):
            optional_characterization = use_UV
            optional_characterization_flag = True
        else:
            optional_characterization = {}
            optional_characterization_flag = False

        if len(order) < len(self.database.columns):
            order = order + [1 for i in range(len(self.database.columns)-len(order))]

        reagents = list(self.database.columns)

        reaction_index_from_group = []
        for value in reagent_group_indices.values():
            for rxn in value:
                if rxn in reaction_index_from_group:
                    raise ValueError('Species at index ' + str(rxn) + ' was used more than once. A chemical species may not be classified as multiple species!' )
                else:
                    reaction_index_from_group.append(rxn)
        
        for cix in range(len(reagents)):
            if cix not in reaction_index_from_group:
                reagent_group_indices['reagents'].append(cix)
        print(reagents)
        print([group for group in reagent_group_indices.keys()])
        print([[rxn for rxn in reagent_group_indices[group]] for group in reagent_group_indices.keys()])
        print({group : [[reagents[rxn]] for rxn in reagent_group_indices[group]] for group in reagent_group_indices.keys()})
        print(dict({group : [[reagents[rxn],self.database.loc[list(self.database.index)[0],reagents[rxn]],
                                             order[rxn]] for rxn in reagent_group_indices[group]] for group in reagent_group_indices.keys()},**{'solvents': [[solvent, 1, 1]],'total_volume': total_volume}))
        #print([ix for ix in self.database.index])
        #print([self.database[:,ix] for ix in self.database.index])

        #print(self.database)
        reaction_grouping = {ix : dict({group : [[reagents[rxn],float(self.database.loc[ix,reagents[rxn]]),
                                             order[rxn]] for rxn in reagent_group_indices[group]  if not self.database.loc[ix,reagents[rxn]] == 0] for group in reagent_group_indices.keys()},**{'solvents': [[solvent, 1, 1]],'total_volume': total_volume}) for ix in self.database.index}
        
        for kn in reaction_grouping.keys():
            print(reaction_grouping[kn]['reactants'])
            if len(reaction_grouping[kn]['reactants']) == 0:
                print(reaction_grouping[kn]['reagents'][-1])
                reaction_grouping[kn]['reactants'].append(reaction_grouping[kn]['reagents'][-1])
                reaction_grouping[kn]['reagents'] = reaction_grouping[kn]['reagents'][:-1]

        self.group = {'campaign_name': name,
                         'workflow_type': 'dual_catalyst_exploration',
                         'wellplate_grouping':{
                             'reaction_plate_1':{
                                 'workflow_details':{
                                     'reaction_temperature' : Temperature,
                                     'wellplate_type': wellplate_type,
                                     'intermediate_purification': intermediate_purification,
                                     'separate_hplc_plate': separate_hplc_plate,
                                     'hplc_analysis': hplc_analysis,
                                     'hplc_processing': 'kinetics',
                                     'aliquot_time_series': time_points,
                                     'aliquot_details':{
                                         'solvent_to_use': aliquot_solvent,
                                         'aliquot_volume': aliquot_volume,
                                         'target_volume': aliquot_target,

                                     },
                                     'optional_characterization': optional_characterization,
                                     'additional_prep_details' : additional_prep_details,
                                     'liquid_class': 'Adjusted Low Vapor Free Dispense',
                                     'transfer_type': 'single_transfer'
                                 },
                                 'reaction_grouping': reaction_grouping

                             }
                         }
                         }
        if not optional_characterization_flag:
            del self.group['wellplate_grouping']['reaction_plate_1']['workflow_details']['optional_characterization']

        return self.group

    def generate_queue(self,output_directory = os.getcwd()):
        """
        if self.group == {}:
            t_name = glob.glob('./test_dual_cat_campaign*.yaml')
            if t_name == []:
                t_num = '001'
            else:
                t_nums = [int(n[-8:-5])]
                mt_num = max(t_nums)
                if mt_num >= 100:
                    t_num = str(mt_num)
                elif mt_num >= 10:
                    t_num = '0' + str(mt_num)
                else:
                    t_num = '00' + str(mt_num)
            self.generate_group([60 * i for i in [5,10,15,20,30,60,120,240,480]],'test_dual_cat_campaign' + t_num)
            warnings.warn('Generated group with standard arguments from generate_queue call. Specify arguments through calling generate_group directly if needed.')
        """
        print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
        print(output_directory)
        print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
        self.queue = generate_queue_document(output_directory, self.group, {})
        return self.queue


if __name__ == '__main__':
    RS = ReactionSet()
    RS.set_db(pd.DataFrame([[0.01,0.01,0.01,0.01,0.01,0.00],[0.02,0.03,0.04,0.05,0.06,0.07]],columns = ['A','B','C','D','E','F']))
    print(RS.database)
    #RS.generate_group([60 * i for i in [5,10,15,20,30,60,120,240,480]],'dual_cat_test')
    #RS.generate_queue()

    RS2 = ReactionSet()
    RS2.load('data_table_3.csv')
    RS2.generate_group([60 * i for i in [0,7.5,15,30,45,60,120,180,0]],'real_Cu_CN_test3',use_UV = False, order = [1,2,5,3,4],reagent_group_indices = {'reactants' : [2,3],'catalysts' : [1],'reagents' : [0,4]})
    RS2.generate_queue()
    print(RS2.group)
    print(RS2.queue)

    #RS3 = ReactionSet()
    #RS3.load('standards.csv')
    #RS3.generate_group([60 * i for i in [0,7.5,15]],'standards',use_UV = False, order = [1,1,1,1,1,2],reagent_group_indices = {'reactants' : [0,1,2,3,4],'reagents' : [5]},total_volume=100)
    #RS3.generate_queue()
