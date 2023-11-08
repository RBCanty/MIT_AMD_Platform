# -*- coding: utf-8 -*-
"""
Created on 11/02/2022

Queue generation using a csv file containing photochem reactions

@author: Nicola Hagger
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
import random
warnings.simplefilter("default")

class Reactions():

    def __init__(self):
        self.queue = ''
        self.list_of_wells = []
        alphabet = ['A','B','C','D','E','F','G','H']
        for j in range(8):
            for i in range(1,12+1):
                self.list_of_wells.append(alphabet[j]+str(i))
        self.random_list_of_wells = random.sample(self.list_of_wells,len(self.list_of_wells))

    def load_csv(self,filepath):
        self.database = pd.read_csv(filepath,sep=';')
        return self.database

    def generate_queue(self,output_directory = os.getcwd()):
        self.queue = generate_queue_document(output_directory,self.group,{})
        return self.queue

    def generate_group(self,name):
        self.name = name
        workflow = 'photochemistry'
        reaction_temperature = 20
        wellplate_type = '96 Well Microplate'
        total_volume = 200 #In micro-liters


        """
        Parameters to think about
        """

        led_power_setting = 0#150 #Value 0-255
        intermediate_purification = True
        separate_hplc_plate = False
        hplc_analysis = True
        hplc_processing = 'photochemistry'
        reaction_time = 60 * 120  #Seconds
        incubation_time = 60 * 10 #Seconds
        stabilization_time = 60 * 10 #Seconds


        mol_reagent = 10E-6   #10micro Mol
        mol_cat = 0.05*mol_reagent #5 Mol percent of reagents
        equivalents = 1

        reaction_grouping_random = {self.random_list_of_wells[ix]: dict({'reactants':[[self.database['First Reactant'][ix],float(mol_reagent),1],
                                        [self.database['Second Reactant'][ix],float(mol_reagent),1],
                                        [self.database['Third Reactant'][ix],float(mol_reagent),1]],
                                        'solvents':[[self.database['Solvent'][ix],equivalents,1]],
                                        'catalysts':[[self.database['Catalyst'][ix],float(mol_cat),2]],
                                        'total_volume': total_volume,}) for ix in self.database.index
                            }
        standard_reaction = {self.random_list_of_wells[-1]: {'reactants':[['C(F)(F)(F)S(=O)(=O)[O-].[Na+]',2E-5,1], ['COc1cc(OC)cc(OC)c1',0.5E-5,1]],
                'solvents': ['CS(C)=O.CC(=O)C(=O)C'],
                'total_volume': total_volume,}
             }
        reaction_grouping = {**reaction_grouping_random, **standard_reaction}


        self.group = {
                        'campaign_name': self.name,
                        'workflow_type': workflow,
                         'wellplate_grouping':{
                             'reaction_plate1':{
                                 'workflow_details': {
                                    'reaction_temperature': reaction_temperature,
                                    'reaction_time': reaction_time,
                                    'incubation_time': incubation_time,
                                    'stabilization_time': stabilization_time,
                                    'wellplate_type': wellplate_type,
                                    'intermediate_purification': intermediate_purification,
                                    'separate_hplc_plate': separate_hplc_plate,
                                    'hplc_analysis': hplc_analysis,
                                    'hplc_processing': hplc_processing,
                                    'led_power_setting': led_power_setting,
                                    'liquid_class': 'Adjusted water free dispense breakoff',
                                    'single_transfer': True},
                                'reaction_grouping': reaction_grouping
                             }
                         }}

    def generate_group_stock_sol_standards(self,name):

        self.name = name
        workflow = 'photochemistry'
        reaction_temperature = 20
        wellplate_type = '96 Well Microplate'
        total_volume = 200 #In micro-liters


        """
        Parameters to think about
        """

        led_power_setting = 0 #Value 0-255
        intermediate_purification = False
        separate_hplc_plate = False
        hplc_analysis = True
        hplc_processing = 'photochemistry'
        reaction_time = 60 * 120  #Seconds
        incubation_time = 60 * 10 #Seconds
        stabilization_time = 60 * 10 #Seconds


        mol_reagent = 10E-6   #10micro Mol
        equivalents = 1

        reaction_grouping_random = {self.list_of_wells[ix]: dict({'reactants':[[self.database['First Reactant'][ix],float(mol_reagent),1],
                                        [self.database['Second Reactant'][ix],float(mol_reagent),1],
                                        [self.database['Third Reactant'][ix],float(mol_reagent),1]],
                                        'solvents':[[self.database['Solvent'][ix],equivalents,1]],
                                        'catalysts':[[self.database['Catalyst'][ix],float(mol_cat),2]],
                                        'total_volume': total_volume,}) for ix in self.database.index
                            }
        standard_reaction = {self.random_list_of_wells[-1]: {'reactants':[['C(F)(F)(F)S(=O)(=O)[O-].[Na+]',2E-5,1], ['COc1cc(OC)cc(OC)c1',0.5E-5,1]],
                'solvents': ['CS(C)=O.CC(=O)C(=O)C'],
                'total_volume': total_volume,}
             }
        reaction_grouping = {**reaction_grouping, **standard_reaction}

        self.group = {
                        'campaign_name': self.name,
                        'workflow_type': workflow,
                         'wellplate_grouping':{
                             'reaction_plate1':{
                                 'workflow_details': {
                                    'reaction_temperature': reaction_temperature,
                                    'reaction_time': reaction_time,
                                    'incubation_time': incubation_time,
                                    'stabilization_time': stabilization_time,
                                    'wellplate_type': wellplate_type,
                                    'intermediate_purification': intermediate_purification,
                                    'separate_hplc_plate': separate_hplc_plate,
                                    'hplc_analysis': hplc_analysis,
                                    'hplc_processing': hplc_processing,
                                    'led_power_setting': led_power_setting,
                                    'liquid_class': 'Adjusted water free dispense breakoff',
                                    'single_transfer': True},
                                'reaction_grouping': reaction_grouping
                             }
                         }}

    def generate_group_multiple_wells(self,name):
        self.name = name
        workflow = 'photochemistry'
        reaction_temperature = 20
        wellplate_type = '96 Well Microplate'
        total_volume = 200 #In micro-liters

        list_of_multi_wells = ['A1','B1','C1','D1',
                    'A2','B2','C2','D2',
                    'A3','B3','C3','D3',
                    'A4','B4','C4','D4',
                    'A5','B5','C5','D5',
                    'A6','B6','C6','D6',
                    'A7','B7','C7','D7']


        """
        Parameters to think about
        """

        led_power_setting = 150 #Value 0-255
        intermediate_purification = True
        separate_hplc_plate = False
        hplc_analysis = True
        hplc_processing = 'photochemistry'
        reaction_time = 60 * 120  #Seconds
        incubation_time = 60 * 10 #Seconds
        stabilization_time = 60 * 10 #Seconds


        mol_reagent = 20E-6   #10micro Mol
        mol_reagent_inosine = 10E-6
        mol_cat = 0.05*mol_reagent #5 Mol percent of reagents
        mol_cat_inosine = 0.05*mol_reagent_inosine
        equivalents = 1
 
        num_wells_per_reaction = 4
        num_wells_per_reaction_inosine = 8
        num_reactions = 6


        reaction_grouping = {}
        index = 0
        for i in range(num_reactions):
            if self.database['First Reactant'][i] != 'O=c1nc[nH]c2c1ncn2[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O' and self.database['Second Reactant'][i] != 'O=c1nc[nH]c2c1ncn2[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O' and self.database['Third Reactant'][i] != 'O=c1nc[nH]c2c1ncn2[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O':
                for j in range(num_wells_per_reaction):
                    reaction_grouping[list_of_multi_wells[index]] = {'reactants':[[self.database['First Reactant'][i],float(mol_reagent),1],
                                            [self.database['Second Reactant'][i],float(mol_reagent),1],
                                            [self.database['Third Reactant'][i],float(mol_reagent),1]],
                                            'solvents':[[self.database['Solvent'][i],equivalents,1]],
                                            'catalysts':[[self.database['Catalyst'][i],float(mol_cat),2]],
                                            'total_volume': total_volume,}
                    index += 1
            else:
                for j in range(num_wells_per_reaction_inosine):
                    reaction_grouping[list_of_multi_wells[index]] = {'reactants':[[self.database['First Reactant'][i],float(mol_reagent_inosine),1],
                                            [self.database['Second Reactant'][i],float(mol_reagent_inosine),1],
                                            [self.database['Third Reactant'][i],float(mol_reagent_inosine),1]],
                                            'solvents':[[self.database['Solvent'][i],equivalents,1]],
                                            'catalysts':[[self.database['Catalyst'][i],float(mol_cat_inosine),2]],
                                            'total_volume': total_volume,}
                    index += 1

        standard_reaction = {'A8': {'reactants':[['C(F)(F)(F)S(=O)(=O)[O-].[Na+]',2E-5,1], ['COc1cc(OC)cc(OC)c1',0.5E-5,1]],
                'solvents': ['CS(C)=O.CC(=O)C(=O)C'],
                'total_volume': total_volume,}
             }
        reaction_grouping = {**reaction_grouping, **standard_reaction}

        self.group = {
                        'campaign_name': self.name,
                        'workflow_type': workflow,
                         'wellplate_grouping':{
                             'reaction_plate1':{
                                 'workflow_details': {
                                    'reaction_temperature': reaction_temperature,
                                    'reaction_time': reaction_time,
                                    'incubation_time': incubation_time,
                                    'stabilization_time': stabilization_time,
                                    'wellplate_type': wellplate_type,
                                    'intermediate_purification': intermediate_purification,
                                    'separate_hplc_plate': separate_hplc_plate,
                                    'hplc_analysis': hplc_analysis,
                                    'hplc_processing': hplc_processing,
                                    'led_power_setting': led_power_setting,
                                    'liquid_class': 'Adjusted water free dispense breakoff',
                                    'single_transfer': True},
                                'reaction_grouping': reaction_grouping
                             }
                         }}

        print(reaction_grouping)

    def large_scale_up(self,name):
        """
        This function can be used to scale up reactions, two reactions are needed, then each of them will be put in 47 wells
        Also the standard reaction will be added in the last well
        """

        self.name = name
        workflow = 'photochemistry'
        reaction_temperature = 20
        wellplate_type = '96 Well Microplate'
        total_volume = 200 #In micro-liters

        """
        Parameters to think about
        """

        led_power_setting = 150 #Value 0-255
        intermediate_purification = True
        separate_hplc_plate = False
        hplc_analysis = False
        hplc_processing = 'photochemistry'
        reaction_time = 60 * 120  #Seconds
        incubation_time = 60 * 10 #Seconds
        stabilization_time = 60 * 10 #Seconds


        mol_reagent = 10E-6   #10micro Mol
        mol_cat = 0.05*mol_reagent #5 Mol percent of reagents
        equivalents = 1
 
        num_wells_per_reaction = 47


        reaction_grouping = {}
        index = 0
    
        #First reaction:
        reactant_1_1 = 'OB(O)c1ccccc1'
        reactant_1_2 = 'O=c1nc[nH]c2c1ncn2[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O'
        reactant_1_3 = 'c1ccc(cc1)P(c2ccccc2)c3ccccc3'
        solvent_1 = 'CS(C)=O'
        cat_1 = 'CC1=CC(=C(C(=C1)C)C2=C3C=CC=CC3=[N+](C4=CC=CC=C42)C)C.[O-]Cl(=O)(=O)=O'
        for j in range(num_wells_per_reaction):
            reaction_grouping[self.list_of_wells[index]] = {'reactants':[[str(reactant_1_1),float(mol_reagent),1],
                                    [str(reactant_1_2),float(mol_reagent),1],
                                    [str(reactant_1_3),float(mol_reagent),1]],
                                    'solvents':[[str(solvent_1),equivalents,1]],
                                    'catalysts':[[str(cat_1),float(mol_cat),2]],
                                    'total_volume': total_volume,}
            index += 1
        index += 1 #Well D12 will be empty

        #reactant_2_1 = 'Oc1ccc(/C=C/c2cc(O)cc(O)c2)cc1'
        #reactant_2_2 = 'Cc1ccc(C(=O)C(C)CN2CCCCC2)cc1.Cl'
        #reactant_2_3 = 'Oc1cccc(O)c1'
        #solvent_2 = 'CS(C)=O'
        #cat_2 = 'CC1=CC(=C(C(=C1)C)C2=C3C=CC=CC3=[N+](C4=CC=CC=C42)C)C.[O-]Cl(=O)(=O)=O'
        reactant_2_1 = reactant_1_1
        reactant_2_2 = reactant_1_2
        reactant_2_3 = reactant_1_3
        solvent_2 = solvent_1
        cat_2 = cat_1
        for j in range(num_wells_per_reaction):
            reaction_grouping[self.list_of_wells[index]] = {'reactants':[[str(reactant_2_1),float(mol_reagent),1],
                                    [str(reactant_2_2),float(mol_reagent),1],
                                    [str(reactant_2_3),float(mol_reagent),1]],
                                    'solvents':[[str(solvent_2),equivalents,1]],
                                    'catalysts':[[str(cat_2),float(mol_cat),2]],
                                    'total_volume': total_volume,}
            index += 1

        standard_reaction = {'H12': {'reactants':[['C(F)(F)(F)S(=O)(=O)[O-].[Na+]',2E-5,1], ['COc1cc(OC)cc(OC)c1',0.5E-5,1]],
                'solvents': ['CS(C)=O.CC(=O)C(=O)C'],
                'total_volume': total_volume,}
             }
        reaction_grouping = {**reaction_grouping, **standard_reaction}

        self.group = {
                        'campaign_name': self.name,
                        'workflow_type': workflow,
                         'wellplate_grouping':{
                             'reaction_plate1':{
                                 'workflow_details': {
                                    'reaction_temperature': reaction_temperature,
                                    'reaction_time': reaction_time,
                                    'incubation_time': incubation_time,
                                    'stabilization_time': stabilization_time,
                                    'wellplate_type': wellplate_type,
                                    'intermediate_purification': intermediate_purification,
                                    'separate_hplc_plate': separate_hplc_plate,
                                    'hplc_analysis': hplc_analysis,
                                    'hplc_processing': hplc_processing,
                                    'led_power_setting': led_power_setting,
                                    'liquid_class': 'Adjusted water free dispense breakoff',
                                    'single_transfer': True},
                                'reaction_grouping': reaction_grouping
                             }
                         }}

    def control_reactions(self,name):

        self.name = name
        workflow = 'photochemistry'
        reaction_temperature = 20
        wellplate_type = '96 Well Microplate'
        total_volume = 200 #In micro-liters

        led_power_setting = 150 #Value 0-255
        intermediate_purification = True
        separate_hplc_plate = False
        hplc_analysis = True
        hplc_processing = 'photochemistry'
        reaction_time = 60 * 120  #Seconds
        incubation_time = 60 * 10 #Seconds
        stabilization_time = 60 * 10 #Seconds

        list_of_wells = ['A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12','B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11','B12','C1','C2','C3','C4','C5','C6','C7']
        random.shuffle(list_of_wells)
        print(list_of_wells)

        mol_reagent = 10E-6   #10micro Mol
        mol_cat = 0.05*mol_reagent #5 Mol percent of reagents
        equivalents = 1

        index = 0

        ['C(F)(F)(F)S(=O)(=O)[O-].[Na+]',2E-5,1]

        reactant_1_1 = 'COC(=O)/C=C/C(=O)OC'
        equiv_1_1 = 1
        reactant_1_2 = 'C1C=CN(C=C1C(=O)N)CC2=CC=CC=C2'
        equiv_1_2 = 2
        reactant_1_3 = 'c1ccncc1'
        equiv_1_3 = 1
        solvent_1 = 'CN(C)C=O'
        cat_1 = 'C1=CC=NC(=C1)C2=CC=CC=N2.C1=CC=NC(=C1)C2=CC=CC=N2.C1=CC=NC(=C1)C2=CC=CC=N2.Cl[Ru]Cl'
        reaction_1 = {}
        for i in range(5):
            reaction_1[list_of_wells[index]] = {'reactants':[[reactant_1_1,mol_reagent*equiv_1_1,1],[reactant_1_2,mol_reagent*equiv_1_2,1],[reactant_1_3,mol_reagent*equiv_1_3,1]],
                'solvents':[[solvent_1,equivalents,1]],
                'catalysts':[[cat_1,float(mol_cat),2]],
                'total_volume': total_volume,}
            index += 1
        
        reactant_2_1 = 'C1=CC=C(C=C1)CBr'
        equiv_2_1 = 5
        reactant_2_2 = 'C1C=CN(C=C1C(=O)N)CC2=CC=CC=C2'
        equiv_2_2 = 0.5
        reactant_2_3 = 'c1ccncc1'
        equiv_2_3 = 1
        solvent_2 = 'CS(C)=O.CC#N'
        cat_2 = 'C1=CC=NC(=C1)C2=CC=CC=N2.C1=CC=NC(=C1)C2=CC=CC=N2.C1=CC=NC(=C1)C2=CC=CC=N2.Cl[Ru]Cl'
        reaction_2 = {}
        for i in range(5):
            reaction_2[list_of_wells[index]] = {'reactants':[[reactant_2_1,mol_reagent*equiv_2_1,1],[reactant_2_2,mol_reagent*equiv_2_2,1],[reactant_2_3,mol_reagent*equiv_2_3,1]],
                'solvents':[[solvent_2,equivalents,1]],
                'catalysts':[[cat_2,float(mol_cat),2]],
                'total_volume': total_volume,}
            index += 1
        
        reactant_3_1 = 'COc1cc(OC)cc(OC)c1'
        equiv_3_1 = 1
        reactant_3_2 = 'c1ccncc1'
        equiv_3_2 = 2
        reactant_3_3 = 'C(=O)(C(F)(F)F)OC(=O)C(F)(F)F'
        equiv_3_3 = 2.1
        solvent_3 = 'CN(C)C=O.CC#N'
        cat_3 = 'C1=CC=NC(=C1)C2=CC=CC=N2.C1=CC=NC(=C1)C2=CC=CC=N2.C1=CC=NC(=C1)C2=CC=CC=N2.Cl[Ru]Cl'
        reaction_3 = {}
        for i in range(5):
            reaction_3[list_of_wells[index]] = {'reactants':[[reactant_3_1,mol_reagent*equiv_3_1,1],[reactant_3_2,mol_reagent*equiv_3_2,1],[reactant_3_3,mol_reagent*equiv_3_3,1]],
                'solvents':[[solvent_3,equivalents,1]],
                'total_volume': total_volume-10,}
            index += 1
        #'catalysts':[[cat_3,float(mol_cat),2]], => Has to be manually added

        reactant_4_1 = 'Sc1ccccc1'
        equiv_4_1 = 1
        reactant_4_2 = 'CCCCCCC=C'
        equiv_4_2 = 1
        solvent_4 = 'CS(C)=O.CC#N'
        cat_4 = 'CC1=CC(=C(C(=C1)C)C2=C3C=CC=CC3=[N+](C4=CC=CC=C42)C)C.[O-]Cl(=O)(=O)=O'
        reaction_4 = {}
        for i in range(5):
            reaction_4[list_of_wells[index]] = {'reactants':[[reactant_4_1,mol_reagent*equiv_4_1,1],[reactant_4_2,mol_reagent*equiv_4_2,1]],
                'solvents':[[solvent_4,equivalents,1]],
                'catalysts':[[cat_4,float(mol_cat),2]],
                'total_volume': total_volume,}
            index += 1

        reactant_5_1 = 'C1=CC=C(C(=C1)C=O)Br'
        equiv_5_1 = 1
        reactant_5_2 = 'C1=CC=C(C(=C1)N)N'
        equiv_5_2 = 1
        solvent_5 = 'CC1=CC=CC=C1'
        cat_5 = 'C1=CC=C2C(=C1)C(=O)OC23C4=C(C=C(C=C4)O)OC5=C3C=CC(=C5)O'
        reaction_5 = {}
        for i in range(5):
            reaction_5[list_of_wells[index]] = {'reactants':[[reactant_5_1,mol_reagent*equiv_5_1,1],[reactant_5_2,mol_reagent*equiv_5_2,1]],
                'solvents':[[solvent_5,equivalents,1]],
                'catalysts':[[cat_5,float(mol_cat),2]],
                'total_volume': total_volume,}
            index += 1

        reactant_6_1 = 'CCOC(=O)C(C(=O)OCC)Br'
        equiv_6_1 = 2
        reactant_6_2 = 'C1=COC=C1'
        equiv_6_2 = 1
        reactant_6_3 = 'COC1=CC=C(C=C1)N(C2=CC=CC=C2)C3=CC=CC=C3'
        equiv_6_3 = 2
        solvent_6 = 'CN(C)C=O'
        cat_6 = 'C1=CC=NC(=C1)C2=CC=CC=N2.C1=CC=NC(=C1)C2=CC=CC=N2.C1=CC=NC(=C1)C2=CC=CC=N2.Cl[Ru]Cl'
        reaction_6 = {}
        for i in range(5):
            reaction_6[list_of_wells[index]] = {'reactants':[[reactant_6_1,mol_reagent*equiv_6_1,1],[reactant_6_2,mol_reagent*equiv_6_2,1],[reactant_6_3,mol_reagent*equiv_6_3,1]],
                'solvents':[[solvent_6,equivalents,1]],
                'catalysts':[[cat_6,float(mol_cat),2]],
                'total_volume': total_volume,}
            index += 1

        standard_reaction = {list_of_wells[index]: {'reactants':[['C(F)(F)(F)S(=O)(=O)[O-].[Na+]',2E-5,1], ['COc1cc(OC)cc(OC)c1',0.5E-5,1]],
                'solvents': ['CS(C)=O.CC(=O)C(=O)C'],
                'total_volume': total_volume,}
             }


        reaction_grouping = {**reaction_1, **reaction_2, **reaction_3, **reaction_4, **reaction_5, **reaction_6, **standard_reaction}

        self.group = {
                        'campaign_name': self.name,
                        'workflow_type': workflow,
                         'wellplate_grouping':{
                             'reaction_plate1':{
                                 'workflow_details': {
                                    'reaction_temperature': reaction_temperature,
                                    'reaction_time': reaction_time,
                                    'incubation_time': incubation_time,
                                    'stabilization_time': stabilization_time,
                                    'wellplate_type': wellplate_type,
                                    'intermediate_purification': intermediate_purification,
                                    'separate_hplc_plate': separate_hplc_plate,
                                    'hplc_analysis': hplc_analysis,
                                    'hplc_processing': hplc_processing,
                                    'led_power_setting': led_power_setting,
                                    'liquid_class': 'Adjusted water free dispense breakoff',
                                    'single_transfer': True},
                                'reaction_grouping': reaction_grouping
                             }
                         }}

    def generate_queue(self,output_directory = os.getcwd()):
        return_statement = generate_queue_document(output_directory,self.group,{})
        print(return_statement)
        print('DONE, generated a queue documente and stored it to '+str(output_directory))



if __name__ == '__main__':
    filepath_csv = r'C:\Users\nicol\Downloads\Dark_Control_2.csv'
    output_directory = os.getcwd()+'\photochem_queues'
    output_directory = r'C:\Users\nicol\Documents\GitHub\queue_generation\photochem_queues'
    queue_name = 'large_scale_up_02'


    rs = Reactions()
    #rs.load_csv(filepath_csv)
    #rs.generate_group(queue_name)
    #rs.generate_group_stock_sol_standards(queue_name)
    #rs.generate_group_multiple_wells(queue_name)
    rs.large_scale_up(queue_name)
    #rs.control_reactions(queue_name)
    rs.generate_queue(output_directory)
