# -*- coding: utf-8 -*-
"""
Created on Sat Dec 14 14:52:06 2019

@author: brent
"""

import pubchempy
import requests
import re
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from data_entry import update

def pubchem_json_lookup(url_molecule):
    url_base = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/'
    url_tail = '/XML/?response_type=display'
    url_total = url_base + str(url_molecule) + url_tail
    web_page = requests.get(url_total)
    return web_page.text

def quick_gen(sublayerlist,name):
    return next((index for (index, d) in enumerate(sublayerlist) if d["TOCHeading"] == name), None)

def record_title(xmldata):
    start_index = xmldata.find('<RecordTitle>')
    end_index = xmldata.find('</RecordTitle>')
    return xmldata[start_index:end_index].split('>')[1].replace("&apos;","'")

def exp_logp_pubchem(xmldata):
    start_index = xmldata.find('<TOCHeading>Octanol/Water Partition Coefficient</TOCHeading>')
    end_index = xmldata[start_index::].find('</Section>')
    value_start_idx =[m.start() for m in re.finditer('<Information>',xmldata[start_index:end_index+start_index])]
    value_end_idx =[m.start() for m in re.finditer('</Information>',xmldata[start_index:end_index+start_index])]
    values_to_keep = []
    references_to_keep = []
    for i in range(0,len(value_start_idx)):
        string_check = xmldata[start_index+value_start_idx[i]:start_index+value_end_idx[i]]
        if '<Reference>' in string_check:
            references_to_keep.append(string_check.split('<Reference>')[1].split('</Reference>')[0])
        else:
            continue
        #if '-' in string_check:
        #    continue
        if '<Number>' in string_check:
            values_to_keep.append(string_check.split('<Number>')[1].split('</Number>')[0])
        elif '<String>' in string_check:
            result = re.sub(r" ?\([^)]+\)", "",string_check.split('<String>')[1].split('</String>')[0])
            values_to_keep.append(result.replace('log Kow = ','').replace('log Kow= ',''))
    return values_to_keep, references_to_keep

def f_to_c_converter(value):
    return (value-32)*5.0/9.0

def melting_point_pubchem(xmldata):
    start_index = xmldata.find('<TOCHeading>Melting Point</TOCHeading>')
    end_index = xmldata[start_index::].find('</Section>')
    value_start_idx =[m.start() for m in re.finditer('<Value>',xmldata[start_index:end_index+start_index])]
    value_end_idx =[m.start() for m in re.finditer('</Value>',xmldata[start_index:end_index+start_index])]
    values_to_keep = []
    for i in range(0,len(value_start_idx)):
        string_check = xmldata[start_index+value_start_idx[i]:start_index+value_end_idx[i]]
        if '<Number>' in string_check:
            try:
                first_value = string_check.split('<Number>')[1].split('</Number>')[0]
                #second_value = string_check.split('<Unit>')[1].split('</Unit>')[0].replace('Â','')
                values_to_keep.append(first_value)
            except IndexError:
                first_value = string_check.split('<Number>')[1].split('</Number>')[0]
                second_value = ''
                values_to_keep.append(first_value + '' + second_value)
        elif '<String>' in string_check:
            string_value = string_check.split('<String>')[1].split('</String>')[0].replace(' Â',' ').replace('Â',' ')
            try:
                if '-' in string_value:
                    continue
                elif 'greater than' in string_value:
                    continue
                elif '°F' in string_value:
                    string_value = "{:.2f}".format((f_to_c_converter(float(re.sub("\D+\.","",string_value.split(' °F')[0])))))
                    values_to_keep.append(string_value)
                elif '°C' in string_value:
                    string_value = string_value.split(' °C')[0]
                    values_to_keep.append(string_value)
            except ValueError:
                continue
    return values_to_keep

def boiling_point_pubchem(xmldata):
    start_index = xmldata.find('<TOCHeading>Boiling Point</TOCHeading>')
    end_index = xmldata[start_index::].find('</Section>')
    value_start_idx =[m.start() for m in re.finditer('<Value>',xmldata[start_index:end_index+start_index])]
    value_end_idx =[m.start() for m in re.finditer('</Value>',xmldata[start_index:end_index+start_index])]
    values_to_keep = []
    for i in range(0,len(value_start_idx)):
        string_check = xmldata[start_index+value_start_idx[i]:start_index+value_end_idx[i]]
        if '<Number>' in string_check:
            try:
                first_value = string_check.split('<Number>')[1].split('</Number>')[0]
                #second_value = string_check.split('<Unit>')[1].split('</Unit>')[0].replace('Â','')
                values_to_keep.append(first_value)
            except IndexError:
                first_value = string_check.split('<Number>')[1].split('</Number>')[0]
                second_value = ''
                values_to_keep.append(first_value + '' + second_value)
        elif '<String>' in string_check:
            string_value = string_check.split('<String>')[1].split('</String>')[0].replace(' Â',' ').replace('Â',' ')
            if '-' in string_value:
                continue
            elif 'greater than' in string_value:
                    continue
            elif '°F' in string_value:
                string_value = "{:.2f}".format((f_to_c_converter(float(re.sub("\D+\.","",string_value.split(' °F')[0])))))
                values_to_keep.append(string_value)
            elif '°C' in string_value:
                string_value = string_value.split(' °C')[0]
                values_to_keep.append(string_value)
            #values_to_keep.append(string_check.split('<String>')[1].split('</String>')[0].replace(' Â',' ').replace('Â',' '))
    return values_to_keep

def cas_number(xmldata):
    start_index = xmldata.find('<TOCHeading>CAS</TOCHeading>')
    end_index = xmldata[start_index::].find('</String>')
    if start_index != -1:
        string_check = xmldata[start_index:end_index+start_index].split("<String>")[1]
    else:
        string_check = []
    return string_check

def parent_compound(xmldata):
    start_index = xmldata.find('<TOCHeading>Parent Compound</TOCHeading>')
    end_index = xmldata[start_index::].find('</Number>') 
    if start_index != -1:
        string_check = xmldata[start_index:end_index+start_index].split("<Number>")[1]
    else:
        string_check = []
    return string_check

def density_grabber(xmldata):
    start_index = xmldata.find('<TOCHeading>Density</TOCHeading>')
    end_index = xmldata[start_index::].find('</Section>')
    value_start_idx =[m.start() for m in re.finditer('<Value>',xmldata[start_index:end_index+start_index])]
    value_end_idx =[m.start() for m in re.finditer('</Value>',xmldata[start_index:end_index+start_index])]
    values_to_keep = []
    for i in range(0,len(value_start_idx)):
        string_check = xmldata[start_index+value_start_idx[i]:start_index+value_end_idx[i]]
        if '<Number>' in string_check:
            print(string_check)
        elif '<String>' in string_check:
            string_value = string_check.split('<String>')[1].split('</String>')[0].replace(' Â',' ').replace('Â',' ')
            if 'at ' in string_value:
                continue
            elif '(water = 1)' in string_value:
                try:
                    string_value = string_value.split(':')[1].strip(' ')
                    values_to_keep.append(string_value)
                except IndexError:
                    continue
            else:
                values_to_keep.append(string_value)
    return values_to_keep

def vapor_pressure_grabber(xmldata):
    start_index = xmldata.find('<TOCHeading>Vapor Pressure</TOCHeading>')
    end_index = xmldata[start_index::].find('</Section>')
    value_start_idx =[m.start() for m in re.finditer('<Value>',xmldata[start_index:end_index+start_index])]
    value_end_idx =[m.start() for m in re.finditer('</Value>',xmldata[start_index:end_index+start_index])]
    values_to_keep = []
    for i in range(0,len(value_start_idx)):
        string_check = xmldata[start_index+value_start_idx[i]:start_index+value_end_idx[i]]
        if '<String>' in string_check:
            string_value = string_check.split('<String>')[1].split('</String>')[0].replace(' Â',' ').replace('Â',' ')
            if 'at ' in string_value:
                continue
            elif '°' in string_value:
                continue
            else:
                values_to_keep.append(string_value)
    return values_to_keep

def refractive_index_grabber(xmldata):
    start_index = xmldata.find('<TOCHeading>Refractive Index</TOCHeading>')
    end_index = xmldata[start_index::].find('</Section>')
    value_start_idx =[m.start() for m in re.finditer('<Value>',xmldata[start_index:end_index+start_index])]
    value_end_idx =[m.start() for m in re.finditer('</Value>',xmldata[start_index:end_index+start_index])]
    values_to_keep = []
    for i in range(0,len(value_start_idx)):
        string_check = xmldata[start_index+value_start_idx[i]:start_index+value_end_idx[i]]
        if '<String>' in string_check:
            string_value = string_check.split('<String>')[1].split('</String>')[0].replace(' Â',' ').replace('Â',' ').replace('°','')
            values_to_keep.append(string_value.strip('Index of refraction: ').strip('°'))
    return values_to_keep

def mesh_entry_grabber(xmldata):
    start_index = xmldata.find('<TOCHeading>MeSH Entry Terms</TOCHeading>')
    end_index = xmldata[start_index::].find('</Section>')
    value_start_idx =[m.start() for m in re.finditer('<String>',xmldata[start_index:end_index+start_index])]
    value_end_idx =[m.start() for m in re.finditer('</String>',xmldata[start_index:end_index+start_index])]
    values_to_keep = []
    for i in range(0,len(value_start_idx)):
        string_check = xmldata[start_index+value_start_idx[i]:start_index+value_end_idx[i]]
        string_value = string_check.split('<String>')[1]
        if ',' in string_value:
            string_value = " ".join(string_value.split(", ")[::-1])
        if string_value not in values_to_keep:
            values_to_keep.append(string_value)
    return values_to_keep

def get_pubchem_data(incoming_string):
    pypubdata = {}
    molecule = pubchempy.get_compounds(incoming_string,'inchikey')
    molecule_cid = molecule[0].cid
    if molecule_cid == None:
        return -1
    outcoming_xml_data = pubchem_json_lookup(molecule_cid)
    cas_num = cas_number(outcoming_xml_data)
    pypubdata['basic_properties'] = {}
    if len(cas_num) != 0:
        pypubdata['basic_properties']['compound_cas'] = cas_num
    pypubdata['basic_properties']['compound_iupac_name'] = molecule[0].iupac_name
    pypubdata['basic_properties']['other_identifiers'] = [record_title(outcoming_xml_data)]
    #mesh_entries = mesh_entry_grabber(outcoming_xml_data)
    #if len(mesh_entries) != 0:
    #    pypubdata['basic_properties']['other_identifiers'] = mesh_entries
    temp = {}
    temp['phys_chem_props'] = {}
    meltings = melting_point_pubchem(outcoming_xml_data)
    if len(meltings) != 0:
        temp['phys_chem_props']['melting_point'] = meltings
    boilings = boiling_point_pubchem(outcoming_xml_data)
    if len(boilings) != 0: 
        temp['phys_chem_props']['boiling_point'] = boilings 
    density = density_grabber(outcoming_xml_data)
    if len(density) != 0:
        temp['phys_chem_props']['density'] = density
    vapor_pressure = vapor_pressure_grabber(outcoming_xml_data)
    if len(vapor_pressure) != 0:
        temp['phys_chem_props']['vapor_pressure'] = vapor_pressure
    refractive_index = refractive_index_grabber(outcoming_xml_data)
    if len(refractive_index) != 0:
        temp['phys_chem_props']['refractive_index'] = refractive_index
    parent = parent_compound(outcoming_xml_data)
    #print(pypubdata)
    if any(temp['phys_chem_props']) != False:
        update(pypubdata, temp)
    #print(pypubdata)
    if len(parent) != 0:
        parent_molecule = pubchempy.get_compounds(parent,'cid')
        psmiles = parent_molecule[0].canonical_smiles
        pmol = Chem.MolFromSmiles(psmiles)
        prdkit_smiles = Chem.MolToSmiles(pmol)
        prdkit_inchi = Chem.MolToInchi(pmol)
        prdkit_inchikey = Chem.MolToInchiKey(pmol)
        prdkit_molformula = rdMolDescriptors.CalcMolFormula(pmol)
        prdkit_molweight = round(Descriptors.MolWt(pmol), 2)
        parent_list = [prdkit_inchikey,prdkit_smiles,prdkit_inchi,prdkit_molformula,prdkit_molweight]
        pypubdata['associated_compounds'] = {}
        pypubdata['associated_compounds']['parent_molecule'] = parent_list
    logPs, references = exp_logp_pubchem(outcoming_xml_data)
    if len(logPs) != 0:
        pypubdata['experimental_data'] = {}
        pypubdata['experimental_data']['partition_coefficient'] = {}
        for element in range(0, len(logPs)):
            data_id_string = 'pubchemvalue' + re.sub("[^0-9]", "", logPs[element])
            pypubdata['experimental_data']['partition_coefficient'][data_id_string] = {}
            pypubdata['experimental_data']['partition_coefficient'][data_id_string]['logp'] = logPs[element]
            pypubdata['experimental_data']['partition_coefficient'][data_id_string]['type_id'] = 'exp_logp_lit'
            pypubdata['experimental_data']['partition_coefficient'][data_id_string]['measurement_metadata'] = references[element]
    return pypubdata

























