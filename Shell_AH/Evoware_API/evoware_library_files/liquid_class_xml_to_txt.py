# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 13:05:24 2020

This script is very infrequently needed

The liquid class files are defined in XML, but we need to have them ready to
process in the schedule generation, so here we are getting the names of them,
the details are not import to what we are doing... so basically we just need to
get the names of everything (this will be used for liquid class error checking)

Anytime there are new liquid classes added this needs to be updated

Towards the bottom there is a "signature" line that results in a parsing error,
just delete the line and this will work just fine!

In the customLC xml file there is also an Evo Child towards the top that is not
needed and yields an error that should be avoided

@author: Brent Koscher
"""

import xml.etree.ElementTree as ET

filepath = r'.\xml_files\DefaultLCs.XML'
outpath = r'.\DefaultLCs.txt'

xml_tree = ET.parse(filepath)
root = xml_tree.getroot()

with open(outpath, 'w') as outfile:
    for child in root:
        outfile.write(child.get('name') + '\t' + child.get('liquidName') + '\n')
        outfile.flush()

filepath = r'.\xml_files\CustomLCs.XML'
outpath = r'.\CustomLCs.txt'

xml_tree = ET.parse(filepath)
root = xml_tree.getroot()

with open(outpath, 'w') as outfile:
    for child in root:
        #print(child.get('name'), child.get('liquidName'))
        outfile.write(child.get('name') + '\t' + child.get('liquidName') + '\n')
        outfile.flush()