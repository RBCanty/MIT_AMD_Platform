# -*- coding: utf-8 -*-
"""
Created on 09/27/2022

Implementation of code the build queue documents, this version will be included in the
prediction pipeline for the AMD platform

@author: Brent Koscher
"""

import sys
sys.path.append(r'.\queue_generation')
sys.path.append(r'.\reaction_grouping\queue_generation')

import glob
import datetime
import os
from queue_workflows import queue_workflows
import yaml


def generate_queue_document(output_directory, grouping_document, grouping_information):
    current_time = datetime.datetime.now()
    # Check to make sure that a campaign name was provided and make the output filename
    if 'campaign_name' not in grouping_document.keys():
        return ['Error', 'Campaign name not provided in the grouping document']
    output_filepath = r'%s\%s_generated_queue.yaml' % (output_directory, grouping_document['campaign_name'])

    # Look for the workflow type that needs to be executed
    workflow_type = grouping_document['workflow_type']
    return_statement, queue_documents, reagents_needed = queue_workflows(workflow_type, grouping_document,
                                                                         grouping_information)
    if return_statement != 'Success':
        return ['Error', queue_documents]
    print(output_filepath)
    with open(output_filepath, 'w') as yaml_file:
        yaml.dump(queue_documents, yaml_file)
    reagents_output = output_filepath.replace('.yaml', '_reagents.yaml')
    with open(reagents_output, 'w') as yaml_file:
        yaml.dump(reagents_needed, yaml_file)
    return ['Success', 'Generated queue documents and reagents needed for %s' % grouping_document['campaign_name']]


if __name__ == "__main__":
    filepath = r'C:\Users\Brent\Desktop\Data Storage\Data Files\20230517\photochem_discovery_20230529.yaml'
    with open(filepath, 'r') as yaml_file:
        grouping_details = yaml.load(yaml_file, Loader=yaml.Loader)

    output_directory = r'.\queue_documents\\'
    returned_statement = generate_queue_document(output_directory, grouping_details, {})
    print(returned_statement)
