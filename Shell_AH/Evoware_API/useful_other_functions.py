# -*- coding: utf-8 -*-
"""
These are supporting functions that end up being useful for the Tecan functions

@author: brent
"""
import copy
import math
import os
import re
import sys
import natsort
from pprint import pprint
from copy import deepcopy
import traceback
import itertools
import numpy as np
import random
import yaml

"""
Not-cleaned functions
"""


def check_for_problems(current_bed, location, labware_type):
    small_labware = ['96 Well Microplate']
    large_labware = ['96 Well DeepWell', 'Paradox Thermal Plate']
    xlarge_labware = ['8 Position Vial Carrier', 'Inert Box Block', 'DiTi SBS Waste', '12 Position Vial Carrier']

    xlarge_locations_to_check = []
    large_locations_to_check = []
    small_locations_to_check = []
    problems = 0
    if location == [32, 1, 96]:
        xlarge_locations_to_check = [[38, 1, 96], [38, 2, 96]]
        large_locations_to_check = [[38, 1, 96], [38, 2, 96]]
        small_locations_to_check = []
    elif location == [32, 2, 96]:
        xlarge_locations_to_check = [[38, 1, 96], [38, 2, 96], [38, 3, 96]]
        large_locations_to_check = [[38, 1, 96], [38, 2, 96], [38, 3, 96]]
        small_locations_to_check = []
    elif location == [38, 1, 96]:
        xlarge_locations_to_check = [[44, 1, 336], [44, 2, 336]]
        large_locations_to_check = [[44, 1, 336], [44, 2, 336]]
        small_locations_to_check = []
    elif location == [38, 2, 96]:
        xlarge_locations_to_check = [[44, 1, 336], [44, 2, 336]]
        large_locations_to_check = [[44, 1, 336], [44, 2, 336]]
        small_locations_to_check = []
    elif location == [38, 3, 96]:
        xlarge_locations_to_check = [[44, 2, 336]]
        large_locations_to_check = [[44, 2, 336]]
        small_locations_to_check = []
    elif location == [44, 1, 336]:
        xlarge_locations_to_check = []
        large_locations_to_check = []
        small_locations_to_check = []
    elif location == [44, 2, 336]:
        xlarge_locations_to_check = []
        large_locations_to_check = []
        small_locations_to_check = []
    elif location == [25, 1, 323]:
        xlarge_locations_to_check = [[32, 1, 96], [32, 2, 96]]
        large_locations_to_check = []
        small_locations_to_check = []
    elif location == [25, 2, 323]:
        xlarge_locations_to_check = [[32, 2, 96], [32, 3, 96]]
        large_locations_to_check = []
        small_locations_to_check = []
    elif location == [35, 9, 84]:
        xlarge_locations_to_check = [[32, 1, 96]]
        large_locations_to_check = [[32, 1, 96]]
        small_locations_to_check = [[32, 1, 96]]
    elif location == [35, 8, 84]:
        xlarge_locations_to_check = [[32, 1, 96]]
        large_locations_to_check = []
        small_locations_to_check = []
    elif location == [39, 1, 84]:
        xlarge_locations_to_check = [[38, 1, 96]]
        large_locations_to_check = [[38, 1, 96]]
        small_locations_to_check = [[38, 1, 96]]
    elif location == [39, 2, 84]:
        xlarge_locations_to_check = [[38, 1, 96]]
        large_locations_to_check = []
        small_locations_to_check = []
    elif location == [18, 5, 333]:
        xlarge_locations_to_check = [[25, 2, 323], [25, 4, 323]]
        large_locations_to_check = [[25, 2, 323], [25, 4, 323]]
    elif location == [18, 2, 333]:
        xlarge_locations_to_check = [[25, 1, 323]]
        large_locations_to_check = [[25, 1, 323]]

    for xlarge_location in xlarge_locations_to_check:
        if current_bed[xlarge_location[0]][xlarge_location[2]]['labware_types'][
            xlarge_location[1] - 1] in xlarge_labware:
            problems += 1
            # print('Found %s in %s' % (current_bed[xlarge_location[0]][xlarge_location[2]]['labware_types'][xlarge_location[1] - 1], xlarge_location))
    for large_location in large_locations_to_check:
        if current_bed[large_location[0]][large_location[2]]['labware_types'][large_location[1] - 1] in large_labware:
            problems += 1
            # print('Found %s in %s' % (current_bed[large_location[0]][large_location[2]]['labware_types'][large_location[1] - 1], large_location))
    for small_location in small_locations_to_check:
        if current_bed[small_location[0]][small_location[2]]['labware_types'][small_location[1] - 1] in small_labware:
            problems += 1
            # print('Found %s in %s' % (current_bed[small_location[0]][small_location[2]]['labware_types'][small_location[1] - 1], small_location))

    # Here are the locations which cannot be accessed by the RoMa

    if location == [44, 3, 336]:  # This is the inert box location on the bed
        problems += 1
    elif location == [22, 7, 84]:  # Not accessible due to the TeVaC
        problems += 1
    elif location == [22, 8, 84]:  # Not accessible due to the TeVaC
        problems += 1
    elif location == [22, 9, 84]:  # Not accessible due to the TeVaC
        problems += 1
    elif location == [28, 2, 329]:  # Not accessible due to a heater-shaker in the way
        problems += 1
    elif location == [65, 6, 84]: # Not accessible due to the heater-shaker
        problems += 1
    elif location == [65, 7, 84]: # Not accessible due to the heater-shaker
        problems += 1
    elif location == [65, 8, 84]: # Not accessible due to the heater-shaker
        problems += 1
    elif location == [65, 9, 84]: # Not accessible due to the heater-shaker
        problems += 1
    if problems > 0:
        return 'Unacceptable'
    else:
        return 'Acceptable'


# ---------------------------------------------------------------------------------------------------------

"""
Cleaned functions
"""


def find_value(dictionary, relevant_keys, relevant_values, current_path=None, dest=None):
    """
    Parameters
    ----------
    dictionary : dict
        The current platform resource library nested dictionary
    relevant_keys : list of strings
        List with the important key to find in the nested dictionary
    relevant_values : list of strings
        List with the important values to find in the nested dictionary

    Returns
    -------
    list
        List with all instances found in the nested library
    """
    if current_path is None:
        current_path = []
    if dest is None:
        dest = []
    for key, value in dictionary.items():
        growing_dest = current_path + [key]
        if isinstance(value, dict):
            find_value(value, relevant_keys, relevant_values, growing_dest, dest)
        elif isinstance(value, list):
            if key in relevant_keys:
                flat_value_list = list(itertools.chain(*value))
                if key in relevant_keys and len(set(flat_value_list).intersection(set(relevant_values))) > 0:
                    try:
                        dest.append(growing_dest)
                    except KeyError:
                        dest.append(growing_dest)
        else:
            if key in relevant_keys and value in relevant_values:
                try:
                    dest.append(growing_dest)
                except KeyError:
                    dest.append(growing_dest)
    return dest


def solvent_finder_origin(library, solvents, solvent_dict):
    """
    Parameters
    ----------
    library : dict
        The current platform resource library dictionary
    solvents : list of lists
        List with solvents to find in the library
    solvent_dict : dict
        Endpoint for future implementation

    Returns
    -------
    list
        Status code : String
            Either 'Error' or 'Success'
        Function return : Dict or String
            Either dictionary of solvent details or exception as a string
    """
    try:
        special_solvent_labware = ['8 Position Vial Carrier']
        solvents_return = {}
        relevant_collections = ['reagents', 'solvents']
        for collection in relevant_collections:
            for container in library[collection].keys():
                pass
        solvent_problems = []
        for solvent in solvents:
            return_object = find_value(library, ['chemical_name', 'chemical_smiles'], [solvent[0]])
            if not return_object:
                solvent_problems.append('Problem: %s not found' % solvent)
                continue
            max_solvent_volume = [0, []]
            for solvent_return in return_object:
                if solvent_return[0] not in relevant_collections:
                    continue
                solvent_volume = library[solvent_return[0]][solvent_return[1]][solvent_return[2]][solvent_return[3]][
                    'volume_ul']
                if solvent_volume > max_solvent_volume[0]:
                    max_solvent_volume = [solvent_volume, solvent_return]
            if not max_solvent_volume[1]:
                solvent_problems.append('Problem: %s not found in suitable location' % solvent[0])
                continue
            best_location = library[max_solvent_volume[1][0]][max_solvent_volume[1][1]][max_solvent_volume[1][2]][
                max_solvent_volume[1][3]]
            solvents_return[solvent[0]] = {'chemical_name': best_location['chemical_name'],
                                           'chemical_smiles': best_location['chemical_smiles'],
                                           'container_id': max_solvent_volume[1][1],
                                           'collection': max_solvent_volume[1][0],
                                           'origin_well_location': max_solvent_volume[1][3],
                                           'container_name':
                                               library[max_solvent_volume[1][0]][max_solvent_volume[1][1]][
                                                   'container_name'],
                                           'solvent_volume': best_location['volume_ul'],
                                           'volume_needed': 0}
    except:
        return ['Error', traceback.format_exc()]
    if solvent_problems:
        return ['Error', solvent_problems]
    else:
        return ['Success', solvents_return]


def solvent_rinse_location(library, solvent_name):
    """
    Parameters
    ----------
    library : dict
        The current platform resource library dictionary
    solvent_name : name
        Dictionary with solvents to find in the library

    Returns
    -------
    list
        Status code : String
            Either 'Error' or 'Success'
        Function return : Dict or String
            Either dictionary of solvent details or exception as a string
    """
    relevant_collections = ['solvents', 'reagents']
    return_object = find_value(library, ['chemical_name', 'chemical_smiles'], [solvent_name])
    if not return_object:
        return ['Error', 'Problem: %s not found' % solvent_name]
    max_solvent_volume = [0, []]
    for solvent_return in return_object:
        if solvent_return[0] not in relevant_collections:
            continue
        solvent_volume = library[solvent_return[0]][solvent_return[1]][solvent_return[2]][solvent_return[3]][
            'volume_ul']
        if solvent_volume > max_solvent_volume[0] and solvent_volume > 0.9 * 5000:
            max_solvent_volume = [solvent_volume, solvent_return]
    if not max_solvent_volume[1]:
        return ['Error', 'Problem: %s not found in suitable location' % solvent_name]
    return_location = library[max_solvent_volume[1][0]][max_solvent_volume[1][1]]['location'][1]
    return ['Success', return_location]


def reagent_finder(library, reagents, reagent_dict):
    """
    Parameters
    ----------
    library : dict
        The current platform resource library dictionary
    reagents : dict
        Dictionary with reagents to find in the library
    reagent_dict : dict
        Endpoint for future implementation

    Returns
    -------
    list
        Status code : String
            Either 'Error' or 'Success'
        Function return : Dict or String
            Either dictionary of reagent details or exception as a string
    """
    try:
        reagents_return = {}
        reagent_problems = []
        relevant_collections = ['reagents']
        if reagents is None:
            return ['Success', reagents_return]
        for reagent in reagents:
            return_object = find_value(library, ['chemical_name', 'chemical_smiles'], [reagent[0]])
            if return_object == []:
                reagent_problems.append('Problem: %s not found' % reagent[0])
                continue
            max_moles_available = [0, []]
            for reagent_return in return_object:
                if reagent_return[0] not in relevant_collections:
                    continue
                reagent_definition = deepcopy(
                    library[reagent_return[0]][reagent_return[1]][reagent_return[2]][reagent_return[3]])
                reagent_moles = reagent_definition['concentration_molar'] * reagent_definition['volume_ul']
                if reagent_moles > max_moles_available[0]:
                    max_moles_available = [reagent_moles, reagent_return]
            if max_moles_available[1] == []:
                reagent_problems.append('Problem: %s not found in suitable container' % reagent[0])
                continue
            best_location = library[max_moles_available[1][0]][max_moles_available[1][1]][max_moles_available[1][2]][
                max_moles_available[1][3]]
            if len(reagent) == 2:
                chemical_sequence = 1
            else:
                chemical_sequence = reagent[2]
            if reagent[0] in reagents_return.keys():
                reagents_return[reagent[0]]['amount_needed'] += float(reagent[1])
                reagents_return[reagent[0]]['volume_needed'] += float(reagent[1]) / best_location[
                    'concentration_molar'] * 1E6
            else:
                reagents_return[reagent[0]] = deepcopy(best_location)
                reagents_return[reagent[0]].update({'amount_needed': float(reagent[1]),
                                                    'volume_needed': float(reagent[1]) / best_location[
                                                        'concentration_molar'] * 1E6,
                                                    'transfer_sequence': chemical_sequence,
                                                    'container_id': max_moles_available[1][1],
                                                    'collection': max_moles_available[1][0],
                                                    'origin_well_location': best_location['plate_well']})
        #if len(reagent_problems) != 0:
        #    return ['Error', reagent_problems]
    except:
        return ['Error', traceback.format_exc()]
    return ['Success', reagents_return]


def previous_product_finder(library, target_product_tree, reactants):
    """
    Parameters
    ----------
    library : dict
        The current platform resource library dictionary
    target_product_tree : string
        String that specifies which reaction tree to find a previous product from
    reactants : list
        List that has all reactants to find in the database, only the ones that
        match previous products and not reagents are identified

    Returns
    -------
    list
        Status code : String
            Either 'Error' or 'Success'
        Function return : Dict or String
            Either dictionary of previous product details or exception as a string
    """
    try:
        products_return = {}
        relevant_collections = ['wellplates']
        product_problems = []
        if reactants is None:
            return ['Success', products_return]
        for reactant in reactants:
            return_object = find_value(library, ['target_product'], [reactant[0]])
            if not return_object:
                continue
            #if len(return_object) != 1:
            #    product_problems.append('Problem: Too many %s instances (%s found)' % (reactant[0], len(return_object)))
            #    continue
            best_final_product = [0, []]
            for product_return in return_object:
                if product_return[0] not in relevant_collections:
                    continue
                product_definition = deepcopy(
                    library[product_return[0]][product_return[1]][product_return[2]][product_return[3]])
                if 'final_product' not in product_definition.keys() or product_definition['final_product'][
                    2] != target_product_tree:
                    continue
                if len(reactant) == 2:
                    container_category_target = 'filtrate_plate'
                else:
                    container_category_target = reactant[2]
                if 'container_category' not in library[product_return[0]][product_return[1]].keys() or \
                        library[product_return[0]][product_return[1]][
                            'container_category'] != container_category_target:
                    continue
                # if product_definition['confirmation'] in ['none', 'Testing'] continue
                # how to decide which to use if multiple?
                best_final_product = [0, product_return]
            if not best_final_product[1]:
                product_problems.append('Problem: No suitable instances of %s remain' % reactant[0])
                continue
            previous_reactant_string = '%s>%s' % (reactant[0], target_product_tree)
            if previous_reactant_string not in products_return.keys():
                if 'product_remaining' in product_definition.keys() and type(product_definition['product_remaining']) == list:
                    previous_mols = product_definition['product_remaining'][0][0][1]
                else:
                    previous_mols = product_definition['target_product'][0][1]
                products_return[previous_reactant_string] = {'target_product': reactant[0],
                                                'collection': best_final_product[1][0],
                                                'container_id': best_final_product[1][1],
                                                'origin_well_location': best_final_product[1][3],
                                                'container_name': library[product_return[0]][product_return[1]][
                                                    'container_name'],
                                                'relevant_field': best_final_product[1][4],
                                                'previous_mols':previous_mols,
                                                'amount_needed': reactant[1],
                                                'previous_solvent': product_definition['solvents'],
                                                'previous_volume': product_definition['total_volume'],
                                                'transfer_sequence': 1,
                                                'reaction_scale': reactant[1] / previous_mols}
    except:
        return ['Error', traceback.format_exc()]
    if product_problems:
        return ['Error', product_problems]
    return ['Success', products_return]


def find_labware(current_bed, library, labware_name, labware_type, target_container='', exceptions=[]):
    """
    Parameters
    ----------
    current_bed : dict
        The configuration and labware on the liquid handler at call in method,
        this checks to see if there is an access issue and ignores wellplates
    library : dict
        The current platform resource library
    labware_name : string or None
        Used to find a labware when specified as string
    labware_type : string
        Used to find an available labware when not equal to 'name_lookup_only'
    target_container : string
        Hold-over from previous version, may be deleted in future code cleanups
    exceptions : list
        Contains all exceptions of wellplates to avoid

    Returns
    -------
    list
        Status code : String
            Either 'Error' or 'Success'
        Function return : Dict or String
            Either dictionary of a suitable labware for the script to use
    """
    lpx_flag = 0
    lowest_index = [[1E6, 1E6, 1E6], []]
    if labware_name is None and labware_type != 'name_only_lookup':
        return_object = find_value(library, ['labware_type'], [labware_type])
    elif labware_name is not None:
        return_object = find_value(library, ['container_name'], [labware_name])
        if len(return_object) != 1:
            return ['Error', 'Bad number of documents (%s) with container name %s' % (len(return_object), labware_name)]
    else:
        return ['Error', 'No search method for wellplate_name "%" and labware_type "%s"' % (labware_name, labware_type)]

    for labware_keys in return_object:
        labware_details = library[labware_keys[0]][labware_keys[1]]
        if labware_type in ['96 Well Filtration Plate', 'Inert Box Block']:
            if labware_keys[0] in ['consumables'] and labware_type == 'Inert Box Block':
                pass
            elif labware_keys[0] not in ['consumables'] or labware_details['consumable_status'] != 'ready':
                continue
        else:
            if labware_keys[0] not in ['wellplates', 'reagents']:
                continue
            if labware_name is None:
                if labware_details['contents'] != 'Empty':
                    continue
            if labware_details['location'][0] == 'lpx':
                lpx_flag = 1
            if labware_details['location'][0] != 'liquid_handler' or labware_details['container_name'] in exceptions:
                continue
            problem_return = check_for_problems(current_bed, labware_details['location'][1], labware_type)
            if problem_return != 'Acceptable':
                continue
        
        if all(lowest_index[0][idx] >= labware_details['location'][1][idx] for idx in range(0, len(lowest_index[0]))):
            lowest_index = [labware_details['location'][1], labware_keys]
    if lowest_index[1] == []:
        return ['Error', 'No suitable labware on the liquid handler: >%s>%s>%s' % (lpx_flag, labware_type, target_container)]
    labware_return = deepcopy(library[lowest_index[1][0]][lowest_index[1][1]])
    labware_return['collection'] = lowest_index[1][0]
    return ['Success', labware_return]


def find_accessible_open_location(current_bed, exceptions=[], labware_type='96 Well Microplate'):
    """
    Parameters
    ----------
    current_bed : dict
        The configuration and labware on the liquid handler at call in method,
        this checks to see if there is an access issue and ignores wellplates
    exceptions : list
        Contains all exceptions of wellplates to avoid
    labware_type : string
        Default: 96 Well Microplate
        Used by check for problems to make sure that transfers can happen

    Returns
    -------
    list
        Status code : String
            Either 'Error' or 'Success'
        Function return : List or String
            Either Bed location list or error message string
    """
    acceptable_carrier_ids = [96, 336]
    acceptable_positions = []
    for grid in current_bed.keys():
        for carrier in current_bed[grid].keys():
            if carrier not in acceptable_carrier_ids:
                continue
            for site_index, site in enumerate(current_bed[grid][carrier]['labware_types']):
                if site != '':
                    continue
                bed_location = [grid, site_index + 1, carrier]
                if bed_location not in exceptions:
                    problem_return = check_for_problems(current_bed, bed_location, labware_type)
                    if problem_return == 'Acceptable':
                        return ['Success', bed_location]
    return ['Error', 'No acceptable bed positions exist for "%s"' % labware_type]


def liha_grouping(destinations, carriers_labware, origin_labware_type, destination_labware_type):
    """
    Parameters
    ----------
    destinations : list of lists
        The transfer destination for LiHa to make [[Well, Volume, sub_group], ...]
    carriers_labware : dict
        Dictionary of the carriers and labware in the Evoware library
    origin_labware_type : string
        Used to determine what the layout of the wells of the origin labware
    destination_labware_type : string
        Used to determine what the layout of the wells of the destination labware

    Returns
    -------
    list
        Status code : String
            Either 'Error' or 'Success'
        Function return : Dict or String
            Either grouped destinations dict or error message string
    """
    origin_labware = carriers_labware['labware'][origin_labware_type]
    destination_labware = carriers_labware['labware'][destination_labware_type]
    try:
        nrows = destination_labware['nrows']
        ncols = destination_labware['ncols']
    except:
        return ['Error', traceback.format_exc()]
    well_groups = {}
    for col in range(0, ncols):
        well_groups['Col%s' % str(col)] = []
        for row in range(0, nrows):
            well_groups['Col%s' % str(col)].append(chr(ord('A') + row) + str(col + 1))
    grouped_destinations = {}
    for destination in destinations:
        if destination[1] <= 0:
            continue
        for group in well_groups.keys():
            if destination[0] in well_groups[group]:
                if group not in grouped_destinations.keys():
                    grouped_destinations[group] = []
                match = re.match(r"([a-z]+)([0-9]+)", destination[0], re.I)
                wellname_split = list(match.groups())
                grouped_destinations[group].append(destination + [ord(wellname_split[0]) - ord('A') + 1])
    group_return = {}
    for group_dest in grouped_destinations.keys():
        if any(item[1] >= 320 for item in grouped_destinations[group_dest]):
            divisions = [math.ceil(item[1] / 320) for item in grouped_destinations[group_dest]]
            for i in range(1, max(divisions) + 1):
                sub_group = []
                for item_idx, item in enumerate(grouped_destinations[group_dest]):
                    if divisions[item_idx] >= i:
                        sub_group.append([item[0], item[1] / divisions[item_idx], item[-1]])
                group_return[group_dest + '_' + str(i)] = sub_group
        else:
            sub_group = []
            for item_idx, item in enumerate(grouped_destinations[group_dest]):
                sub_group.append([item[0], item[1], item[-1]])
            group_return[group_dest + '_' + str(item_idx + 1)] = sub_group
    return ['Success', group_return]


def priority_sort(incoming_dict, sorting_type='Normal'):
    """
    Parameters
    ----------
    incoming_dict : dict
        Dictionary of the chemicals to transfer
    sorting_type : string
        String to specifiy the type of sorting

    Returns
    -------
    list
        Status code : String
            Either 'Error' or 'Success'
        Function return : Dict or String
            Either chemical sequence dict or error message string
    """

    """
    ##### OTHER SORITNG METHOD FROM PREVIOUS VERSIONS #####
    def priority_sort_maintain_order(incoming_list):
        holder_sort = incoming_list
        holder_dict = {}
        for item in holder_sort:
            if item[1] not in holder_dict.keys():
                holder_dict[item[1]] = [item]
            else:
                holder_dict[item[1]].append(item)
        holder_list = [[key, holder_dict[key]] for key in holder_dict.keys()]
        return ['Success', holder_list]
    """

    try:
        chemical_sequence = {}
        for chemical in incoming_dict.keys():
            destinations = incoming_dict[chemical]['destinations']
            for transfer_sequence in destinations.keys():
                if transfer_sequence not in chemical_sequence.keys():
                    chemical_sequence[transfer_sequence] = {}
                container_id = incoming_dict[chemical]['container_id']
                if container_id not in chemical_sequence[transfer_sequence].keys():
                    chemical_sequence[transfer_sequence][container_id] = {}
                if chemical not in chemical_sequence[transfer_sequence][container_id].keys():
                    chemical_sequence[transfer_sequence][container_id][chemical] = {'destinations': [],
                                                                                    'container_name':
                                                                                        incoming_dict[chemical][
                                                                                            'container_name']}
                for target_destination in incoming_dict[chemical]['destinations'][transfer_sequence]:
                    chemical_sequence[transfer_sequence][container_id][chemical]['destinations'].append(
                        target_destination)
        return ['Success', chemical_sequence]
    except Exception:
        return ['Error', traceback.format_exc()]


def find_open_heater_shaker(current_bed, temperature, plate_type, platform_locations, location_type='standard',
                            exceptions=[]):
    """
        Parameters
        ----------
        current_bed : dict
            Dictionary of the chemicals to transfer
        temperature : float
            Target temperature for the platform to find
        plate_type : string
            Specifies the wellplate type
        platform_locations : dict
            Contains locations for the heater shakers present on the platform
        location_type : string
            Specifies a special string of the specific type of site, default is standard
        exceptions : list
            Contains list exceptions to sites to ignore

        Returns
        -------
        list
            Status code : String
                Either 'Error' or 'Success'
            Function return : List or String
                Either acceptable location string or error message string
    """
    if location_type == 'low_temperature_prep':
        acceptable_pos = platform_locations['low_temperature_prep_locations']
    elif location_type == 'standard' and temperature in range(5, 79):
        acceptable_pos = platform_locations['thermoshake_locations']
    elif location_type == 'standard' and temperature in range(25, 120):
        acceptable_pos = platform_locations['teleshake_locations']
    elif temperature not in range(4, 120):
        return ['Error', 'Problem: Temperature (%s) outside of the valid ranges for heater-shakers' % str(temperature)]
    else:
        return ['Error', 'Problem: A location of %s with temperature %s not implemented' % (location_type, temperature)]
    for bed_pos in acceptable_pos:
        if bed_pos in exceptions:
            continue
        carrier_location = current_bed[bed_pos[0]][bed_pos[2]]['labware_types'][bed_pos[1] - 1]
        if carrier_location == '':
            problem_return = check_for_problems(current_bed, bed_pos, plate_type)
            if problem_return != 'Acceptable':
                continue
            else:
                acceptable_location = bed_pos
                return ['Success', acceptable_location]
    else:
        return ['Error', 'Problem: No valid locations open for Temperature = %s' % str(temperature)]


def find_accessible_hotel_location(current_bed, hotel_type, plate_type, platform_locations, exceptions=[]):
    """
        Parameters
        ----------
        current_bed : dict
            Dictionary of the chemicals to transfer
        hotel_type : string
            Specifies the type of hotel (storage, transfer, etc.)
        plate_type : string
            Specifies the wellplate type
        platform_locations : dict
            Contains locations for the heater shakers present on the platform
        exceptions : list
            Contains list exceptions to sites to ignore

        Returns
        -------
        list
            Status code : String
                Either 'Error' or 'Success'
            Function return : List or String
                Either acceptable location string or error message string
    """
    acceptable_destinations = ['transfer', 'storage']
    if hotel_type not in acceptable_destinations:
        return ['Error', 'Problem: Destination type %s not implemented for hotel transfers' % hotel_type]

    storage_hotel_labware = {'96 Well Microplate': ['Hotel 9Pos Microplate'],
                             '96 Well Filtration Plate': ['Hotel 3Pos Filterplate'],
                             '96 Well DeepWell': ['Hotel 3Pos DeepWell'],
                             'Paradox Thermal Plate': ['Hotel 2Pos Paradox'], }

    transfer_hotel_labware = {'standard_transfer': {'plate_types': ['Paradox Thermal Plate', '96 Well DeepWell',
                                                                    '96 Well Microplate Half Area',
                                                                    '96 Well PCR Plate'],
                                                    'locations': platform_locations['standard_transfer_locations']},
                              'special_transfer': {'plate_types': ['96 Well Microplate'],
                                                   'locations': platform_locations['special_transfer_locations']}}

    # We need to find acceptable positions that we could transfer the wellplate
    # into that are available on the liquid handler platform
    if hotel_type == 'transfer':
        if plate_type in transfer_hotel_labware['standard_transfer']['plate_types']:
            acceptable_bed_pos = transfer_hotel_labware['standard_transfer']['locations']
        elif plate_type in transfer_hotel_labware['special_transfer']['plate_types']:
            acceptable_bed_pos = transfer_hotel_labware['standard_transfer']['locations']
        else:
            return ['Error', 'Problem: Plate type %s is not pre-specified' % plate_type]
    elif hotel_type == 'storage':
        acceptable_bed_pos = []
        if plate_type in storage_hotel_labware.keys():
            acceptable_hotels = storage_hotel_labware[plate_type]
            for bed_grid in current_bed.keys():
                for carrier_id in current_bed[bed_grid].keys():
                    if current_bed[bed_grid][carrier_id]['carrier_name'] in acceptable_hotels:
                        acceptable_bed_pos.extend([[bed_grid, site_index + 1, carrier_id] for site_index, element
                                                   in enumerate(current_bed[bed_grid][carrier_id]['labware_labels'])
                                                   if element == ''])
            if len(acceptable_bed_pos) == 0:
                return ['Error', 'Problem: No acceptable sites exist for %s' % plate_type]
        else:
            return ['Error', 'Problem: Plate type %s is not pre-specified' % plate_type]
    else:
        return ['Error', 'Problem: Destination type not fully implemented for hotel transfers']

    # With the potential transfer locations found, we can then check to find the
    # first one that does not have a transfer conflict
    for bed_pos in acceptable_bed_pos:
        if bed_pos in exceptions:
            continue
        problem_return = check_for_problems(current_bed, bed_pos, plate_type)
        if problem_return == 'Acceptable':
            return ['Success', bed_pos]
    return ['Error', 'Problem: No acceptable location exists for %s, %s' % (plate_type, hotel_type)]


def check_spark_availability(current_bed):
    """
        Parameters
        ----------
        current_bed : dict
            Dictionary of the chemicals to transfer

        Returns
        -------
        list
            Status code : String
                Either 'Error' or 'Success'
            Function return : List or String
                Either acceptable location string or error message string
    """
    spark_carrier_id = 334
    for grid in current_bed.keys():
        for carrier in current_bed[grid].keys():
            if carrier == spark_carrier_id:
                position = current_bed[grid][carrier]['labware_types'][0]
                if position == '':
                    return ['Success', [grid, 1, spark_carrier_id]]
    return ['Error', 'Spark is currently occupied by another wellplate']


def evaporation_estimation(solvents, labware_type, elapsed_time, temperature, model='epa'):
    """
        Parameters
        ----------
        solvents : list
            Contains all the solvents and their amounts
        labware_type : str
            Exists to change the evaporation rate for different well-plate types
        elapsed_time : number
            Seconds since the last volume reading for the incoming well
        temperature : number
            Temperature experienced since last well update
        model : str
            Allows for switching between different evaporation models (currently EPA is the only one)

        Returns
        -------
        list
            Status code : str
                Either 'Error' or 'Success'
            Function return : list or str
                Either updated solvents field or error message from the function
    """
    solvent_details = {'dmf': {'mw': 73.09, 'dhvap': 43.6, 'density': 944, 'vp': 3.87},
                       'dmso': {'mw': 78.13, 'dhvap': 52.9, 'density': 1101, 'vp': 0.593},
                       'isopropanol': {'mw': 60.1, 'dhvap': 39.85, 'density': 786, 'vp': 32.93},
                       'ethanol': {'mw': 46.07, 'dhvap': 38.56, 'density': 789, 'vp': 44.629},
                       'methanol': {'mw': 32.04, 'dhvap': 35.21, 'density': 792, 'vp': 97.658},
                       'chloroform': {'mw': 119.4, 'dhvap': 29.4, 'density': 1490, 'vp': 196.79},
                       'dioxane': {'mw': 88.11, 'dhvap': 34.16, 'density': 1030, 'vp': 38.103},
                       'water': {'mw': 43.9, 'dhvap': 43.09, 'density': 997, 'vp': 23.702},
                       'acetonitrile': {'mw': 29.8, 'dhvap': 29.8, 'density': 786, 'vp': 32.252},
                       'acetone': {'mw': 31.27, 'dhvap': 31.27, 'density': 784, 'vp': 230.269},
                       'ethyl_acetate': {'mw': 88.11, 'dhvap': 31.94, 'density': 902, 'vp': 93.2}}

    updated_solvents = []
    for solvent in solvents:
        if solvent[0] == 'water':
            updated_solvents.append(solvent)
        else:
            if solvent[0] in solvent_details.keys():
                relevant_solvent = solvent_details[solvent[0]]
            else:
                relevant_solvent = {'mw': 80, 'dhvap': 35, 'density': 800, 'vp': 100}
            if model == 'epa':
                vapor_pressure = relevant_solvent['vp'] * np.exp((relevant_solvent['dhvap'] * 1000 / 8.3145) *
                                                                 ((1 / (273 + 25)) - (1 / (273 + temperature))))
                r_value = 8.205736E-5 * 100**3 * 760
                wind_speed = 0.001
                gas_phase_transfer_coef = 0.25 * wind_speed**0.78 * (18 / relevant_solvent['mw'])**(1/3)
                surface_area = 0.5
                rate = (1000 / relevant_solvent['density']) * 1E3 * relevant_solvent['mw'] * gas_phase_transfer_coef * \
                    surface_area * vapor_pressure / (r_value * (273 + temperature))

                # Now check to see how much volume might be left in the wells
                volume_remaining = solvent[1] - rate * elapsed_time
                if volume_remaining < 10:
                    continue
                updated_solvents.append([solvent[0], volume_remaining, solvent[2]])

    return ['Success', updated_solvents]


def chemical_characterization_lookup(library, product_name, product_dict):
    """
        Parameters
        ----------
        library : dict
            Contains reagents that need to be found
        product_name : str
            The chemical name (for reagents) or chemical smiles (for both reagents and products)
        product_dict : dict
            Details about the product from the queue, gives basic guiding information

        Returns
        -------
        list
            Status code : str
                Either 'Error' or 'Success'
            Function return : list or str
                Either pointer list of keys or error message from the function
    """
    try:
        if product_dict['container_category'] == 'reagent_tray':
            relevant_collections = ['reagents']
            relevant_fields = ['chemical_smiles', 'chemical_name']
        else:
            relevant_collections = ['wellplates']
            relevant_fields = ['target_product']
        return_object = find_value(library, relevant_fields, [product_name])
        best_return = [0, []]
        for reagent_return in return_object:
            if reagent_return[0] not in relevant_collections:
                continue
            container_details = library[reagent_return[0]][reagent_return[1]]
            if container_details['location'][0] != 'liquid_handler':
                continue
            if product_dict['container_category'] == 'reagent_tray':
                reagent_concentration = container_details[reagent_return[2]][reagent_return[3]]['concentration_molar']
                reagent_volume = container_details[reagent_return[2]][reagent_return[3]]['volume_ul']
                material_amount = reagent_volume * reagent_concentration
                if material_amount > best_return[0]:
                    best_return = [material_amount, reagent_return]
            else:
                if 'container_category' not in container_details.keys() or container_details['container_category'] != \
                                                                                    product_dict['container_category']:
                    continue
                target_products = container_details[reagent_return[2]][reagent_return[3]][reagent_return[4]]
                for target_product in target_products:
                    if target_product[0] != product_name:
                        continue
                    material_amount = target_product[1]
                    if material_amount > best_return[0]:
                        best_return = [material_amount, reagent_return]
        if not best_return[1]:
            return ['Missing Chemical', 'No product found']
        return ['Success', best_return[1]]
    except:
        return ['Error', traceback.format_exc()]


def multi_transfer_grouping(target_destinations, max_single_transfer_volume):
    """
        Parameters
        ----------
        target_destinations : list
            Contains reagents that need to be found
        max_single_transfer_volume : number
            The chemical name (for reagents) or chemical smiles (for both reagents and products)

        Returns
        -------
        list
            Status code : str
                Either 'Error' or 'Success'
            Function return : dict or str
                Either updated destination groups or error message from the function
    """
    try:
        best_grouping = {}
        best_details = {}
        grouping_attempts = 20
        for attempt in range(grouping_attempts):
            current_grouping = {}
            if best_grouping:
                random.shuffle(target_destinations)
            for destination in target_destinations:
                if not current_grouping:
                    group_number = 1
                    current_grouping[group_number] = [destination]
                else:
                    for group_number in current_grouping.keys():
                        current_members = current_grouping[group_number]
                        if sum([member[2] for member in current_members]) + destination[2] <= max_single_transfer_volume:
                            current_grouping[group_number].append(destination)
                            current_grouping[group_number] = natsort.natsorted(current_grouping[group_number],
                                                                               key = lambda x: x[0])
                            break
                    else:
                        group_number = max(list(current_grouping.keys())) + 1
                        current_grouping[group_number] = [destination]

            number_of_groups = len(list(current_grouping.keys()))
            group_volumes = []
            for group_number in current_grouping:
                group_volumes.append(sum(item[2] for item in current_grouping[group_number]))

            if not best_grouping:
                best_grouping = current_grouping
                best_details = {'number_of_groups': number_of_groups, 'group_volumes': group_volumes}
            else:
                if number_of_groups <= best_details['number_of_groups']:
                    if np.std(group_volumes) <= np.std(best_details['group_volumes']):
                        best_grouping = current_grouping
                        best_details = {'number_of_groups': number_of_groups, 'group_volumes': group_volumes}
    except:
        return ['Error', traceback.format_exc()]
    return ['Success', best_grouping]


# -------------------
# Reference Functions
# -------------------


def previous_product_well_origin(library, product_name, product_dict, lookup_type='normal'):
    relevant_collections = ['wellplates']
    relevant_fields = ['target_product', 'chemical_smiles', 'chemical_name']
    return_loc = {}
    for collection in relevant_collections:
        for container in library[collection].keys():
            # pprint(library[collection][container])
            if lookup_type == 'characterization':
                if 'container_category' not in library[collection][container].keys():
                    continue
                if library[collection][container]['container_category'] != product_dict['container_category']:
                    continue
            elif lookup_type == 'normal':
                if library[collection][container]['container_name'] != product_dict['origin_plate']:
                    continue
            else:
                return ['Error', 'Lookup type %s is not implemented' % lookup_type]
            if library[collection][container]['location'][0] != 'liquid_handler':
                continue
            if type(library[collection][container]['contents']) == str:
                continue
            for well in library[collection][container]['contents'].keys():
                if lookup_type == 'normal':
                    if well != product_dict['origin_well_location']:
                        continue
                for relevant_field in relevant_fields:
                    # print(collection, container, well, relevant_field)
                    if relevant_field in library[collection][container]['contents'][well].keys():
                        if 'remaining' in library[collection][container]['contents'][well].keys():
                            if library[collection][container]['contents'][well]['remaining'] == 0:
                                continue
                        if relevant_field == 'target_product':
                            for sub_item in library[collection][container]['contents'][well][relevant_field]:
                                # pprint(sub_item)
                                if product_name in sub_item:
                                    key_list = [container, well, library[collection][container]['container_name'],
                                                relevant_field, sub_item]
                                    break
                            else:
                                continue
                            break
                        else:
                            if product_name in library[collection][container]['contents'][well][relevant_field]:
                                key_list = [container, well, library[collection][container]['container_name'],
                                            relevant_field]
                                break
                else:
                    continue
                break
            else:
                continue
            break
        else:
            continue
        break
    else:
        return ['Error', 'There was an issue finding the target product %s' % product_name]
    try:
        return_loc['collection'] = collection
        return_loc['relevant_field'] = key_list[3]
        return_loc['destination'] = product_dict['destination']
        return_loc['target_product'] = library[collection][key_list[0]]['contents'][key_list[1]]['target_product']
        return_loc['previous_mols'] = key_list[4]
        return_loc['container_id'] = key_list[0]
        return_loc['container_name'] = key_list[2]
        return_loc['origin_well_location'] = key_list[1]
        return_loc['previous_solvent'] = library[collection][key_list[0]]['contents'][key_list[1]]['solvents']
        return_loc['previous_volume'] = library[collection][key_list[0]]['contents'][key_list[1]]['total_volume']
        if 'properties' in library[collection][key_list[0]]['contents'][key_list[1]].keys():
            return_loc['properties'] = library[collection][key_list[0]]['contents'][key_list[1]]['properties']
    except Exception as e:
        return ['Error', traceback.format_exc()]
    return ['Success', return_loc]


# ---------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------
# -------------------------------------Potenially cut functions below--------------------------------------
# ---------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------

def directory_check(path_to_check):
    if not os.path.exists(path_to_check):
        os.makedirs(path_to_check)
        return ['Success', 'Made a new directory: %s' % path_to_check]
    return ['Success', 'Directory already exists: %s' % path_to_check]


def solvent_finder(library, solvent_name, solvent_dict):
    volume_check = [item[1] < 0 for item in solvent_dict['destination']]
    if any(volume_check):
        indices = [i for i, val in enumerate(volume_check) if val]
        well_indexes = [solvent_dict['destination'][i] for i in indices]
        return ['Error', 'Negative solvent volume, check wells: %s' % well_indexes]
    total_volume_needed = sum([float(item[1]) for item in solvent_dict['destination']])
    return_loc = {}
    relevant_collections = ['reagents', 'solvents']
    for collection in relevant_collections:
        for container in library[collection].keys():
            if library[collection][container]['location'][0] != 'liquid_handler':
                continue
            if collection == 'reagents' and library[collection][container]['labware_type'] != '8 Position Vial Carrier':
                continue
            for well in library[collection][container]['contents'].keys():
                if library[collection][container]['contents'][well]['chemical_smiles'] == 'empty':
                    continue
                elif library[collection][container]['contents'][well]['chemical_smiles'] == solvent_name or \
                        library[collection][container]['contents'][well]['chemical_name'] == solvent_name:
                    if total_volume_needed > library[collection][container]['contents'][well]['volume_ul'] * 0.95:
                        return ['Error', 'There is an issue with the solvent volume remaining: %s -> %s of %s' %
                                (solvent_name, total_volume_needed,
                                 library[collection][container]['contents'][well]['volume_ul'])]
                    return_loc['initial_solvent_details'] = copy.deepcopy(
                        library[collection][container]['contents'][well])
                    return_loc['final_solvent_details'] = copy.deepcopy(
                        library[collection][container]['contents'][well])
                    key_list = [container, well, library[collection][container]['container_name']]
                    break
            else:
                continue
            break
        else:
            continue
        break
    else:
        return ['Error', 'There was an issue finding the solvent on the liquid handler: %s' % solvent_name]
    try:
        return_loc['collection'] = collection
        return_loc['destination'] = []
        for dest_loc in solvent_dict['destination']:
            return_loc['destination'].append([dest_loc[0], float(dest_loc[1]), '', dest_loc[2]])
            return_loc['final_solvent_details']['volume_ul'] -= float(dest_loc[1])
        return_loc['container_id'] = key_list[0]
        return_loc['container_name'] = key_list[2]
        return_loc['well_location'] = key_list[1]
    except Exception as e:
        return ['Error', traceback.format_exc()]
    return ['Success', return_loc]


def chemical_well_origin(library, chemical_name, chemical_dict):
    total_moles_needed = sum([float(item[1]) for item in chemical_dict['destination']])
    return_loc = {}
    for well_plate in library['reagents'].keys():
        if library['reagents'][well_plate]['location'][0] != 'liquid_handler':
            continue
        for well in library['reagents'][well_plate]['contents'].keys():
            if library['reagents'][well_plate]['contents'][well]['chemical_smiles'] == chemical_name or \
                    library['reagents'][well_plate]['contents'][well]['chemical_name'] == chemical_name:
                initial_well_volume = library['reagents'][well_plate]['contents'][well]['volume_ul']
                concentration = library['reagents'][well_plate]['contents'][well]['concentration_molar']
                if ((float(total_moles_needed) / concentration) * 1E6) > (initial_well_volume * 0.95):
                    return ['Error', 'There is an issue with the volume remaining: %s -> %s of %s' %
                            (chemical_name, '{:.2f}'.format((float(total_moles_needed) / concentration) * 1E6),
                             str(initial_well_volume))]
                return_loc['final_well_details'] = copy.deepcopy(library['reagents'][well_plate]['contents'][well])
                key_list = [well_plate, well, library['reagents'][well_plate]['container_name']]
                break
        else:
            continue
        break
    else:
        return ['Error', 'Issue finding reagent (%s) on the liquid handler' % chemical_name]

    try:
        return_loc['destination'] = []
        return_loc['final_well_details']['volume_ul'] = 0
        for dest_loc in chemical_dict['destination']:
            volume_needed = (float(dest_loc[1]) / return_loc['final_well_details']['concentration_molar']) * 1E6
            return_loc['destination'].append([dest_loc[0], dest_loc[1], volume_needed, dest_loc[2]])
            return_loc['final_well_details']['volume_ul'] += volume_needed
        return_loc['plate_location'] = library['reagents'][well_plate]['location'][1]
        return_loc['plate_id'] = key_list[0]
        return_loc['container_name'] = key_list[2]
        return_loc['well_location'] = return_loc['final_well_details']['plate_well']
        return_loc['initial_well_details'] = copy.deepcopy(library['reagents'][well_plate]['contents'][well])
        return_loc['collection'] = 'reagents'
    except Exception as e:
        return ['Error', traceback.format_exc()]
    return ['Success', return_loc]
