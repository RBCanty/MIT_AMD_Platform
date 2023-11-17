# -*- coding: utf-8 -*-
"""
Updated on 4/8/2022

@author: Brent Koscher
"""

import json
import re
import traceback
from pprint import pprint


def initial_worktable_prep(carriers_labware_filepath, empty_table_filepath, mongo_library):
    """
    Parses the empty Evoware Worktable to start building a schedule around

    Parameters
    ----------
    carriers_labware_filepath : string
        Path to the processed carriers and labware library (JSON file)
    empty_table_filepath : string
        Path to the empty Evoware worktable (EWT file)
    mongo_library : dictionary
        Dictionary of the current platform library

    Returns
    -------
    return_statement : string
        Short string to check later for error handling purposes
    start_carriers : dictionary
        Processed platform carrier dictionary of grid -> carrier -> details
    carriers_labware : dictionary
        Carrier and labware dictionary that is used by some other functions

    """
    try:
        # We will need the carriers and labware dictionary to reference for validating
        # the worktable carriers and labware as we build it
        with open(carriers_labware_filepath, 'r') as jsonfile:
            carriers_labware = json.load(jsonfile)

        # We will also need an empty worktable to reference as we start constructing
        # the output generated schedule
        with open(empty_table_filepath, 'r') as infile:
            empty_worktable = infile.readlines()

        # Now we can grab the line that specifies all of the carriers on the bed
        carrier_wt_list = empty_worktable[8].strip('\n').split(';')[0:-1]
        carrier_wt_list.pop(0)

        # With the carriers list we can look them all up in the carrier library
        start_carriers = {}
        for grid, carrier in enumerate(carrier_wt_list):
            # Carrier 14 is the "System" carrier, not important to process
            if carrier == '14':
                continue

            # Carrier -1 is an empty carrier grid, not important to process
            if carrier == '-1':
                continue

            # The number of positions and carrier name will be used when constructing
            # the carrier dictionary to build
            num_positions = carriers_labware['carriers'][carrier]['positions']
            carrier_name = carriers_labware['carriers'][carrier]['name']

            # Special exception for 2 Pos + Waste Carrier for MCA Tips, the tips
            # are managed by the database in a more straightforward method
            if carrier == '328':
                labware_labels = [''] * num_positions
                labware_types = ['DiTi 200ul Nested MCA96']*16 + ['DiTi Nested Waste'] + ['']*7 + ['DiTi SBS Waste']*2 + ['']

            # Special exception for wash station for the LiHa
            elif carrier == '30':
                labware_labels = [''] * num_positions
                labware_types = ['Wash Station Cleaner shallow', 'Wash Station Waste', 'Wash Station Cleaner deep']

            # Special exception for the TeVaC Custom Carrier (modified Tecan default)
            elif carrier == '333':
                labware_labels = ['SepPlateR', 'FiltrateR', 'BlockR', 'SepPlateF',
                                  'FiltrateF', 'BlockF', 'BlockFtoR', 'BlockRtoF']
                labware_types = ['96 Well Separation Plate', '96 Well Microplate', 'TeVac Block', '96 Well Separation Plate',
                                '96 Well Microplate', 'TeVac Block', 'TeVac Block', 'TeVac Block']

            # These are the other non-special carriers
            else:
                labware_labels = [''] * num_positions
                labware_types = [''] * num_positions

            # We add the carrier to a dict to keep track of locations for scripts
            if grid not in start_carriers.keys():
                start_carriers[grid] = {}
            start_carriers[grid][int(carrier)] = {'carrier_id': int(carrier), 'num_positions': num_positions,
                                                  'labware_labels': labware_labels, 'labware_types': labware_types,
                                                  'carrier_name': carrier_name}

        # Now we can move onto the off-grid carriers like hotels. These happen after
        # a line has been given to every carrier grid on the worktable and this changes
        # depending on the model of Tecan that is being used
        num_carrier_slots = len(carrier_wt_list)
        line_index = 1
        off_grid_start_index = None
        for sub_index, element in enumerate(empty_worktable):
            if line_index == num_carrier_slots:
                off_grid_start_index = sub_index + 1
                break
            if re.match(r'998;\d{1,3};', element):
                line_index += 1
        else:
            return ['Error', 'Empty worktable file is incorrectly structured']

        # The location of the off-grid carriers need to be parsed, with the
        # off-grid occupied carriers appearing first in the worktable
        hotel_number = int(empty_worktable[off_grid_start_index].split(';')[1])
        off_grid_end_index = off_grid_start_index + hotel_number + 1
        off_grid_hotels = [line.strip('\n') for line in empty_worktable[off_grid_start_index+1:off_grid_end_index]]

        # Then all of the off-grid carriers appear in a block with their names,
        # followed by a block with their carrier grid numbers
        off_grid_number = int(empty_worktable[off_grid_end_index].split(';')[1])
        off_grid_carriers = [line.strip('\n') for line in empty_worktable[off_grid_end_index+1:off_grid_end_index + off_grid_number + 1]]
        off_grid_carrier_divider = off_grid_end_index + off_grid_number + 2
        off_grid_locations = [line.strip('\n') for line in empty_worktable[off_grid_carrier_divider:off_grid_carrier_divider + off_grid_number]]

        # Then with those two blocks in hand we can parse them into our carrier
        # dictionary that the script will use to generate schedules
        for i in range(1, len(off_grid_locations)):
            location = int(off_grid_locations[-i].split(';')[1])
            carrier_name = off_grid_carriers[-i].split(';')[3]

            # Since we only have the name of the carrier we need to find it in
            # the carriers and labware dictionary, these are uniquely defined
            # because of the needs of Evoware
            for carrier in carriers_labware['carriers'].keys():
                if carrier_name == carriers_labware['carriers'][carrier]['name']:
                    carrier_key = carrier
                    break
            else:
                return ['Error', 'Labware and carrier values are not defined yet for %s' % off_grid_carriers[-i], '']
            if location not in start_carriers.keys():
                start_carriers[location] = {}

            # With the carrier key, we have access to the carrier information
            carrier_details = carriers_labware['carriers'][carrier_key]
            start_carriers[location][int(carrier_key)] = {'carrier_id': int(carrier_key),
                                                          'num_positions': carrier_details['positions'],
                                                          'labware_labels': [''] * carrier_details['positions'],
                                                          'labware_types': [''] * carrier_details['positions'],
                                                          'carrier_name': carrier_name}

        # With the carriers defined and in place, the bed can be populated
        relevant_collections = ['solvents', 'reagents', 'wellplates', 'consumables']
        for relevant_collection in relevant_collections:
            for container_key in mongo_library[relevant_collection].keys():
                relevant_container = mongo_library[relevant_collection][container_key]
                if relevant_container['container_name'] == 'DiTi_200uL_Nested_328': # Special exception
                    continue
                container_location = relevant_container['location']
                if container_location[0] != 'liquid_handler':
                    continue
                bed_location = container_location[1]
                grid = bed_location[0]
                carrier_site = bed_location[1]
                carrier_id = bed_location[2]
                start_carriers[grid][carrier_id]['labware_labels'][carrier_site - 1] = relevant_container['container_name']
                start_carriers[grid][carrier_id]['labware_types'][carrier_site - 1] = relevant_container['labware_type']

    except Exception:
        pprint(start_carriers)
        return ['Error', traceback.format_exc(), '']
    return ['Success', start_carriers, carriers_labware]


def worktable_export(carriers_labware_filepath, empty_table_filepath, updated_carriers, full_library_mongo, simplified_schedule):
    """
    Reparses the empty worktable to be consistent with any changes to the worktable
    that arise during the method execution and generates the strings to write

    Parameters
    ----------
    carriers_labware_filepath : string
        Path to the processed carriers and labware library (JSON file)
    empty_table_filepath : string
        Path to the empty Evoware worktable (EWT file)
    updated_carriers : dictionary
        Dictionary of the current carriers on the bed
    full_library_mongo : dictionary
        Dictionary of the current platform library
    simplified_schedule : list of lists
        Building blocks needed to build method strings for populating locations

    Returns
    -------
    return_statement : string
        Short string to check later for error handling purposes
    output_schedule_block : list
        List of strings that define the carriers and labware block of schedules

    """
    try:
        # We will need the carriers and labware dictionary to reference for
        # validating the worktable carriers and labware as we build it
        with open(carriers_labware_filepath, 'r') as jsonfile:
            carriers_labware = json.load(jsonfile)

        # We will also need an empty worktable to reference as we start
        # constructing the output generated schedule
        with open(empty_table_filepath, 'r') as infile:
            empty_worktable = infile.readlines()

        # The updated carriers does not include the liquid handler locations
        # that are transiently occupied by labware on the bed, these locations
        # need to be populated to avoid Evoware errors
        for task in simplified_schedule:
            method_name = task[0]
            if method_name in ['TRANSFER_LABWARE', 'SPECIAL_LID_TRANSFER']:
                source = task[2]
                source_location = updated_carriers[source[0]][source[2]]
                source_location['labware_labels'][source[1] - 1] = 'location_%s_%s_%s' % (source[0], source[1], source[2])
                source_location['labware_types'][source[1] - 1] = full_library_mongo[task[-1]][str(task[-2])]['labware_type']
                destination = task[3]
                destination_location = updated_carriers[destination[0]][destination[2]]
                destination_location['labware_labels'][destination[1] - 1] = 'location_%s_%s_%s' % (destination[0], destination[1], destination[2])
                destination_location['labware_types'][destination[1] - 1] = full_library_mongo[task[-1]][str(task[-2])]['labware_type']

            # This catches non-database defined transfers that are needed
            if method_name == 'SPECIAL_TRANSFER_LABWARE':
                source = task[2]
                source_location = updated_carriers[source[0]][source[2]]
                source_location['labware_labels'][source[1] - 1] = 'location_%s_%s_%s' % (source[0], source[1], source[2])
                source_location['labware_types'][source[1] - 1] = task[1]
                destination = task[3]
                destination_location = updated_carriers[destination[0]][destination[2]]
                destination_location['labware_labels'][destination[1] - 1] = 'location_%s_%s_%s' % (destination[0], destination[1], destination[2])
                destination_location['labware_types'][destination[1] - 1] = task[1]

        # We need the admin lines at the start and end of the worktable
        admin_lines = empty_worktable[0:9]
        end_admin_lines = [line.strip('\n') for line in empty_worktable[-2:]]
        
        # Now we can grab the line that specifies all of the carriers on the bed
        carrier_wt_list = empty_worktable[8].strip('\n').split(';')[0:-1]
        carrier_wt_list.pop(0)

        # With the carriers list we can look them all up in the carrier library
        on_plat_carriers = []
        for line_index, carrier in enumerate(carrier_wt_list):
            if carrier == '14': # Carrier 14 is the "System" carrier, not important
                continue

            if carrier == '-1': # Empty carrier grid, add the filler line
                on_plat_carriers.append(['998;0;'])
                continue

            if carrier == '328': # 2 Pos + Waste Carrier for MCA Tips
                on_plat_carriers.append([''.join(['998;27;'] + ['DiTi 200ul Nested MCA96;'] * 16
                                                 + ['DiTi Nested Waste;'] + [';'] * 7 + ['DiTi SBS Waste;'] * 2
                                                 + [';']), ''.join(['998;'] + [';'] * 27)])

            elif carrier == '30': # Wash station for the LiHa
                on_plat_carriers.append(['998;3;Wash Station Cleaner shallow;Wash Station Waste;Wash Station Cleaner deep;', '998;;;;'])

            elif carrier == '333': # TeVaC Custom Carrier (modified Tecan default)
                on_plat_carriers.append(['998;10;96 Well Separation Plate;96 Well Microplate;TeVac Block;96 Well Separation Plate;96 Well Microplate;TeVac Block;TeVac Block;TeVac Block;;;',
                                         '998;SepPlateR;FiltrateR;BlockR;SepPlateF;FiltrateF;BlockF;BlockFtoR;BlockRtoF;;;'])
                
            else: # These are the other non-special carriers
                num_positions = carriers_labware['carriers'][carrier]['positions']
                on_plat_carriers.append([['998;%s;' % str(num_positions)] + [';'] * num_positions,
                                         ['998;'] + [';'] * num_positions])

        # Now we can move onto the off-grid carriers like hotels. These happen after
        # a line has been given to every carrier grid on the worktable and this changes
        # depending on the model of Tecan that is being used
        num_carrier_slots = len(carrier_wt_list)
        line_index = 1
        off_grid_start_index = None
        for sub_index, element in enumerate(empty_worktable):
            if line_index == num_carrier_slots:
                off_grid_start_index = sub_index + 1
                break
            if re.match(r'998;\d{1,3};', element):
                line_index += 1
        else:
            return ['Error', 'Empty worktable file is incorrectly structured']

        # The location of the off-grid carriers need to be parsed, with the
        # off-grid occupied carriers appearing first in the worktable
        hotel_number = int(empty_worktable[off_grid_start_index].split(';')[1])
        off_grid_end_index = off_grid_start_index + hotel_number + 1
        off_grid_hotels = [line.strip('\n') for line in empty_worktable[off_grid_start_index+1:off_grid_end_index]]

        # Then all of the off-grid carriers appear in a block with their names,
        # followed by a block with their carrier grid numbers
        off_grid_number = int(empty_worktable[off_grid_end_index].split(';')[1])
        off_grid_carriers = [line.strip('\n') for line in empty_worktable[off_grid_end_index+1:off_grid_end_index + off_grid_number + 1]]
        off_grid_carrier_divider = off_grid_end_index + off_grid_number + 2
        off_grid_locations = [line.strip('\n') for line in empty_worktable[off_grid_carrier_divider:off_grid_carrier_divider + off_grid_number]]

        # Anything that is occupied because of the script gets updated starting
        # with the on platform carriers first
        for grid in updated_carriers.keys():
            for carrier_type in updated_carriers[grid].keys():
                if carrier_type in [30, 328, 333] : # Ignore the special carriers
                    continue
                if str(carrier_type) in carrier_wt_list:
                    on_plat_carriers[grid][0] = ';'.join(['998;%s' % str(len(updated_carriers[grid][carrier_type]['labware_types']))] +
                        updated_carriers[grid][carrier_type]['labware_types']) + ';'
                    on_plat_carriers[grid][1] = ';'.join(['998'] + updated_carriers[grid][carrier_type]['labware_labels']) + ';'

        # Now we need to keep track of the hotels that are occupied by labware,
        # this can only be a single type of labware, first we
        occupied = []
        hotel_problems = []
        for hotel in off_grid_hotels:
            hotel_split = hotel.split(';')
            occupied.append([int(hotel_split[2]), int(hotel_split[1]), ''])
        for location_index, location in enumerate(occupied):
            # Hotels that are double occupied are problematic to reliable operation,
            # so we check for more than 2 set items (one unique labware and empty string)
            if len(list(set(updated_carriers[location[0]][location[1]]['labware_types']))) > 2:
                hotel_problems.append('Too many types in hotel %s, %s types' % (location[0:2],
                                      len(list(set(updated_carriers[location[0]][location[1]]['labware_types'])))))
            for wellplate_site in updated_carriers[location[0]][location[1]]['labware_types']:
                # Here any non-empty string is a type of labware
                if wellplate_site != '':
                    wellplate_type = wellplate_site
                    occupied[location_index][2] = wellplate_site
                    break
            else:
                continue
            # Now we need to check all of the hotels that have the same hotel type
            # due to Evoware they must hold the same type of wellplate
            for sub_location_index, sub_location in enumerate(occupied):
                if location_index != sub_location_index:
                    if sub_location[1] == location[1]:
                        if sub_location[2] not in ['', wellplate_type]:
                            hotel_problems.append('Mismanaged wellplate types in %s and %s for plate type %s' %(
                                location, sub_location, wellplate_type))
                        occupied[sub_location_index][2] = wellplate_type
        # Chedk to make sure that there are no hotel based problems in the worktable,
        # these will cause un-predictable errors in Evoware
        if len(hotel_problems) != 0:
            return ['Error', hotel_problems]

        # Now make the worktable lines to get ready to write the schedule
        off_grid_lines = []
        for grid in updated_carriers.keys():
            for carrier_type in updated_carriers[grid].keys():
                if carrier_type not in [322, 323]:
                    continue
                wellplate_type = None
                for wellplate in updated_carriers[grid][carrier_type]['labware_types']:
                    if wellplate != '':
                        wellplate_type = wellplate
                if wellplate_type is None:
                    continue
                off_grid_lines.append('998;%s;%s;' % (carrier_type, wellplate_type))
                
        off_grid_lines.extend(['998;%s;%s;' % (off_grid_site[1], off_grid_site[2]) for off_grid_site
                          in occupied if off_grid_site[2] != ''])
        on_plat_carrier_lines = []
        for carrier in on_plat_carriers:
            if len(carrier) == 1:
                on_plat_carrier_lines.append(carrier[0])
            else:
                on_plat_carrier_lines.append(''.join(carrier[0]).strip('\n'))
                on_plat_carrier_lines.append(''.join(carrier[1]).strip('\n'))

        # The output is then ordered for Evoware to read: admin lines -> on platform
        # carriers -> off grid hotels -> off grid carriers -> off grid occupied
        # lines -> off grid locations -> end admin lines
        output_schedule_block = [line.strip('\n') for line in admin_lines] + \
                    on_plat_carrier_lines + ['998;%s;' % len(off_grid_hotels)] + off_grid_hotels + \
                    ['998;%s;' % len(off_grid_carriers)] + off_grid_carriers + \
                    ['998;%s;' % len(off_grid_lines)] + off_grid_lines + \
                    off_grid_locations + end_admin_lines

    except Exception:
        return ['Error', traceback.format_exc()]
    return ['Success', output_schedule_block]
