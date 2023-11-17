
import random
import yaml
import datetime
import sys

from copy import deepcopy
from pprint import pprint
if __name__ == '__main__':
    from queue_generation.groups_to_queues import generate_queue_document
    sys.path.append(r'.\queue_generation')
else:
    from reaction_grouping.queue_generation.groups_to_queues import generate_queue_document


def reaction_plate_grouping(job_information, job_name, temperature_sequence, candidate_reactions, molecule_values=None):
    print('\t\tPutting reactions into the wellplate sequence')
    if molecule_values is None:
        molecule_values = {}
    print('length of molecule values is %s' % len(molecule_values))
    print('length of candidate reactions is %s' % len(candidate_reactions))
    
    # There are some key plating details we need to get from the job document
    grouping_information = job_information['grouping_information']
    temperature_spacing = grouping_information['temperature_spacing']
    if 'max_plating_iterations' in grouping_information.keys():
        plating_iterations = grouping_information['max_plating_iterations']
    else:
        plating_iterations = 100
    if 'number_of_reactions_per_plate' in grouping_information['uncertainty'].keys():
        n_reactions = grouping_information['uncertainty']['number_of_reactions_per_plate']
    else:
        if grouping_information['uncertainty']['method'] == 'exploitation':
            n_reactions = 12
        else:
            n_reactions = 96

    # Next we can start the grouping
    best_grouping = {'grouping': {}, 'n_products': 0, 'product_value': 0, 'products': set()}
    for index, item in enumerate(temperature_sequence):
        best_grouping['grouping'][index] = {'reactions': [], 'temperature': item}
    products_to_consider = list(candidate_reactions.keys())

    reaction_templates = set()
    for iteration in range(0, plating_iterations):
        random_product_sequence = deepcopy(products_to_consider)
        random.shuffle(random_product_sequence)
        current_grouping = {'grouping': {}, 'n_products': 0, 'product_value': 0, 'products': set()}
        for index, item in enumerate(temperature_sequence):
            current_grouping['grouping'][index] = {'reactions': [], 'temperature': item}
        for product_key in random_product_sequence:
            for reaction_tree_key in candidate_reactions[product_key].keys():
                reaction_tree = candidate_reactions[product_key][reaction_tree_key]
                product_temperatures = []
                tree_products = set()
                bad_tree = 0
                for reaction_index, _ in enumerate(reaction_tree.keys()):
                    reaction_key = 'reaction_%s' % (reaction_index + 1)
                    if reaction_key not in reaction_tree.keys():
                        bad_tree = 1
                        break
                    tree_products.update(reaction_tree[reaction_key]['products'])
                if bad_tree == 1:
                    print(product_key, reaction_tree_key, list(reaction_tree.keys()))
                    continue
                for reaction_index, _ in enumerate(reaction_tree.keys()):
                    reaction_key = 'reaction_%s' % (reaction_index + 1)
                    for context_key in reaction_tree[reaction_key]['contexts'].keys():
                        context = reaction_tree[reaction_key]['contexts'][context_key]
                        context_temperature = temperature_spacing * round(context['temperature'] / temperature_spacing)
                        if context_temperature in temperature_sequence:
                            break
                    else:
                        continue
                    product_temperatures.append([product_key, reaction_tree_key, reaction_key,
                                                 len(list(reaction_tree.keys())),
                                                 context_key, context_temperature])
                if len(product_temperatures) != len(list(reaction_tree.keys())):
                    print(product_key, reaction_tree_key)
                    continue
                space_available = []
                reaction_index = 0
                for group_index in range(0, len(temperature_sequence)):
                    if len(product_temperatures) <= reaction_index:
                        break
                    if len(product_temperatures[reaction_index]) == 0:
                        continue
                    current_group_temperature = current_grouping['grouping'][group_index]['temperature']
                    if product_temperatures[reaction_index][-1] == current_group_temperature:
                        if len(current_grouping['grouping'][group_index]['reactions']) <= n_reactions - 1:
                            space_available.append([group_index, product_temperatures[reaction_index]])
                            reaction_index += 1
                        if reaction_index == len(product_temperatures):
                            break
                if len(space_available) != len(product_temperatures):
                    continue
                for space in space_available:
                    current_grouping['grouping'][space[0]]['reactions'].append(space[1])
                current_grouping['n_products'] += 1
                if molecule_values:
                    if product_key in molecule_values.keys():
                        prediction_value = molecule_values[product_key]['prediction_value']
                        uncertainty_value = molecule_values[product_key]['uncertainty_value']
                        current_grouping['product_value'] += (prediction_value * uncertainty_value)
                else:
                    current_grouping['product_value'] += 1
                current_grouping['products'].add(product_key)
        print_update = False
        if current_grouping['product_value'] > best_grouping['product_value']:
            best_grouping = deepcopy(current_grouping)
            print_update = True
        # There was an idea to use the maximum number of wells on the well-plate, removed this for now
        #elif current_grouping['product_value'] == best_grouping['product_value']:
        #    best_used_wells = 0
        #    for key in best_grouping['grouping'].keys():
        #        best_used_wells += len(best_grouping['grouping'][key]['reactions'])
        #    test_used_wells = 0
        #    for key in current_grouping['grouping'].keys():
        #        test_used_wells += len(best_grouping['grouping'][key]['reactions'])
        #    if test_used_wells > best_used_wells:
        #        best_grouping = deepcopy(current_grouping)
        #        print_update = True
        if print_update:
            print('\t\t\tIteration: %s, %s products, %s best products, %s/%s value' % (iteration,
                                                                                       current_grouping['n_products'],
                                                                                       best_grouping['n_products'],
                                                                                       current_grouping['product_value'],
                                                                                       best_grouping['product_value']))
        reaction_templates = set()
        for group in best_grouping['grouping'].keys():
            for reaction in best_grouping['grouping'][group]['reactions']:
                if 'reaction_details' not in best_grouping['grouping'][group].keys():
                    best_grouping['grouping'][group]['reaction_details'] = {}
                group_reaction_details = best_grouping['grouping'][group]['reaction_details']
                current_reaction = candidate_reactions[reaction[0]]
                current_reaction_details = current_reaction[reaction[1]][reaction[2]]
                if reaction[4] not in current_reaction_details['contexts'].keys():
                    pprint(current_reaction_details)
                current_reaction_context = current_reaction_details['contexts'][reaction[4]]
                reaction_templates.update(current_reaction_details['templates'])
                group_reaction_key = '_'.join([str(item) for item in reaction])
                group_reaction_details[group_reaction_key] = deepcopy(current_reaction_details)
                group_reaction_details[group_reaction_key]['contexts'] = deepcopy(current_reaction_context)
                tree_depth = len(list(current_reaction[reaction[1]].keys()))
                group_reaction_details[group_reaction_key]['tree_depth'] = tree_depth

    print('------------------------------')
    print('The queue will use these reaction templates: ')
    pprint(reaction_templates)
    print('------------------------------')
    return ['Success', best_grouping]


def reaction_grouping_to_queue_input(job_information, job_name, reaction_grouping):
    grouping_information = job_information['grouping_information']
    queue_details = grouping_information['queue_details']
    workflow_type = queue_details['workflow_type']
    reaction_scale = queue_details['preparation_scale']
    current_time = datetime.datetime.now()
    campaign_name = '%s' % job_name
    if workflow_type == 'standard':
        workflow_details = {'reaction_time': 60 * 60 * 8,
                            'intermediate_purification': True,
                            'hplc_analysis': True,
                            'separate_hplc_plate': True,
                            'workflow_type': 'reaction'}

        queue_input = {'campaign_name': campaign_name,
                       'workflow_type': workflow_type,
                       'wellplate_type': '',
                       'wellplate_grouping': {}}
        collected_products = set()
        for group_key in reaction_grouping['grouping'].keys():
            reaction_group = reaction_grouping['grouping'][group_key]
            reaction_temperature = reaction_group['temperature']
            plated_reaction_details = {}
            plate_chemicals = {'reagents': set(), 'catalysts': set(), 'solvents': set()}
            for reaction_index, reaction in enumerate(reaction_group['reactions']):
                group_reaction_details = reaction_group['reaction_details']['_'.join([str(item) for item in reaction])]
                reaction_details = {'reactants': [],
                                    'reaction_smiles': group_reaction_details['reaction'],
                                    'templates': group_reaction_details['templates']}
                tree_depth = group_reaction_details['tree_depth']
                reaction_number = int(reaction[2].split('_')[1])
                reactions_remaining = tree_depth - reaction_number
                reactant_scale = (tree_depth - reaction_number + 1) * reaction_scale
                for reactant in group_reaction_details['reactants']:
                    reaction_details['reactants'].append([reactant, reactant_scale, 1])
                    plate_chemicals['reagents'].add(reactant)
                reaction_context = group_reaction_details['contexts']
                if 'reagent' in reaction_context.keys() and len(reaction_context['reagent']) > 0:
                    reaction_details['reagents'] = []
                    for reagent in reaction_context['reagent']:
                        reaction_details['reagents'].append([reagent, reactant_scale, 2])
                        plate_chemicals['reagents'].add(reagent)
                if 'catalyst' in reaction_context.keys() and len(reaction_context['catalyst']) > 0:
                    reaction_details['catalysts'] = []
                    catalyst_scale = reactant_scale * 0.1
                    for catalyst in reaction_context['catalyst']:
                        reaction_details['catalysts'].append([catalyst, catalyst_scale, 3])
                        plate_chemicals['catalysts'].add(catalyst)
                reaction_details['solvents'] = []
                if 'solvent' in reaction_context.keys():
                    for solvent in reaction_context['solvent']:
                        reaction_details['solvents'].append([solvent, 1, 1])
                        plate_chemicals['solvents'].add(solvent)
                else:
                    solvent = 'dmso'
                    reaction_details['solvents'].append([solvent, 1, 1])
                    plate_chemicals['solvents'].add(solvent)
                if reactions_remaining == 0:
                    collection = 'yes'
                else:
                    collection = 'no'
                if collection == 'yes':
                    collected_products.add(reaction[0])
                reaction_details['final_product'] = [reactions_remaining, collection, reaction[0], tree_depth]
                reaction_details['total_volume'] = 500

                reaction_key = 'reaction_%s' % (reaction_index + 1)
                plated_reaction_details[reaction_key] = reaction_details

            reaction_plate_name = 'reaction_plate_%s' % (len(queue_input['wellplate_grouping'].keys()) + 1)
            group_details = deepcopy(workflow_details)
            group_details['reaction_temperature'] = reaction_temperature
            extra_conditions = ''
            for catalyst in plate_chemicals['catalysts']:
                if 'pd' in catalyst.lower() or 'cu' in catalyst.lower():
                    extra_conditions = 'inert_atmosphere'
                    break
            else:
                if reaction_temperature < 30:
                    extra_conditions = 'low_temperature'
            if extra_conditions != '':
                group_details['extra_conditions'] = extra_conditions
            if extra_conditions == 'inert_atmosphere' or reaction_temperature > 100:
                wellplate_type = 'Paradox Thermal Plate'
            else:
                wellplate_type = '96 Well DeepWell'
            group_details['wellplate_type'] = wellplate_type
            queue_input['wellplate_grouping'][reaction_plate_name] = {'workflow_details': group_details,
                                                                      'reaction_grouping': plated_reaction_details}
        if 'characterization_needed' in job_information.keys():
            characterization_details = job_information['characterization_needed']
            for characterization_plate in characterization_details.keys():
                plate_details = characterization_details[characterization_plate]
                workflow_details = {'characterizations': deepcopy(plate_details['measurements']),
                                    'workflow_type': 'characterization',
                                    'wellplate_type': '96 Well Microplate'}
                group_details = {}
                for collected_product in collected_products:
                    well_details = {'solvents': plate_details['plate_solvent'],
                                    'target_molecules': [[collected_product, 1, 'OD', 'fraction_plate', 1]],
                                    'total_volume': 300}
                    product_key = len(list(group_details.keys())) + 1
                    group_details[product_key] = well_details
                queue_input['wellplate_grouping'][characterization_plate] = {'workflow_details': workflow_details,
                                                                             'reaction_grouping': group_details}
    else:
        return ['Error', 'Grouping type %s is not implemented' % workflow_type]
    return ['Success', queue_input]


def build_queue_documents(job_information, job_name, reaction_grouping, output_directory):

    return_statement, grouping_document = reaction_grouping_to_queue_input(job_information, job_name, reaction_grouping)
    if return_statement != 'Success':
        return ['Error', grouping_document]

    return_statement, return_object = generate_queue_document(output_directory, grouping_document, {})
    if return_statement != 'Success':
        return ['Error', return_object]

    return ['Success', 'Queue documents were generated for %s' % job_name]


if __name__ == '__main__':
    print('Loading reactions!')
    with open(r'D:\Python\prediction-pipeline-development\reaction_grouping.yaml', 'r') as yaml_file:
        candidate_reactions = yaml.load(yaml_file, Loader=yaml.Loader)
    print('Loading Job Information!')
    with open(r'D:\Python\prediction-pipeline-development\input_jobs\tztz_scaffold_20230406_test.yaml', 'r') as yaml_file:
        job_information = yaml.load(yaml_file, Loader=yaml.Loader)
    job_name = 'tztz_scaffold_20230406_test'

    return_statement, grouping_document = reaction_grouping_to_queue_input(job_information, job_name, candidate_reactions)
    print(return_statement)
    #pprint(grouping_document)
    output_directory = r'D:\Python\prediction-pipeline-development\output_jobs\tztz_scaffold_20230406_test'

    return_statement, return_object = generate_queue_document(output_directory, grouping_document, {})
    print(return_statement)
    print(return_object)
