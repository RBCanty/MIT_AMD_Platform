
import yaml
import json
import os
import sys
import queue
import subprocess
import threading
import itertools
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from pprint import pprint
from functions.supporting_functions import print_progress_bar, directory_check


def run_chemprop_models(job_information, job_name, candidate_filepath, candidate_type):
    chemprop_details = job_information['chemprop_details']
    if len(list(chemprop_details['constraining_models'])) == 0:
        return ['Success', 'No Chemprop models were selected']

    # Start by checking/making an output Chemprop folder
    chemprop_model_output = r'.\output_jobs\%s\property_predictions' % job_name
    return_statement, print_statement = directory_check(chemprop_model_output)
    if return_statement != 'Success':
        return ['Error', print_statement]
    # Next check/make a Chemprop status file
    chemprop_status_filepath = r'%s\%s_chemprop_status.yaml' % (chemprop_model_output, job_name)
    if os.path.exists(chemprop_status_filepath):
        with open(chemprop_status_filepath, 'r') as yaml_file:
            chemprop_status = yaml.load(yaml_file, Loader=yaml.Loader)
    else:
        chemprop_status = {}
    chemprop_smiles_filepath = r'%s\%s_chemprop_smiles.yaml' % (chemprop_model_output, job_name)
    if os.path.exists(chemprop_smiles_filepath):
        with open(chemprop_smiles_filepath, 'r') as yaml_file:
            chemprop_smiles = yaml.load(yaml_file, Loader=yaml.Loader)
    else:
        chemprop_smiles = {}

    # Next we need to figure out which SMILES that have input files
    if candidate_type == 'reactions':
        with open(candidate_filepath, 'r') as json_file:
            reactions = json.load(json_file)
        candidates = list(reactions.keys())
    elif candidate_type == 'smiles':
        with open(candidate_filepath, 'r') as yaml_file:
            smiles_status = yaml.load(yaml_file, Loader=yaml.Loader)
        candidates = list(smiles_status['smiles_set'])
    else:
        candidates = []

    # Now we will prepare the Chemprop input files for each of the models
    chemprop_predictions_folder = r'%s\chemprop_files' % chemprop_model_output
    return_statement, print_statement = directory_check(chemprop_predictions_folder)
    if return_statement != 'Success':
        return ['Error', print_statement, '', '']
    group_size = 10000
    for constraining_model in chemprop_details['constraining_models'].keys():
        model_details = chemprop_details['constraining_models'][constraining_model]
        if constraining_model not in chemprop_smiles.keys():
            chemprop_smiles[constraining_model] = []
        if constraining_model not in chemprop_status.keys():
            chemprop_status[constraining_model] = {'input_files': [], 'completed_files': [], 'output_files': []}
        new_smiles = list(set(candidates).difference(set(chemprop_smiles[constraining_model])))
        print('\t\t\tNeed to make input files for %s new SMILES' % len(new_smiles))
        if len(new_smiles) != 0:
            smiles_groups = [new_smiles[i:i + group_size] for i in range(0, len(new_smiles), group_size)]
            for file_index, smiles_group in enumerate(smiles_groups):
                number_of_input_files = len(chemprop_status[constraining_model]['input_files'])
                filepath = r'%s\%s_%s_%s_input.csv' % (chemprop_predictions_folder, job_name, constraining_model,
                                                       number_of_input_files+1)
                with open(filepath, 'w') as outfile:
                    # Write the header line
                    if constraining_model == 'uv_vis':
                        outfile.write('smiles,solvent\n')
                    else:
                        outfile.write('smiles\n')
                    # Next write all SMILES to the input files
                    for smiles in smiles_group:
                        if constraining_model == 'uv_vis':
                            outfile.write('%s,%s\n' % (smiles, model_details['solvent']))
                        else:
                            outfile.write('%s\n' % smiles)
                # Next update the Chemprop status tracker
                chemprop_status[constraining_model]['input_files'].append(filepath)
                chemprop_smiles[constraining_model].extend(smiles_group)
            with open(chemprop_status_filepath, 'w') as yaml_file:
                yaml.dump(chemprop_status, yaml_file)
            with open(chemprop_smiles_filepath, 'w') as yaml_file:
                yaml.dump(chemprop_smiles, yaml_file)

    # With the input files written we can start to run Chemprop models
    new_chemprop_data = False
    unavailable_models = [False, []]
    chemprop_queue = queue.Queue()
    chemprop_queue.put({'running_jobs': [], 'error_state': []})
    for constraining_model in chemprop_status.keys():
        if constraining_model not in chemprop_details['constraining_models'].keys():
            continue
        model_files = chemprop_status[constraining_model]
        files_to_run = list(set(model_files['input_files']).difference(set(model_files['output_files'])))
        model_details = chemprop_details['constraining_models'][constraining_model]
        for file_to_run in files_to_run:
            output_filepath = file_to_run.replace('input', 'output')
            if output_filepath in model_files['output_files']:
                continue
            model_version = model_details['model_version']
            print(file_to_run)
            # Depending on the model, the Chemprop argument is slightly different
            if constraining_model == 'uv_vis':
                if 'special_attribute' in model_details.keys() and model_details['special_attribute'] == 'azo_molecules':
                    model_dir = r'.\chemprop_predicting\chemprop_models\uv_vis'
                    model_checkpoint_dir = model_dir
                    if not os.path.exists(model_dir):
                        unavailable_models[0] = True
                        unavailable_models[1].append([constraining_model, model_version])
                        continue
                    chemprop_input_arguments = {'function': 'bash',
                                                'arguments': [os.path.abspath(r'%s\uv_vis_run_azo.sh' % model_dir),
                                                              os.path.abspath(file_to_run),
                                                              os.path.abspath(output_filepath),
                                                              model_version]}
                else:
                    model_dir = r'.\chemprop_predicting\chemprop_models\uv_vis\%s' % model_version
                    model_checkpoint_dir = r'%s\fold_0' % model_dir
                    chemprop_input_arguments = {'function': 'chemprop_predict',
                                                'arguments': {'test_path': file_to_run,
                                                              'preds_path': output_filepath,
                                                              'checkpoint_dir': model_checkpoint_dir,
                                                              'number_of_molecules': 2,
                                                              'ensemble_variance': ''}}
            elif constraining_model == 'logp':
                model_dir = r'.\chemprop_predicting\chemprop_models\logp'
                model_checkpoint_dir = model_dir
                if not os.path.exists(model_dir):
                    unavailable_models[0] = True
                    unavailable_models[1].append([constraining_model, model_version])
                    continue
                chemprop_input_arguments = {'function': 'bash',
                                            'arguments': [os.path.abspath(r'%s\predict_from_model.sh' % model_dir),
                                                          os.path.abspath(file_to_run),
                                                          os.path.abspath(r'%s\%s' % (model_dir, model_version)),
                                                          os.path.abspath(output_filepath)]}
            elif constraining_model == 'photodeg':
                model_dir = r'.\chemprop_predicting\chemprop_models\photodeg_rf'
                model_checkpoint_dir = r'%s\%s' % (model_dir, model_details['model_version'])
                bash_file = r'%s\predict_from_model.sh' % model_dir
                chemprop_input_arguments = {'function': 'bash',
                                            'arguments': [os.path.abspath(bash_file),
                                                          os.path.abspath(file_to_run),
                                                          os.path.abspath(output_filepath),
                                                          model_details['model_version']]}
            elif constraining_model == 'hyperpolarizability':
                model_dir = r'.\chemprop_predicting\chemprop_models\%s' % model_details['model_version']
                model_checkpoint_dir = r'%s\fold_0' % model_dir
                chemprop_input_arguments = {'function': 'chemprop_predict',
                                            'arguments': {'test_path': file_to_run,
                                                          'preds_path': output_filepath,
                                                          'checkpoint_dir': model_checkpoint_dir,
                                                          'ensemble_variance': ''}}
            elif constraining_model == 'aromatase_inhibition':
                model_dir = r'.\chemprop_predicting\chemprop_models\%s' % model_details['model_version']
                model_checkpoint_dir = model_dir
                bash_file = r'%s\aromatase_inhibs_predict.sh' % model_dir
                chemprop_input_arguments = {'function': 'bash',
                                            'arguments': [os.path.abspath(bash_file),
                                                          os.path.abspath(file_to_run),
                                                          os.path.abspath(output_filepath)]}
            elif constraining_model == 'aromatase_isomer_docking':
                model_dir = r'.\chemprop_predicting\chemprop_models\%s' % model_details['model_version']
                model_checkpoint_dir = model_dir
                bash_file = r'%s\isomer_docking_difference.sh' % model_dir
                chemprop_input_arguments = {'function': 'bash',
                                            'arguments': [os.path.abspath(bash_file),
                                                          os.path.abspath(file_to_run),
                                                          os.path.abspath(output_filepath)]}
            elif constraining_model == 'biotoxicity':
                model_dir = r'.\chemprop_predicting\chemprop_models\%s' % model_details['model_version']
                model_checkpoint_dir = model_dir
                bash_file = r'%s\run_biotoxicity_model.sh' % model_dir
                chemprop_input_arguments = {'function': 'bash',
                                            'arguments': [os.path.abspath(bash_file),
                                                          os.path.abspath(file_to_run),
                                                          os.path.abspath(output_filepath)]}
            else:
                print('\t\t\tConstraining model: %s not valid' % constraining_model)
                unavailable_models[0] = True
                unavailable_models[1].append([constraining_model, model_version])
                continue
            if not os.path.exists(model_checkpoint_dir):
                print('\t\t\tConstraining model: %s is unavailable' % constraining_model)
                unavailable_models[0] = True
                unavailable_models[1].append([constraining_model, model_version])
                continue

            if chemprop_input_arguments['function'] == 'chemprop_predict':
                chemprop_arguments = chemprop_input_arguments['arguments']
                chemprop_input_command = 'chemprop_predict ' + ' '.join(['--%s %s' % (key, chemprop_arguments[key]) for
                                                                         key in chemprop_arguments.keys()])
            elif chemprop_input_arguments['function'] == 'bash':
                chemprop_arguments = chemprop_input_arguments['arguments']
                chemprop_input_command = 'bash ' + ' '.join(chemprop_arguments)
            else:
                unavailable_models[0] = True
                unavailable_models[1].append([constraining_model, model_version])
                continue

            chemprop_thread = threading.Thread(target=chemprop_thread_function,
                                               args=([chemprop_input_command, file_to_run, chemprop_queue]))
            chemprop_thread.start()
            chemprop_thread.join()
            chemprop_status[constraining_model]['output_files'].append(output_filepath)
            with open(chemprop_status_filepath, 'w') as yaml_file:
                yaml.dump(chemprop_status, yaml_file)
            new_chemprop_data += 1

    for constraining_model in chemprop_status.keys():
        model_details = chemprop_status[constraining_model]
        if len(list(set(model_details['input_files']).difference(set(model_details['completed_files'])))) != 0:
            new_chemprop_data += 1

    if new_chemprop_data == 0:
        return ['Success', 'All molecules have been predicted!', new_chemprop_data, unavailable_models]
    print('\tWorking on compiling the ChemProp predictions')
    compiled_chemprop_datasets = {}
    for constraining_model in chemprop_status.keys():
        print(constraining_model)
        chemprop_status[constraining_model]['completed_files'] = []
        for output_file in chemprop_status[constraining_model]['output_files']:
            with open(output_file, 'r') as infile:
                next(infile)
                for line in infile:
                    if len(line.strip('\n')) < 5:
                        continue
                    line_split = line.strip('\n').split(',')
                    if 'Invalid SMILES' in line_split:
                        continue
                    if constraining_model == 'uv_vis':
                        prediction = {'predicted_value': float(line_split[2]),
                                      'prediction_solvent': line_split[1],
                                      'prediction_uncertainty': float(np.sqrt(float(line_split[3])))}
                    else:
                        prediction = {'predicted_value': float(line_split[1]),
                                      'prediction_uncertainty': float(np.sqrt(float(line_split[2])))}
                    if line_split[0] not in compiled_chemprop_datasets.keys():
                        compiled_chemprop_datasets[line_split[0]] = {}
                    compiled_chemprop_datasets[line_split[0]][constraining_model] = prediction
            chemprop_status[constraining_model]['completed_files'].append(output_file)

    chemprop_compiled_output_filepath = r'%s\%s_compiled_chemprop_predictions.json' % (chemprop_model_output, job_name)
    with open(chemprop_compiled_output_filepath, 'w') as json_file:
        json.dump(compiled_chemprop_datasets, json_file)
    with open(chemprop_status_filepath, 'w') as yaml_file:
        yaml.dump(chemprop_status, yaml_file)

    return_statement, return_object = plot_chemprop_predictions(chemprop_model_output,
                                                                compiled_chemprop_datasets,
                                                                job_name,
                                                                chemprop_details['constraining_models'])
    if return_statement != 'Success':
        return ['Error', return_object, new_chemprop_data, unavailable_models]

    return ['Success', chemprop_compiled_output_filepath, new_chemprop_data, unavailable_models]


def chemprop_thread_function(chemprop_command, job_name, chemprop_queue):
    print('Received command: %s' % chemprop_command)
    current_chemprop_status = chemprop_queue.get()
    if job_name not in current_chemprop_status['running_jobs']:
        current_chemprop_status['running_jobs'].append(job_name)
    chemprop_queue.put(current_chemprop_status)
    process = subprocess.run(chemprop_command, shell=False, capture_output=True)
    formatted_stderr = []
    for line in process.stderr.decode().split('\r'):
        if 'WARNING' in line:
            continue
        clean_line = line.replace('\r', '').strip('\n')
        if clean_line.strip() == '':
            continue
        if '#' not in clean_line:
            continue
        formatted_stderr.append(clean_line)
    print(process.returncode)
    if process.returncode != 0:
        pprint(formatted_stderr)
    print('Finishing!: %s' % chemprop_command)
    current_chemprop_status = chemprop_queue.get()
    if job_name in current_chemprop_status['running_jobs']:
        current_chemprop_status['running_jobs'].remove(job_name)
    if process.returncode != 0:
        current_chemprop_status['error_state'].append(job_name)
    chemprop_queue.put(current_chemprop_status)


def plot_chemprop_predictions(chemprop_model_output, compiled_chemprop_datasets, job_name, constraining_models):
    model_data_plotting = {}
    for product_key in compiled_chemprop_datasets.keys():
        product_predictions = compiled_chemprop_datasets[product_key]
        if any(model not in product_predictions.keys() for model in constraining_models):
            # sys.exit('Something prevented Chemprop from making predictions... %s' % constraining_models)
            continue
        for model in constraining_models.keys():
            if model not in model_data_plotting.keys():
                model_data_plotting[model] = {'data': [], 'uncertainty': []}
                if type(constraining_models[model]['center']) == int:
                    model_data_plotting[model]['targets'] = [constraining_models[model]['center']]
                    model_data_plotting[model]['widths'] = [constraining_models[model]['range']]
                else:
                    model_data_plotting[model]['targets'] = constraining_models[model]['center']
                    model_data_plotting[model]['widths'] = constraining_models[model]['range']
            if model not in product_predictions.keys():
                continue
            model_predictions = product_predictions[model]
            model_data_plotting[model]['data'].append(model_predictions['predicted_value'])
            model_data_plotting[model]['uncertainty'].append(model_predictions['prediction_uncertainty'])
    model_combinations = list(itertools.permutations(list(constraining_models.keys()), 2))
    for model_combination in model_combinations:
        xdata = model_data_plotting[model_combination[0]]
        ydata = model_data_plotting[model_combination[1]]
        if 'uncertainty' in xdata.keys():
            x_error = xdata['uncertainty']
        else:
            x_error = None
        if 'uncertainty' in ydata.keys():
            y_error = ydata['uncertainty']
        else:
            y_error = None
        fig, ax = plt.subplots()
        ax.errorbar(xdata['data'], ydata['data'], xerr=x_error, yerr=y_error,
                    fmt='bo', ecolor='k', markeredgecolor='k', alpha=0.8)
        x_targets = model_data_plotting[model_combination[0]]['targets']
        y_targets = model_data_plotting[model_combination[1]]['targets']
        if len(x_targets) > 1:
            for x_target_index, x_target_value in enumerate(x_targets):
                for y_target_index, y_target_value in enumerate(y_targets):
                    x_width = model_data_plotting[model_combination[0]]['widths'][x_target_index]
                    y_width = model_data_plotting[model_combination[1]]['widths'][y_target_index]
                    anchor_point = (x_target_value - x_width, y_target_value - y_width)
                    ax.add_patch(Rectangle(anchor_point, 2 * x_width, 2 * y_width,
                                           edgecolor='black', facecolor=(1, 0.5, 0, 0.2), zorder=3))
        plt.xlabel(model_combination[0])
        plt.ylabel(model_combination[1])
        figure_name = '%s_%s_%s_figure.png' % (job_name, model_combination[0], model_combination[1])
        plt.savefig(r'%s\%s' % (chemprop_model_output, figure_name), dpi=300)

    return ['Success', '']
