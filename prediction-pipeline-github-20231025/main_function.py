
# Normal Module Imports
import yaml
import glob
import sys
import os

# Function imports
from functions.supporting_functions import directory_check
from molecular_generation.molecular_generation_functions import generate_molecules
from reaction_planning.reaction_pathway_planning import plan_reaction_pathways
from reaction_planning.pathway_expansion import reaction_pathway_expansion
from reaction_grouping.reaction_filtering import reaction_filtering
from chemprop_predicting.chemprop_functions import run_chemprop_models
from reaction_grouping.product_grouping import group_products
from reaction_grouping.groups_to_queues import reaction_plate_grouping
from reaction_grouping.groups_to_queues import build_queue_documents


input_files = glob.glob(r'.\input_jobs\*.yaml')
output_files = glob.glob(r'.\output_jobs\*.yaml')
for input_file in input_files:
    # Here we are going to check each of the input files to see what jobs still need to be executed
    job_name = input_file.split('\\')[-1].split('.yaml')[0]
    print(job_name)
    print(input_file.split('\\')[-1].strip('.yaml'))
    folder_output = r'.\output_jobs\%s' % job_name
    print('Working on job: %s' % job_name)
    return_statement, return_object = directory_check(folder_output)
    if return_statement != 'Success':
        sys.exit('\tError encountered: %s' % return_object)
    else:
        print('\t%s' % return_object)

    # Start by loading the input file
    # The input file contains options that change the way predictions and planning are done in the prediction
    # pipeline, these can be adjusted. Refer to the input_job example for details on the input formatting
    print('\tOpening the input file: %s' % input_file)
    with open(input_file, 'r') as yaml_file:
        job_information = yaml.load(yaml_file, Loader=yaml.Loader)

    # Check to see if the job has been completed
    if not job_information['job_status']['in_progress']:
        continue

    # Construct a job-status file to keep track of SMILES
    smiles_status_filepath = r'%s\%s_smiles_set.yaml' % (folder_output, job_name)
    if os.path.exists(smiles_status_filepath):
        with open(smiles_status_filepath, 'r') as yaml_file:
            smiles_status = yaml.load(yaml_file, Loader=yaml.Loader)
    else:
        smiles_status = {'smiles_set': set()}

    # Check to see if there are extra SMILES in the job information document
    if 'smiles' in job_information.keys():
        smiles_status['smiles_set'].update(job_information['smiles'])
    # Then move to molecular generation
    if 'molecular_generation' in job_information.keys():
        return_statement, return_object = generate_molecules(job_information, job_name)
        if return_statement != 'Success':
            sys.exit('\tError encountered: %s' % return_object)
        # if return_object != '':
        mol_gen_filepath = r'.\output_jobs\%s\molecular_generation\%s_gen_smiles.txt' % (job_name, job_name)
        with open(mol_gen_filepath, 'r') as infile:
            for line in infile:
                smiles_status['smiles_set'].add(line.strip())
    # Now re-write the SMILES status file
    with open(smiles_status_filepath, 'w') as yaml_file:
        yaml.dump(smiles_status, yaml_file)

    # There is the option to just run Chemprop models if not using reaction planning, this is not the standard
    # operation of the prediction pipeline but exists for testing and special uses
    if 'reaction_planning' in job_information.keys():
        return_statement, return_object, askcos_reaction_plans = plan_reaction_pathways(job_information, job_name,
                                                                                        smiles_status_filepath)
        if return_statement != 'Success':
            print(return_object)
            sys.exit()
        reactions_filepath = r'%s\%s_compiled_askcos.json' % (folder_output, job_name)
        expansion_reaction_plans = False
        if 'pathway_expansion' in job_information.keys():
            expansion_details = job_information['pathway_expansion']
            if 'pathway_expansion' in expansion_details and expansion_details['pathway_expansion']:
                return_statement, return_object, expansion_reaction_plans = reaction_pathway_expansion(job_information,
                                                                                                       job_name)
                if return_statement != 'Success':
                    print(return_object)
                    sys.exit()
                reactions_filepath = r'%s\reaction_planning\%s_compiled_products.json' % (folder_output, job_name)

        if askcos_reaction_plans or expansion_reaction_plans:
            new_reaction_plans = True
        else:
            new_reaction_plans = False
        return_statement, return_object, filtered_reactions_filepath = reaction_filtering(job_information,
                                                                                          job_name,
                                                                                          reactions_filepath,
                                                                                          new_reaction_plans)
        if return_statement != 'Success':
            print(return_object)
            sys.exit()
        candidate_type = 'reactions'
        candidate_filepath = filtered_reactions_filepath
    else:
        candidate_type = 'smiles'
        candidate_filepath = smiles_status_filepath

    return_statement, chemprop_data_path, new_chemprop_data, unavailable_models = run_chemprop_models(job_information,
                                                                                                      job_name,
                                                                                                      candidate_filepath,
                                                                                                      candidate_type)
    if return_statement != 'Success':
        print(chemprop_data_path)
        sys.exit()

    # When the prediction pipeline attempts to run chemprop models that are not currently available (such as a
    # future version of the model that has not been trained yet, unavailable models are returned
    if unavailable_models[0]:
        print('Missing models: %s' % unavailable_models)
        continue




    # With the Chemprop Predictions and ASKCOS reaction pathways we can select and group candidates
    if 'reaction_planning' in job_information.keys():
        return_statement, temperature_sequence, candidate_reactions, molecule_values = group_products(job_information,
                                                                                                      job_name,
                                                                                                      chemprop_data_path,
                                                                                                      candidate_filepath)
        if return_statement != 'Success':
            print(return_object)
            sys.exit()
        return_statement, reaction_grouping = reaction_plate_grouping(job_information,
                                                                      job_name,
                                                                      temperature_sequence,
                                                                      candidate_reactions,
                                                                      molecule_values=molecule_values)
        if return_statement != 'Success':
            print(reaction_grouping)
            sys.exit()
        return_statement, temp_object = build_queue_documents(job_information,
                                                              job_name,
                                                              reaction_grouping,
                                                              folder_output)
        if return_statement != 'Success':
            print(temp_object)
            sys.exit()

