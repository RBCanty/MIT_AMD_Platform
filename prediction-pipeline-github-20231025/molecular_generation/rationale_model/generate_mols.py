# -*- coding: utf-8 -*-
"""
Adapted from Camille's implementation of the Rationale Model
"""

import sys
import os
if os.path.abspath(r'.\molecular_generation\rationale_model') not in sys.path:
    sys.path.append(os.path.abspath(r'.\molecular_generation\rationale_model'))

from rationale.generate_mols import generate_mols
import traceback

from tap import Tap  # pip install typed-argument-parser (https://github.com/swansonk14/typed-argument-parser)
from rationale.fuseprop.vocab import common_atom_vocab

def run_rationale_model(input_file, output_file, num_decodes, pretrained_model = ''):
    class CommonArgs(Tap):
    
        # Files (can be changed)
        save_dir: str = './molecular_generation/rationale_model/rationale/utils/'
        """
        Location to store standard output files. These include files for model training (formatted model.x),
        the generated list (gen_list).csv, and the filtered list (filtered_list.smi). 
        """
        running_smiles_list: str = None
        """
        Path to pre-existing list of molecules for consideration during selection. If this is None, the model will consider the list empty.
        """
        out_smiles_list: str = output_file
        """
        Full list of molecules generated in this iteration or included in running_smiles_list.
        """
        selected_mols_list: str = 'test/test_selected_list.csv'
        """
        List of molecules selected for next round of active learning.
        """
        unselected_mols_list: str = 'test/test_unselected_list.csv'
        """
        Lit of molecules not selected for next round of active learning. This list should be used for running_smiles_list for the next iteration.
        """
        predictor_dir: str = 'chemprop_models/tddft_ensemble'
        """
        Path to ChemProp folder containing model folders that predicts the results of tddft simulations (must use an ensemble).
        """   
        init_model: str = './molecular_generation/rationale_model/rationale/utils/pretrained_model.19'
        """
        (Doesn't need to be changed) Starting model to be used for finetuning. Defaults to model trained on ChEMBl and Dye datasets.
        """

        rationale: str = input_file
        """
        (Doesn't need to be changed) Location of rationale file (currently at rationale/utils/scaffold.txt, but this can be changed for other applications)
        """
        
        # Settings to Change:
        num_decode: int = int(num_decodes)
        """
        Number of times to decode each scaffold (Increasing this number will cause the model to try to generate a larger number of molecules).
        """
        epoch: int = 20
        """
        Number of epochs to train the finetuned model for.
        """
        cutoff: float = 200.0
        """
        Threshold ensemble variance considered for finetuning. Decreasing/increasing this will reduce/increase the constraint on the model resulting
        in increased/decreasd diversity.
        """
        active_learning_batch_size: int = 30
        """
        Number of molecules to be selected for active learning iteration (based on top N epistemic uncertainty).
        """

        # Settings Not to Change:
        atom_vocab = common_atom_vocab
        prop: str = 'tddft_ensemble'
        load_epoch: int = -1
        constraint_func: str = 'tddft_func'
        decode_batch_size: int = 20
            
        # Obsolete:
        model: str = None

        # Model Training Parameters:
        seed: int = 1
        rnn_type: str = 'LSTM'
        hidden_size: int = 400
        embed_size: int = 400
        batch_size: int = 20
        latent_size: int = 20
        depth: int = 10
        diter: int = 3
        lr: float = 5e-4
        clip_norm: float = 20.0
        alpha: float = 1.0
        anneal_rate: float = 1.0
        print_iter: int = 50
        beta: float = 0.3
        
        def __init__(self, *args, **kwargs):
            super(CommonArgs, self).__init__(*args, **kwargs)
    
    
    try:
        common_args = CommonArgs()
        args = common_args.parse_args()
        generate_mols(args)
        return ['Success', '']
    except:
        return ['Error', traceback.format_exc()]

