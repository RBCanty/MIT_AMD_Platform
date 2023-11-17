#!/usr/bin/env python
from __future__ import print_function, division
import numpy as np
from rdkit import Chem
from rdkit import rdBase
from rdkit.Chem import AllChem
from rdkit import DataStructs
import rdkit.Chem.QED as QED
import rationale.scripts.sascorer as sascorer
import os
import pickle
import joblib

from chemprop.train import predict
from chemprop.data import MoleculeDataset, MoleculeDataLoader
from chemprop.data.utils import get_data, get_data_from_smiles
from chemprop.utils import load_args, load_checkpoint, load_scalers

from rationale.utils.tddft_predictions_inverse_transform_convert import get_wavelength

rdBase.DisableLog('rdApp.error')


class gsk3_model():
    """Scores based on an ECFP classifier for activity."""

    kwargs = ["clf_path"]
    clf_path = 'data/gsk3/gsk3.pkl'

    def __init__(self):
        with open(self.clf_path, "rb") as f:
            self.clf = pickle.load(f)

    def __call__(self, smiles_list):
        fps = []
        mask = []
        for i,smiles in enumerate(smiles_list):
            mol = Chem.MolFromSmiles(smiles)
            mask.append( int(mol is not None) )
            fp = gsk3_model.fingerprints_from_mol(mol) if mol else np.zeros((1, 2048))
            fps.append(fp)

        fps = np.concatenate(fps, axis=0)
        scores = self.clf.predict_proba(fps)[:, 1]
        scores = scores * np.array(mask)
        return np.float32(scores)

    @classmethod
    def fingerprints_from_mol(cls, mol):  # use ECFP4
        features_vec = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        features = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(features_vec, features)
        return features.reshape(1, -1)
    

class jnk3_model():
    """Scores based on an ECFP classifier for activity."""

    kwargs = ["clf_path"]
    clf_path = 'data/jnk3/jnk3.pkl'

    def __init__(self):
        with open(self.clf_path, "rb") as f:
            self.clf = pickle.load(f)

    def __call__(self, smiles_list):
        fps = []
        mask = []
        for i,smiles in enumerate(smiles_list):
            mol = Chem.MolFromSmiles(smiles)
            mask.append( int(mol is not None) )
            fp = jnk3_model.fingerprints_from_mol(mol) if mol else np.zeros((1, 2048))
            fps.append(fp)

        fps = np.concatenate(fps, axis=0)
        scores = self.clf.predict_proba(fps)[:, 1]
        scores = scores * np.array(mask)
        return np.float32(scores)

    @classmethod
    def fingerprints_from_mol(cls, mol):  # use ECFP4
        features_vec = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        features = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(features_vec, features)
        return features.reshape(1, -1)


class qed_func():

    def __call__(self, smiles_list):
        scores = []
        for smiles in smiles_list:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                scores.append(0)
            else:
                scores.append(QED.qed(mol))
        return np.float32(scores)


class sa_func():

    def __call__(self, smiles_list):
        scores = []
        for smiles in smiles_list:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                scores.append(100)
            else:
                scores.append(sascorer.calculateScore(mol))
        return np.float32(scores)
    
class tddft_func():    
    def __init__(self, checkpoint_dir='tddft_predictor'):
        self.checkpoints = []
        for root, _, files in os.walk(checkpoint_dir):
            for fname in files:
                if fname.endswith('.pt'):
                    fname = os.path.join(root, fname)
                    #self.scaler, self.features_scaler = load_scalers(fname)
                    scalers =load_scalers(fname)
                    self.scaler, self.features_scaler = scalers[0], scalers[1]
                    self.train_args = load_args(fname)
                    model = load_checkpoint(fname) #, cuda=True)
                    self.checkpoints.append(model)

    def __call__(self, smiles, batch_size=500):
        test_data = get_data_from_smiles(smiles=[[x] for x in smiles], skip_invalid_smiles=False) #, args=self.train_args)
        test_data_loader = MoleculeDataLoader(
            dataset=test_data,
            batch_size=batch_size,
            num_workers=8)
        
        valid_indices = [i for i in range(len(test_data)) if test_data[i].mol is not None]
        full_data = test_data
        test_data = MoleculeDataset([test_data[i] for i in valid_indices])

        if self.train_args.features_scaling:
            test_data.normalize_features(self.features_scaler)

        sum_preds = np.zeros((len(test_data), 1))
        for model in self.checkpoints:
            
            model_preds = predict(
                model=model,
                data_loader=test_data_loader,
                scaler=self.scaler
            )
            sum_preds += np.array(model_preds)

        # Ensemble predictions
        avg_preds = sum_preds / len(self.checkpoints)
        avg_preds = avg_preds.squeeze(-1).tolist()

        # Put zero for invalid smiles
        full_preds = [0.0] * len(full_data)
        for i, si in enumerate(valid_indices):
            full_preds[si] = avg_preds[i]
        vacuum_tddft_peak_scaler = joblib.load('rationale/utils/vacuum_tddft_peaks_std_scaler.bin')
        return get_wavelength(np.array(full_preds, dtype=np.float32),vacuum_tddft_peak_scaler)

class tddft_epistempic_unc_func():    
    # Function for manipulating the ensemble variance from a chemprop model
    def __init__(self, args):
        self.checkpoints = []
        for root, _, files in os.walk(args.predictor_dir):
            for fname in files:
                if fname.endswith('.pt'):
                    fname = os.path.join(root, fname)
                    scalers =load_scalers(fname)
                    self.scaler, self.features_scaler = scalers[0], scalers[1]
                    self.train_args = load_args(fname)
                    model = load_checkpoint(fname) #, cuda=True)
                    self.checkpoints.append(model)

    def __call__(self, smiles, batch_size=500):
        # Takes as input a list of smiles and outputs a list of ensemble variances
        test_data = get_data_from_smiles(smiles=[[x] for x in smiles], skip_invalid_smiles=True) #, args=self.train_args)
        test_data_loader = MoleculeDataLoader(
            dataset=test_data,
            batch_size=batch_size,
            num_workers=8)
        
        valid_indices = [i for i in range(len(test_data)) if test_data[i].mol is not None]
        full_data = test_data
        test_data = MoleculeDataset([test_data[i] for i in valid_indices])

        if self.train_args.features_scaling:
            test_data.normalize_features(self.features_scaler)

        all_preds = []
        for model in self.checkpoints:
            
            model_preds = predict(
                model=model,
                data_loader=test_data_loader,
                scaler=self.scaler
            )
            all_preds.append(np.array(model_preds))

        return np.var(all_preds,axis=0)    

# Not updated! Problems may occur related to multiple molecules:
class chemprop_model():
    
    def __init__(self, checkpoint_dir):
        self.checkpoints = []
        for root, _, files in os.walk(checkpoint_dir):
            for fname in files:
                if fname.endswith('.pt'):
                    fname = os.path.join(root, fname)
                    #self.scaler, self.features_scaler = load_scalers(fname)
                    scalers =load_scalers(fname)
                    self.scaler, self.features_scaler = scalers[0], scalers[1]
                    self.train_args = load_args(fname)
                    model = load_checkpoint(fname) #, cuda=True)
                    self.checkpoints.append(model)

    def __call__(self, smiles, batch_size=500):
        test_data = get_data_from_smiles(smiles=smiles, skip_invalid_smiles=False) #, args=self.train_args)
        test_data_loader = MoleculeDataLoader(
            dataset=test_data,
            batch_size=batch_size,
            num_workers=8)
        
        valid_indices = [i for i in range(len(test_data)) if test_data[i].mol is not None]
        full_data = test_data
        test_data = MoleculeDataset([test_data[i] for i in valid_indices])

        if self.train_args.features_scaling:
            test_data.normalize_features(self.features_scaler)

        sum_preds = np.zeros((len(test_data), 1))
        for model in self.checkpoints:
            
            model_preds = predict(
                model=model,
                data_loader=test_data_loader,
                scaler=self.scaler
            )
            sum_preds += np.array(model_preds)

        # Ensemble predictions
        avg_preds = sum_preds / len(self.checkpoints)
        avg_preds = avg_preds.squeeze(-1).tolist()

        # Put zero for invalid smiles
        full_preds = [0.0] * len(full_data)
        for i, si in enumerate(valid_indices):
            full_preds[si] = avg_preds[i]

        return np.array(full_preds, dtype=np.float32)


def get_scoring_function(prop_name,args):
    """Function that initializes and returns a scoring function by name"""
    if prop_name == 'jnk3':
        return jnk3_model()
    elif prop_name == 'gsk3':
        return gsk3_model()
    elif prop_name == 'qed':
        return qed_func()
    elif prop_name == 'sa':
        return sa_func()
    elif prop_name =='tddft':
        return tddft_func()
    elif prop_name =='tddft_ensemble':
        return tddft_epistempic_unc_func(args)
    else:
        return chemprop_model(prop_name)

if __name__ == "__main__":
    import sys
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('--prop', required=True)

    args = parser.parse_args()
    funcs = [get_scoring_function(prop) for prop in args.prop.split(',')]

    data = [line.split()[:2] for line in sys.stdin]
    all_x, all_y = zip(*data)
    props = [func(all_y) for func in funcs]

    col_list = [all_x, all_y] + props
    for tup in zip(*col_list):
        print(*tup)
