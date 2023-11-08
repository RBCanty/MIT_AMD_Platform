import torch
import torch.nn as nn
import torch.optim as optim
import torch.optim.lr_scheduler as lr_scheduler
from torch.utils.data import DataLoader

import argparse
import rdkit
import math, random, sys, os
import numpy as np
from tqdm import tqdm

from rationale.fuseprop import *
from .properties import get_scoring_function
from .constraint_function import *

def remove_order(s):
    for x in range(15):
        s = s.replace(":%d]" % (x,), ":1]")
    return s

# Decode molecules
def decode_rationales(model, rationale_dataset):
    loader = DataLoader(rationale_dataset, batch_size=1, shuffle=False, num_workers=0, collate_fn=lambda x:x[0])
    model.eval()
    cand_mols = []
    with torch.no_grad():
        for init_smiles in tqdm(loader):
            final_smiles = model.decode(init_smiles)
            mols = [(x,y) for x,y in zip(init_smiles, final_smiles) if y and '.' not in y]
            mols = [(x,y) for x,y in mols if Chem.MolFromSmiles(y).HasSubstructMatch(Chem.MolFromSmiles(x))]
            cand_mols.extend(mols)
    return cand_mols

# Predict properties and filter 
def property_filter(cand_mols, scoring_function, epoch, args):
    rationales, smiles_list = zip(*cand_mols)
    cand_props = scoring_function(smiles_list)
    new_data = []
    rationale_dist = {remove_order(r) : [0,0] for r in rationales}
    with open(args.save_dir + '/valid.' + str(epoch), 'w') as f:
        for (init_smiles, final_smiles), prop in zip(cand_mols, cand_props):
            rationale = remove_order(init_smiles)
            rationale_dist[rationale][1] += 1
            print(init_smiles, final_smiles, prop, file=f)
            if args.compare_func(prop):
                rationale_dist[rationale][0] += 1
                new_data.append( (init_smiles, final_smiles) )

    print('property filter: %d -> %d' % (len(cand_mols), len(new_data)))
    rationale_dist = {r : x / n for r,(x,n) in rationale_dist.items() if x / n >= args.alpha}
    return new_data, rationale_dist

def finetune(args):
    
    if args.prop == "tddft_ensemble":
        args.compare_func = tddft_ensemble_func(args.cutoff)
        print("Cutoff is: ",args.cutoff)
    else:
        raise Exception('Model does not currently support options other than tddft_ensemble')
        
    prop_funcs = [get_scoring_function(prop,args) for prop in args.prop.split(',')] #Right now this has to be tddft_ensemble
    scoring_function = lambda x : list( zip(*[func(x) for func in prop_funcs]) )

    with open(args.rationale) as f:
        rationales = [line.split()[1] for line in f]
        rationales = unique_rationales(rationales)
        rationale_dataset = SubgraphDataset(rationales, args.atom_vocab, args.decode_batch_size, args.num_decode)

    model = AtomVGNN(args).cuda()
    if args.load_epoch >= 0:
        path = os.path.join(args.save_dir, f"model.{args.load_epoch}")
        model.load_state_dict(torch.load(path)[1])
    else:
        model.load_state_dict(torch.load(args.init_model))

    print("Model #Params: %dK" % (sum([x.nelement() for x in model.parameters()]) / 1000,))

    optimizer = optim.Adam(model.parameters(), lr=args.lr)
    scheduler = lr_scheduler.ExponentialLR(optimizer, args.anneal_rate)

    param_norm = lambda m: math.sqrt(sum([p.norm().item() ** 2 for p in m.parameters()]))
    grad_norm = lambda m: math.sqrt(sum([p.grad.norm().item() ** 2 for p in m.parameters() if p.grad is not None]))

    for epoch in range(args.load_epoch + 1, args.epoch):
        print('epoch', epoch)
        cand_mols = decode_rationales(model, rationale_dataset)
        cand_mols, rationale_dist = property_filter(cand_mols, scoring_function, epoch, args)

        model_ckpt = (rationale_dist, model.state_dict())
        torch.save(model_ckpt, os.path.join(args.save_dir, f"model.{epoch}"))

        cand_mols = list(set(cand_mols))
        random.shuffle(cand_mols)

        # Update model
        dataset = MoleculeDataset(cand_mols, args.atom_vocab, args.batch_size)
        dataloader = DataLoader(dataset, batch_size=1, shuffle=True, num_workers=0, collate_fn=lambda x:x[0])
        model.train()

        meters = np.zeros(5)
        for total_step, batch in enumerate(dataset):
            if batch is None: continue

            model.zero_grad()
            loss, kl_div, wacc, tacc, sacc = model(*batch, beta=args.beta)
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), args.clip_norm)
            optimizer.step()

            meters = meters + np.array([kl_div, loss.item(), wacc * 100, tacc * 100, sacc * 100])

            if (total_step + 1) % args.print_iter == 0:
                meters /= args.print_iter
                print("[%d] Beta: %.3f, KL: %.2f, loss: %.3f, Word: %.2f, Topo: %.2f, Assm: %.2f, PNorm: %.2f, GNorm: %.2f" % (total_step + 1, args.beta, meters[0], meters[1], meters[2], meters[3], meters[4], param_norm(model), grad_norm(model)))
                sys.stdout.flush()
                meters *= 0
            
        scheduler.step()
        print("learning rate: %.6f" % scheduler.get_lr()[0])

