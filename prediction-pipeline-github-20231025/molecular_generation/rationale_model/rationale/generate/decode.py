import torch
import torch.nn as nn
from torch.utils.data import DataLoader

import math, random, sys
import numpy as np
import argparse
import rdkit
from tqdm import tqdm

from rationale.fuseprop import *

def decode(args):
    random.seed(1)
    model = AtomVGNN(args).cuda()
    model_ckpt = torch.load(args.save_dir+'/pretrained_model.'+str(args.epoch-1))
#     if type(model_ckpt) is tuple:
#         print('loading model with rationale distribution', file=sys.stderr)
#         testdata = list(model_ckpt[0].keys())
#         model.load_state_dict(model_ckpt[1])
#     else:
#         print('loading pre-trained model', file=sys.stderr)
    testdata = [line.split()[1] for line in open(args.rationale)] 
    testdata = unique_rationales(testdata)
    model.load_state_dict(model_ckpt)

    print('total # rationales:', len(testdata), file=sys.stderr)
    model.eval()
    dataset = SubgraphDataset(testdata, args.atom_vocab, args.batch_size, args.num_decode)

    loader = DataLoader(dataset, batch_size=1, shuffle=False, num_workers=0, collate_fn=lambda x:x[0])
    torch.manual_seed(args.seed)
    torch.cuda.manual_seed(args.seed)

    with torch.no_grad():
        with open(args.out_smiles_list, 'w') as f:
            for init_smiles in tqdm(loader, ncols=50):
                final_smiles = model.decode(init_smiles)
                for x,y in zip(init_smiles, final_smiles):
                    f.write('{}, {}\n'.format(x, y))

