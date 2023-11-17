import numpy as np
import pandas as pd
import rdkit
from rationale.generate.rdfilters import filter_mols

def filter_parse(infile,outfile):
    data = pd.read_csv(infile,header=None,sep=' ')
    data['index'] = [x for x in range(len(data))]
    data['mol_num'] = data['index'].apply(lambda x: "MOL{:04d}".format(x))
    data[[1,'mol_num']].to_csv(outfile,index=False,header=None,sep=' ')

def run_filter(args):
    in_file = args.out_smiles_list
    middle_file = args.save_dir+'/gen_list_clean.csv'
    out_file = args.out_smiles_list.replace('.txt', '_filtered.txt')
    filter_parse(in_file,middle_file)
    filter_mols(middle_file,out_file,
                alert_file_name='rationale/utils/rdfilters_settings/alert_collection.csv',
                rules_file_name= '/home/gridsan/cbilod/rd_filters/rd_filters/data/rules_main.json')#'rationale/utils/rdfilters_settings/rules.json')
    
    
    
    
    