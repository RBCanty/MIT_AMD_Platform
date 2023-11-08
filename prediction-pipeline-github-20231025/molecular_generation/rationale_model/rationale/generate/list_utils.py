import pandas as pd
from rationale.generate.properties import tddft_epistempic_unc_func

def join_mol_list(args):
    generated_mols = pd.read_csv(args.save_dir+'/filtered_list.smi',header=None,sep=' ')[0].values
    if args.running_smiles_list is not None:
        prev_mols = pd.read_csv(args.running_smiles_list,header=None,sep=' ')[0].values
    else:
        prev_mols = []
    all_mols = list(set(list(generated_mols)+list(prev_mols)))
    pd.DataFrame(all_mols).to_csv(args.out_smiles_list,index=False,header=None)
    
def select_mols(args):
    print('Selecting Molecules for Active Learning Loop of Batch Size {}'.format(args.active_learning_batch_size))
    all_mols = pd.read_csv(args.out_smiles_list,header=None,sep=' ')[0].values
    func = tddft_epistempic_unc_func(args)
    variances = func(all_mols)
    
    df = pd.DataFrame()
    df['SMILES'] = all_mols
    df['variance'] = variances
#     print(df.sort_values('variance',ascending=False).iloc[0:10])
    selected_mols = df.sort_values('variance',ascending=False).iloc[0:args.active_learning_batch_size][['SMILES']]
    unselected_mols = df.sort_values('variance',ascending=False).iloc[args.active_learning_batch_size:][['SMILES']]
    selected_mols.to_csv(args.selected_mols_list,index=False,header=None)
    unselected_mols.to_csv(args.unselected_mols_list,index=False,header=None)