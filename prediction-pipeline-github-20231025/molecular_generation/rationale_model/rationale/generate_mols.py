from rationale.generate.finetune import finetune
from rationale.generate.decode import decode
from rationale.generate.filter import run_filter
from rationale.generate.list_utils import join_mol_list, select_mols

def generate_mols(args):
    #finetune(args)  # Train fine-tuning model
    decode(args) # Decode molecules
    #run_filter(args) # Filter molecules (using rdfilter, settings stored in utils)
    #join_mol_list(args) # Join generated molecule list with previous list
    #select_mols(args) # Predict uncertainty and select molecules with highest uncertainty
    