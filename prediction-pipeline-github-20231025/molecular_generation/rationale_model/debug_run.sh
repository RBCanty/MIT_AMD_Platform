#!/bin/bash

#SBATCH -o rationale_dev-%j-%a
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --gres=gpu:volta:1
#SBATCH --exclusive


init_model=ckpt/dyes_chembl_pretrained/model.19
save_dir=ckpt/debug
# save_dir=ckpt/dyes_chembl_pretrained_sa
rationale_path=data/amd/scaffold1.txt
num_decode=200
epochs=1
alpha=0.5



python finetune.py --init_model $init_model --save_dir $save_dir --rationale $rationale_path --num_decode $num_decode --prop tddft_ensemble --epoch $epochs --alpha $alpha

# python finetune.py --init_model $init_model --save_dir $save_dir --rationale $rationale_path --num_decode $num_decode --prop sa,/home/gridsan/cbilod/g2g_optimization/predictors/amd/20201127_joung_main_noshare/fold_0 --epoch $epochs --alpha $alpha

# python finetune.py --init_model $init_model --save_dir $save_dir --rationale $rationale_path --num_decode $num_decode --prop sa --epoch $epochs --alpha $alpha