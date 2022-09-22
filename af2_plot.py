# script adapted from https://github.com/jasperzuallaert/VIBFold/blob/main/visualize_alphafold_results.py
# generates plddt and pae plots from AF2 predictions. only plddt for monomers

import os
import numpy as np
from matplotlib import pyplot as plt
import argparse
import pickle

def get_pae_plddt(model_names, is_multimer):
    out = {}
    for i,name in enumerate(model_names):
        d = pickle.load(open(name,'rb'))
        #print(d)
        out[f'model_{i+1}'] = {'plddt': d['plddt'], 'pae':(d['predicted_aligned_error'] if is_multimer == True else False)}

    return out

def generate_output_images(feature_dict, out_dir, name, pae_plddt_per_model,is_multimer):
    msa = feature_dict['msa']
    seqid = (np.array(msa[0] == msa).mean(-1))
    seqid_sort = seqid.argsort()
    non_gaps = (msa != 21).astype(float)
    non_gaps[non_gaps == 0] = np.nan
    final = non_gaps[seqid_sort] * seqid[seqid_sort, None]

    ##################################################################
    # builds plddt and sequence coverage plots
    ##################################################################
    plt.figure(figsize=(28, 8), dpi=400)
    plt.subplot(1, 2, 1)
    plt.title("Sequence coverage")
    plt.imshow(final,
               interpolation='nearest', aspect='auto',
               cmap="rainbow_r", vmin=0, vmax=1, origin='lower')
    plt.plot((msa != 21).sum(0), color='black')
    plt.xlim(-0.5, msa.shape[1] - 0.5)
    plt.ylim(-0.5, msa.shape[0] - 0.5)
    plt.colorbar(label="Sequence identity to query", )
    plt.xlabel("Positions")
    plt.ylabel("Sequences")

    ##################################################################
    plt.subplot(1, 2, 2)
    plt.title("Predicted LDDT per position")
    for model_name, value in pae_plddt_per_model.items():
        plt.plot(value["plddt"], label=model_name)
    plt.legend()
    plt.ylim(0, 100)
    plt.ylabel("Predicted LDDT")
    plt.xlabel("Positions")
    plt.savefig(f"{out_dir}/{name+('_' if name else '')}coverage_LDDT.png")
    ##################################################################
    # builds PAE plot
    ##################################################################
    if is_multimer == True: 
        num_models = 25

        plt.figure(figsize=(5 * num_models, 5), dpi=200)
        for n, (model_name, value) in enumerate(pae_plddt_per_model.items()):
            plt.subplot(1, num_models, n + 1)
            plt.title(model_name)
            plt.imshow(value["pae"], label=model_name, cmap="bwr", vmin=0, vmax=30)
            plt.colorbar()
        plt.savefig(f"{out_dir}/{name+('_' if name else '')}PAE.png")
    else:
        num_models = 5

    ##################################################################

# reads data from pkl files, determines if multimer, then generates list of file names
def main():
    # print(os.listdir(args.input_dir))
    feature_dict = pickle.load(open(f'{args.input_dir}/features.pkl','rb'))
    is_multimer = ('result_model_1_multimer_v2_pred_0.pkl' in [os.path.basename(f) for f in os.listdir(path=args.input_dir)])
    #is_ptm = ('result_model_1_ptm.pkl' in [os.path.basename(f) for f in os.listdir(path=args.input_dir)])
    if is_multimer == True:
        lst1 = [1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5]
        lst2 = [0,0,0,0,0,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4]
        model_names = [f'{args.input_dir}/result_model_{f}_multimer_v2_pred_{g}.pkl' for f,g in zip(lst1,lst2)]
    else: 
        model_names = [f'{args.input_dir}/result_model_{f}_pred_0.pkl' for f in range(1,6)]  
    #print(model_names)

    pae_plddt_per_model = get_pae_plddt(model_names,is_multimer)
    generate_output_images(feature_dict, args.output_dir if args.output_dir else args.input_dir, args.name, pae_plddt_per_model, is_multimer)

# user inputs (ex:  python af2_plot.py --input_dir test/test01 --output_dir . --name test01 )
parser = argparse.ArgumentParser()
parser.add_argument('--input_dir',dest='input_dir',required=True)
parser.add_argument('--name',dest='name')
parser.set_defaults(name='')
parser.add_argument('--output_dir',dest='output_dir')
parser.set_defaults(output_dir='')
args = parser.parse_args()

#invoke main fxn 
main()


