import numpy as np
import argparse
import sys

#-----------------------------------------------------------------------------------------------------------------------

from hmmix_additional_functions_v082 import *

# skov helper functions from git clone https://github.com/LauritsSkov/Introgression-detection
skov_dir = "resources/skov_scripts/Introgression-detection/src/"
sys.path.append(skov_dir)

from make_mutationrate import make_mutation_rate
#from make_test_data import create_test_data
#from hmm_functions import TrainModel, DecodeModel, HMMParam, read_HMM_parameters_from_file, get_default_HMM_parameters, write_HMM_to_file

from make_test_data import simulate_path, write_data
from hmm_functions import TrainModel, HMMParam, get_default_HMM_parameters, write_HMM_to_file, read_HMM_parameters_from_file, Write_Decoded_output, Calculate_Posterior_probabillities, Emission_probs_poisson

from helper_functions import Load_observations_weights_mutrates, handle_individuals_input, handle_infiles, combined_files
from bcf_vcf import make_out_group, make_ingroup_obs
#-----------------------------------------------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description="input params from smk")
parser.add_argument("--inref_list", type=str, help="input.ref_list")
parser.add_argument("--intgt_list", type=str, help="input.tgt_list")
parser.add_argument("--in_vcf", type=str, help="input.vcf")
parser.add_argument("--in_true_tracts", type=str, help="input.true_tracts")

parser.add_argument("--skov_output_dir", type=str, help="skov_output_dir")
parser.add_argument("--seed", type=int, help="wildcards.seed")
parser.add_argument("--out_outgroup_file", type=str, help="output.outgroup_file")
parser.add_argument("--out_prob_file", type=str, help="output.prob_file")

parser.add_argument("--ref_set", type=list, help="ref_set")
parser.add_argument("--ref_id", type=str, help="ref_id")
parser.add_argument("--src_id", type=str, help="src_id")

parser.add_argument("--output_prefix", type=str, help="output_prefix")
parser.add_argument("--cutoff_num", type=int, help="cutoff_num")
parser.add_argument("--demog_id", type=str, help="demog_id")
parser.add_argument("--nref", type=int, help="nref")
parser.add_argument("--ntgt", type=int, help="ntgt")

args = parser.parse_args()

#-----------------------------------------------------------------------------------------------------------------------
# inputs: vcf, ref_list, tgt_list, true_tracts
inref_list = args.inref_list
intgt_list = args.intgt_list
in_vcf = args.in_vcf
in_true_tracts = args.in_true_tracts

#-----------------------------------------------------------------------------------------------------------------------
skov_output_dir = args.skov_output_dir
seed = args.seed

ref_set = args.ref_set
ref_id = args.ref_id
src_id = args.src_id
demog_id = args.demog_id
#-----------------------------------------------------------------------------------------------------------------------
# outputs : outgroup_file, prob_file
out_outgroup_file = args.out_outgroup_file
out_prob_file = args.out_prob_file
output_prefix = args.output_prefix
nref = args.nref
ntgt = args.ntgt
#-----------------------------------------------------------------------------------------------------------------------
cutoff_num = args.cutoff_num
cutoff_list = np.round(np.linspace(0, 1, cutoff_num, endpoint=False), 2)
cutoff_list = np.append(cutoff_list, [0.99, 0.999])
#-----------------------------------------------------------------------------------------------------------------------

with open(inref_list, 'r') as file:
    ref_lines = file.readlines()

with open(intgt_list, 'r') as file:
    tgt_lines = file.readlines()

ref_lines_comma = ','.join([line.strip() for line in ref_lines])

out_folder = os.path.join(skov_output_dir, f"{seed}")

with open(os.path.join(out_folder,'comma_separated_refinds.txt'), 'w') as output_file:
    output_file.write(ref_lines_comma)

comma_separated = [line.strip() for line in ref_lines]
comma_separated_tgt = [line.strip() for line in tgt_lines]

make_out_group_custom_mut(comma_separated, None, [in_vcf], out_outgroup_file, [None], [None], ref_set=ref_set)
make_ingroup_custom_mut(comma_separated_tgt, None, [in_vcf], os.path.join(out_folder , "output_tgt_vs"),  out_outgroup_file, [None], ref_set=ref_set)

make_mutation_rate(out_outgroup_file , os.path.join(out_folder, "mutrates_outgroup1.out"), None, 100000)

new_HMM = HMMParam([ref_id, src_id], [0.5, 0.5], [[0.99,0.01],[0.02,0.98]], [0.03, 0.3])
write_HMM_to_file(new_HMM, os.path.join(out_folder,'hmm_guesses.json'))

print(new_HMM) # check has got here

infiles, trained_files = train_hmm_individuals(comma_separated_tgt,  os.path.join(out_folder,'hmm_guesses.json'), os.path.join(out_folder, "mutrates_outgroup1.out"), out_folder=out_folder, window_size=1000, haploid=False, weights=None)

all_segments = decode_hmm_individuals(comma_separated_tgt, trained_files, os.path.join(out_folder, "mutrates_outgroup1.out"),  out_folder=out_folder, window_size=1000, haploid=False, weights=None)

print(all_segments) # check

all_segments.to_csv(out_prob_file, index=False)

inferred_tract_files = process_output(all_segments, out_folder, output_prefix, src_id, cutoff_list=cutoff_list, return_filenames=True)

print(inferred_tract_files) # check

for inferred_tracts in inferred_tract_files:
    precision, recall = cal_accuracy(in_true_tracts, inferred_tracts[1])
    cutoff = inferred_tracts[0]
    prec_rec_file = inferred_tracts[1].rsplit('.', 1)[0] + ".accuracy"
    with open(prec_rec_file, 'w') as o:
        o.write(f'{demog_id}\tnref_{nref}_ntgt_{ntgt}\t{cutoff}\t{precision}\t{recall}\n')
