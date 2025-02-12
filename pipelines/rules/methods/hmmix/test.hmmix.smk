import numpy as np
import os

import make_mutationrate
from make_mutationrate import make_mutation_rate
from make_test_data import create_test_data
from hmm_functions import TrainModel, DecodeModel, HMMParam, read_HMM_parameters_from_file
from helper_functions import Load_observations_weights_mutrates
from helper_functions import Load_observations_weights_mutrates, handle_individuals_input, handle_infiles, combined_files
from bcf_vcf import make_out_group, make_ingroup_obs
from hmm_functions import TrainModel, DecodeModel, HMMParam, read_HMM_parameters_from_file
from hmm_functions import HMMParam, get_default_HMM_parameters, write_HMM_to_file
from hmmix_additional_functions import *


rule process_test_data:
    input:
        vcf = output_dir + "/{seed}/" + output_prefix + ".vcf.gz",
        ref_list = output_dir + "/{seed}/" + output_prefix + ".ref.ind.list",
        tgt_list = output_dir + "/{seed}/" + output_prefix + ".tgt.ind.list",
        true_tracts = output_dir + "/{seed}/" + output_prefix + ".truth.tracts.bed",
    output:
        outgroup_file = skov_output_dir + "/{seed}/output_ref_vs.txt",
        prob_file = skov_output_dir + "/{seed}/probabilities.txt",
    resources:
        partition = "basic",
        mem_gb = 32,
        cpus = 16,
    log:
        "logs/hmmix/proctestdata.{seed}.log",
    benchmark:
        "benchmarks/hmmix/proctestdata.{seed}.benchmark.txt",
    run:
        with open(input.ref_list, 'r') as file:

            ref_lines = file.readlines()

        with open(input.tgt_list, 'r') as file:

            tgt_lines = file.readlines()

        ref_lines_comma = ','.join([line.strip() for line in ref_lines])

        out_folder = os.path.join(skov_output_dir, wildcards.seed)
        with open(os.path.join(out_folder,'comma_separated_refinds.txt'), 'w') as output_file:
            output_file.write(ref_lines_comma)

        comma_separated = [line.strip() for line in ref_lines]
        comma_separated_tgt = [line.strip() for line in tgt_lines]

        make_out_group_custom_mut(comma_separated, None, [input.vcf], output.outgroup_file, [None], [None], ref_set=ref_set)
        make_ingroup_custom_mut(comma_separated_tgt, None, [input.vcf], os.path.join(out_folder , "output_tgt_vs"),  output.outgroup_file, [None], ref_set=ref_set)

        make_mutation_rate(output.outgroup_file , os.path.join(out_folder, "mutrates_outgroup1.out"), None, 100000)
        new_HMM = HMMParam([ref_id, src_id], [0.5, 0.5], [[0.99,0.01],[0.02,0.98]], [0.03, 0.3])
        write_HMM_to_file(new_HMM, os.path.join(out_folder,'hmm_guesses.json'))

        infiles, trained_files = train_hmm_individuals(comma_separated_tgt,  os.path.join(out_folder,'hmm_guesses.json'), os.path.join(out_folder, "mutrates_outgroup1.out"), out_folder=out_folder, window_size=1000, haploid=False, weights=None)

        all_segments = decode_hmm_individuals(comma_separated_tgt, trained_files, os.path.join(out_folder, "mutrates_outgroup1.out"),  out_folder=out_folder, window_size=1000, haploid=False, weights=None)

        all_segments.to_csv(output.prob_file, index=False)

        inferred_tract_files = process_output(all_segments, out_folder, output_prefix, src_id, cutoff_list=cutoff_list, return_filenames=True)

        for inferred_tracts in inferred_tract_files:
            precision, recall = cal_accuracy(input.true_tracts, inferred_tracts[1])
            cutoff = inferred_tracts[0]

            prec_rec_file = inferred_tracts[1].rsplit('.', 1)[0] + ".accuracy"

            with open(prec_rec_file, 'w') as o:
                o.write(f'{demog_id}\tnref_{nref}_ntgt_{ntgt}\t{cutoff}\t{precision}\t{recall}\n')

