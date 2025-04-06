import os
import numpy as np
import random
#-----------------------------------------------------------------------------------------------------------------------
#configfile for msprime simulations
configfile: "config/scenarios/config_simulate_humnea.yaml"

params_set = "test"

nref = config["nref"]
ntgt = config["ntgt"]
ploidy = config["ploidy"]
win_step = config["win_step"]

output_prefix = config["output_prefix"]
nrep = config["nrep"]
seq_len = config["seq_len"]
demog_id = config["demog_id"]
demes_file = config["demes_file"]
mut_rate = config["mut_rate"]
rec_rate = config["rec_rate"]
ref_id = config["ref_id"]
tgt_id = config["tgt_id"]
src_id = config["src_id"]

#-----------------------------------------------------------------------------------------------------------------------
output_dir = f'results/data/{params_set}/{demog_id}/nref_{nref}/ntgt_{ntgt}'

seed_list = np.random.random_integers(1, 2**31, nrep)
# write out the seeds to file *
np.save('results/seeds/seed_list.npy', seed_list)

# then in the method workflows read in the relevant seeds:
# seed_list = np.load('results/seeds/seed_list.npy')
#-----------------------------------------------------------------------------------------------------------------------
#if binary==True in config file, do simulations with binary mutation model
try:
    binary = config["binary"]
except KeyError:
    binary = False
#-----------------------------------------------------------------------------------------------------------------------
rule all:
    input:
        expand("{output_dir}/{seed}/{output_prefix}.vcf.gz", output_dir=output_dir, seed=seed_list, output_prefix=output_prefix),
        expand("{output_dir}/{seed}/{output_prefix}.truth.tracts.bed", output_dir=output_dir, seed=seed_list, output_prefix=output_prefix),
        expand("{output_dir}/{seed}/{output_prefix}.ref.ind.list", output_dir=output_dir, seed=seed_list, output_prefix=output_prefix),
        expand("{output_dir}/{seed}/{output_prefix}.tgt.ind.list", output_dir=output_dir, seed=seed_list, output_prefix=output_prefix)

include: "../rules/commons/simulation_mp_extended_short.smk"

