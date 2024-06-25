import os
import numpy as np

### CONFIG ###
#configfile for sprime
configfile: "config/sprime/config_sprime.yaml"

#use binary_model
binary_model = False
#remove all multiallelic sites (in case of Jukes-Cantor-model)
only_biallelic = True

#there is no training part for SPrime, so we use the test params
params_set = "test"

nref = config["nref"]
ntgt = config["ntgt"]
ploidy = config["ploidy"]
win_step = config["win_step"]
cutoff_num = config["cutoff_num"]
cutoff_list = np.round(np.linspace(0, 1, cutoff_num, endpoint=False), 2)
cutoff_list = np.append(cutoff_list, [0.99, 0.999])

output_prefix = config["output_prefix"]
nrep = config["nrep"][params_set]
seq_len = config["seq_len"][params_set]
demog_id = config["demog_id"][params_set]
demes_file = config["demes"][params_set]
mut_rate = config["mut_rate"][params_set]
rec_rate = config["rec_rate"][params_set]
ploidy = config["ploidy"]
ref_id = config["ref_id"][params_set]
tgt_id = config["tgt_id"][params_set]
src_id = config["src_id"][params_set]

map_file_config = None

np.random.seed(config["seed"])
seed_list = np.random.random_integers(1, 2**31, nrep)

#thresholds for sprime
threshold_list = config["threshold_list"]
#sprime executable
sprime_exec = config["sprime_exec"]


output_dir = f'results/data/{params_set}/{demog_id}/nref_{nref}/ntgt_{ntgt}'
sprime_output_dir = f'results/sprime/{params_set}/{demog_id}/nref_{nref}/ntgt_{ntgt}'


##### Target rules #####


rule all:
    input:
        expand(sprime_output_dir + "/{demog}/nref_{nref}/ntgt_{ntgt}/{seed}/{threshold}/sprime.1src.out.{threshold}.bed", demog=demog_id, nref=nref, ntgt=ntgt, seed = seed_list, threshold = threshold_list)


##### Modules #####


include: "../rules/commons/simulation_mp_extended.smk"
include: "../rules/methods/SPrime/test.SPrime.smk"
#include: "../rules/commons/plot.smk"
