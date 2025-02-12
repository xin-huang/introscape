import os
import numpy as np


### CONFIG ###
#configfile: "config/test_gorilla.config.yaml"

#there is no training part for hmmix, so we use the test params
params_set = "test"

nref = config["nref"]
ntgt = config["ntgt"]
ploidy = config["ploidy"]
win_step = config["win_step"]
cutoff_num = config["cutoff_num"]
cutoff_list = np.round(np.linspace(0, 1, cutoff_num, endpoint=False), 2)
cutoff_list = np.append(cutoff_list, [0.99, 0.999])

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

#if binary==True in config file, do simulations with binary mutation model
try:
    binary = config["binary"]
except KeyError:
    binary = False

if binary:
    ref_set = ["0", "1"]
else:
    ref_set = ["A", "C", "G", "T"]


np.random.seed(config["seed"])
seed_list = np.random.random_integers(1, 2**31, nrep)


output_dir = f'results/data/{params_set}/{demog_id}/nref_{nref}/ntgt_{ntgt}'
skov_output_dir = f'results/SkovHMM/{params_set}/{demog_id}/nref_{nref}/ntgt_{ntgt}'



##### Target rules #####


rule all:
    input:
        expand(skov_output_dir + "/{seed}/output_ref_vs.txt", demog_id=demog_id, nref=nref, ntgt=ntgt, seed = seed_list)


##### Modules #####


#include: "../rules/commons/simulation_mp_extended.smk"
#include: "../rules/methods/hmmix/test.hmmix.smk"
#include: "../rules/commons/plot.smk"

#-----------------------------------------------------------------------------------------------------------------------
# testing
include: "../rules/commons/simulation_mp_extended_short.smk"
include: "../rules/methods/hmmix/test.hmmix_short.smk"

