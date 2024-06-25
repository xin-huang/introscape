import os
import numpy as np


try:
    libr_dir = os.environ["CONDA_PREFIX"] + "/lib/R/library"
except:
    libr_dir = os.environ["CONDA"] + "/lib/R/library"


### CONFIG ###
#configfile for sstar
configfile: "config/sstar/config_sstar.yaml"

#there is no (real) training part for sstar, so we use the test params
params_set = "test"

nref = config["nref"]
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

#np.random.seed(config["seed"])
#seed_list = np.random.random_integers(1, 2**31, nrep)

snp_num_list = np.arange(25,705,5)
quantile_list = config["quantiles"]

output_dir = f'results/data/{params_set}/{demog_id}/nref_{nref}/ntgt_{ntgt}'
sstar_output_dir = f'results/sstar/{params_set}/{demog_id}/nref_{nref}/ntgt_{ntgt}'
#output folder for the simulations used to compute the sstar model
sstar_output_dir_simulation = os.path.join(sstar_output_dir, "simulation")

#only in case that no new data is simulated
#seed_list = list_subdirectories(output_dir)
seed_list = np.random.random_integers(1, 2**31, nrep)


new_params = config["sstar_ms_params"]
scenario_list = config["scenarios"]
score_seqlen = config["score_seqlen"]

#if binary==True in config file, do simulations with binary mutation model
try:
    binary = config["binary"]
except KeyError:
    binary = False


##### Target rules #####


rule all:
    input:
        sstar_output_dir + "/sstar_1src_accuracy.txt",


##### Modules #####


include: "../rules/commons/simulation_mp_extended.smk"
include: "../rules/methods/sstar/test.sstar.smk"
#include: "../rules/commons/plot.smk"
