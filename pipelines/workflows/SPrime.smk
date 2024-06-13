import os
import numpy as np


##### Configs #####


feature_id = config["feature_id"]
feature_config = config["feature_config"]
nfeature = config["nfeature"]
ploidy = config["ploidy"]
output_prefix = config["output_prefix"]
win_step = config["win_step"]
nref = config["nref"]
ntgt = config["ntgt"]
geno_state_list = config["geno_states"]
dataset_state_list = config["dataset_states"]
cutoff_num = config["cutoff_num"]
cutoff_list = np.round(np.linspace(0, 1, cutoff_num, endpoint=False), 2)
cutoff_list = np.append(cutoff_list, [0.99, 0.999])

nrep = {}
seq_len = {}
demog_id = {}
demes = {}
mut_rate = {}
rec_rate = {}
ref_id = {}
tgt_id = {}
src_id = {}
seed_list = {}
output_dir = {}

for k in ["train", "test"]: 
    np.random.seed(config["seed"])
    nrep[k] = config["nrep"][k]
    seq_len[k] = config["seq_len"][k]
    demog_id[k] = config["demog_id"][k]
    demes[k] = config["demes"][k]
    mut_rate[k] = config["mut_rate"][k]
    rec_rate[k] = config["rec_rate"][k]
    ref_id[k] = config["ref_id"][k]
    tgt_id[k] = config["tgt_id"][k]
    src_id[k] = config["src_id"][k]
    seed_list[k] = np.random.randint(1, 2**31, nrep[k])
    output_dir[k] = f"results/data/{k}/{demog_id[k]}/nref_{nref}/ntgt_{ntgt}"

train_demog = demog_id["train"]
test_demog = demog_id["test"]
performance_dir = f"results/performance/train_{train_demog}_test_{test_demog}"


##### Target rules #####


rule all:
    input:
        expand("results/plots/train_{train_demog}_test_{test_demog}/{output_prefix}.performance.pdf", 
               train_demog=train_demog, test_demog=test_demog, output_prefix=output_prefix),


##### Modules #####


include: "../rules/commons/simulation.smk"
include: "../rules/methods/SPrime/test.SPrime.smk"
include: "../rules/commons/plot.smk"
