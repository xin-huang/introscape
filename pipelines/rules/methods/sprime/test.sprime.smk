import numpy as np
import os

from sprime_additional_functions import *


### CONFIG ###
#configfile for sprime
configfile: "config/sprime/config_sprime.yaml"

binary_model = False
only_biallelic = True

params_set = "test"

feature_id = config["feature_id"]
feature_config = config["feature_config"]
nfeature = config["nfeature"]
nref = config["nref"]
ntgt = config["ntgt"]
ploidy = config["ploidy"]
geno_state_list = config["geno_states"]
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

output_dir = f'results/data/{params_set}/{demog_id}/nref_{nref}/ntgt_{ntgt}'
sprime_output_dir = f'results/sprime/{params_set}/{demog_id}/nref_{nref}/ntgt_{ntgt}'


threshold_list = config["threshold_list"]
sprime_exec = config["sprime_exec"]

rule all:
    input:
        expand(sprime_output_dir + "/{demog}/nref_{nref}/ntgt_{ntgt}/{seed}/{threshold}/sprime.1src.out.{threshold}.bed", demog=demog_id, nref=nref, ntgt=ntgt, seed = seed_list, threshold = threshold_list)

rule sprime_run:
    input:
        vcf = output_dir + "/{seed}/" + output_prefix + ".vcf.gz",
        ref_list = output_dir + "/{seed}/" + output_prefix + ".ref.ind.list",
    output:
        log = os.path.join(sprime_output_dir, "{demog}/nref_{nref}/ntgt_{ntgt}/{seed}/{threshold}/sprime.1src.out.{threshold}.log"),
        score = os.path.join(sprime_output_dir, "{demog}/nref_{nref}/ntgt_{ntgt}/{seed}/{threshold}/sprime.1src.out.{threshold}.score"),
    params:
        sprime_exec = config["sprime_exec"],

        threshold = lambda wildcards: wildcards.threshold,
        output_prefix = os.path.join(sprime_output_dir, "{demog}/nref_{nref}/ntgt_{ntgt}/{seed}/{threshold}/sprime.1src.out.{threshold}")

    resources: time_min=60, mem_mb=2000, cpus=1,
    threads: 1,
    run:
        if map_file_config is None:
            recomb_value = seq_len * rec_rate * 100
            create_map_file(recomb_value, seq_len, map_file = os.path.join(sprime_output_dir, "sim.map"))
            map_file = os.path.join(sprime_output_dir, "sim.map")
        else:
            map_file = map_file_config

        shell("java -Xmx2g -jar {params.sprime_exec} gt={input.vcf} outgroup={input.ref_list} map={map_file} out={params.output_prefix} minscore={params.threshold} mu={mut_rate}")


rule sprime_process_output:
    input:
        scores = rules.sprime_run.output.score,
        true_tracts = output_dir + "/{seed}/" + output_prefix + ".truth.tracts.bed",
    output:
        inferred_tracts = os.path.join(sprime_output_dir, "{demog}/nref_{nref}/ntgt_{ntgt}/{seed}/{threshold}/sprime.1src.out.{threshold}.bed"),
        accuracy = os.path.join(sprime_output_dir, "{demog}/nref_{nref}/ntgt_{ntgt}/{seed}/{threshold}/sprime.1src.out.{threshold}.accuracy"),
    resources: time_min=60, mem_mb=2000, cpus=1,
    threads: 1,
    run:
       process_sprime_output(input.scores, output.inferred_tracts)
       precision, recall = cal_accuracy(input.true_tracts, output.inferred_tracts)
       with open(output.accuracy, 'w') as o:
           o.write(f'{wildcards.demog}\tnref_{wildcards.nref}_ntgt_{wildcards.ntgt}\t{wildcards.threshold}\t{precision}\t{recall}\n')
