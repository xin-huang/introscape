import os
import demes
import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy.stats import nbinom

from sstar_additional_functions import *

try:
    libr_dir = os.environ["CONDA_PREFIX"] + "/lib/R/library"
except:
    libr_dir = os.environ["CONDA"] + "/lib/R/library"

### CONFIG ###
#configfile for sstar
configfile: "config/sstar/config_sstar.yaml"

params_set = "test"

feature_id = config["feature_id"]
feature_config = config["feature_config"]
nfeature = config["nfeature"]
nref = config["nref"]
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
seq_len = config["seq_len"]
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
sstar_output_dir_simulation = os.path.join(sstar_output_dir, "simulation")

seed_list = list_subdirectories(output_dir)

new_params = config["sstar_ms_params"]
scenario_list = config["scenarios"]
score_seqlen = config["score_seqlen"]


rule all:
    input:
        sstar_output_dir + "/sstar_1src_accuracy.txt",


rule mut_rec_combination:
    output:
        rates = sstar_output_dir_simulation + "/{scenario}/snps/rates.combination",
    resources: time_min=60, mem_mb=5000, cpus=1,
    threads: 1,
    run:
        seq_len = 50000

        N0 = 1000

    	mut_rate_mean = float(mut_rate)
     	rec_rate_mean = float(rec_rate)

        scaled_mut_rate_mean = 4*N0*mut_rate_mean*seq_len
        scaled_mut_rate_sdv = 0.233

        scaled_rec_rate_mean = 4*N0*rec_rate_mean*seq_len
        mut_rate_list = norm.rvs(loc=scaled_mut_rate_mean, scale=scaled_mut_rate_sdv, size=20000)
        rec_rate_list = nbinom.rvs(n=0.5, p=0.5/(0.5+scaled_rec_rate_mean), size=20000)

        with open(output.rates, 'w') as o:
            for i in range(len(mut_rate_list)):
                if mut_rate_list[i] < 0.001: mut_rate_list[i] = 0.001
                if rec_rate_list[i] < 0.001: rec_rate_list[i] = 0.001
                mut_rate_new = mut_rate_list[i]
                rec_rate_new = rec_rate_list[i]
                o.write(f'{mut_rate_new}\t{rec_rate_new}\n')


rule simulate_glm_data:
    input:
        rates = rules.mut_rec_combination.output.rates,
    output:
        ms = 'results/sstar/{params_set}/{demog}/nref_{nref}/ntgt_{ntgt}' +"/simulation/{scenario}/snps/{snp_num}/sim1src.ms",

    params:
        nsamp = lambda wildcards: 2*(int(1)+int(wildcards.ntgt)),
        #nsamp = lambda wildcards: 2*(int(wildcards.nref)+int(wildcards.ntgt)),
        #nreps = 10,
        nreps = 20000,

        seq_len = 50000,

        ms_exec = config["ms_exec"],

        ms_params = lambda wildcards: demes.to_ms(demes.load(new_params[wildcards.scenario]["yaml"]), N0=1000, samples=new_params[wildcards.scenario]["samples"]),

    resources: time_min=60, mem_mb=5000, cpus=1,
    threads: 1,
    shell:
        """
        cat {input.rates} | {params.ms_exec} {params.nsamp} {params.nreps} -t tbs -r tbs {params.seq_len} -s {wildcards.snp_num} {params.ms_params} > {output.ms}
        """


rule ms2vcf:
    input:
        ms = rules.simulate_glm_data.output.ms,
    output:
        vcf = 'results/sstar/{params_set}/{demog}/nref_{nref}/ntgt_{ntgt}' +"/simulation/{scenario}/snps/{snp_num}/sim1src.vcf",
        #create lists for sstar simulations
        ss_ref = 'results/sstar/{params_set}/{demog}/nref_{nref}/ntgt_{ntgt}' +"/simulation/{scenario}/snps/{snp_num}/sstarsim.ref.ind.list",
        ss_tgt = 'results/sstar/{params_set}/{demog}/nref_{nref}/ntgt_{ntgt}' +"/simulation/{scenario}/snps/{snp_num}/sstarsim.tgt.ind.list",

    params:
        nsamp = lambda wildcards: 2*(int(1)+int(wildcards.ntgt)),

        seq_len = 50000,
    resources: time_min=60, mem_mb=5000, cpus=1,
    threads: 1,
    run:
        ms2vcf_create_ind_lists(input.ms, output.vcf, params.nsamp, params.seq_len, output.ss_ref, output.ss_tgt, ind_prefix="tsk_")


rule cal_score:
    input:
        vcf = rules.ms2vcf.output.vcf,
        #ref_list = output_dir + "/" + str(seed_list[0]) + "/" + output_prefix + ".ref.ind.list",
        #tgt_list = output_dir + "/" + str(seed_list[0]) + "/" + output_prefix + ".tgt.ind.list",
        ref_list = rules.ms2vcf.output.ss_ref,
        tgt_list = rules.ms2vcf.output.ss_tgt,
    output:
        score = 'results/sstar/{params_set}/{demog}/nref_{nref}/ntgt_{ntgt}' +"/simulation/{scenario}/snps/{snp_num}/sim1src.sstar.scores",
    params:
        seq_len = 50000,
    resources: time_min=3000, mem_mb=10000, cpus=1,
    threads: 1,
    shell:
        """
        sstar score --vcf {input.vcf} --ref {input.ref_list} --tgt {input.tgt_list} --output {output.score} --thread {threads} --win-len {params.seq_len} --win-step {params.seq_len}
        """


rule cal_quantile:
    input:
        score = rules.cal_score.output.score,
    output:
        quantile = 'results/sstar/{params_set}/{demog}/nref_{nref}/ntgt_{ntgt}' + "/simulation/{scenario}/snps/{snp_num}/sim1src.sstar.quantile",
    params:
        sim_quantiles = quantile_list,
    resources: time_min=3000, mem_mb=5000, cpus=1,
    threads: 1,
    run:
        df = pd.read_csv(input.score, sep="\t").dropna()
        mean_df = df.groupby(['chrom', 'start', 'end'], as_index=False)['S*_score'].mean().dropna()
        scores = np.quantile(mean_df['S*_score'], params.sim_quantiles)
        with open(output.quantile, 'w') as o:
            for i in range(len(scores)):
                o.write(f'{scores[i]}\t{wildcards.snp_num}\t{params.sim_quantiles[i]}\n')


rule quantile_summary:
    input:
        res = expand('results/sstar/{params_set}/{demog}/nref_{nref}/ntgt_{ntgt}' +  "/simulation/{scenario}/snps/{snp_num}/sim1src.sstar.quantile", snp_num=snp_num_list, demog=demog_id, params_set=params_set, nref=nref, ntgt=ntgt, scenario=scenario_list),
    output:
        output_res = 'results/sstar/{params_set}/{demog}/nref_{nref}/ntgt_{ntgt}' +  "/simulation/{scenario}/quantile.1src.summary.txt",
    resources: time_min=3000, mem_mb=5000, cpus=1,
    threads: 1,
    shell:
        """
        cat {input.res} | sort -nk 2,2 | sed '1iS*_score\\tSNP_number\\tquantile\n' > {output.output_res}

        """


rule sstar_score:
    input:
        vcf = output_dir + "/{seed}/" + output_prefix + ".vcf.gz",
        ref_ind = output_dir + "/{seed}/" + output_prefix + ".ref.ind.list",
        tgt_ind = output_dir + "/{seed}/" + output_prefix + ".tgt.ind.list"
    output:
        score = "results/sstar/{params_set}/{demog}/nref_{nref}/ntgt_{ntgt}/{seed}/sstar.1src.out.score",
    params:
        seq_len = score_seqlen,
    resources: time_min=3000, mem_mb=10000, cpus=1,
    threads: 1,
    shell:
        """
        sstar score --vcf {input.vcf} --ref {input.ref_ind} --tgt {input.tgt_ind} --output {output.score} --thread {threads} --win-len {params.seq_len} --win-step 10000
        """


rule sstar_threshold:
    input:
        score = rules.sstar_score.output.score,
        summary = 'results/sstar/{params_set}/{demog}/nref_{nref}/ntgt_{ntgt}' +  "/simulation/{scenario}/quantile.1src.summary.txt",
    output:
        quantiles = "results/sstar/{params_set}/{demog}/nref_{nref}/ntgt_{ntgt}" + "/{seed}/{scenario}/sstar.1src.quantile.{quantile}.out",

        #quantiles = output_dir + "inference/sstar/{demog}/nref_{nref}/ntgt_{ntgt}/{seed}/{scenario}/sstar.1src.quantile.{quantile}.out",
    params:
        R_LIBS = libr_dir,
    resources: time_min=3000, mem_mb=10000, cpus=1,
    threads: 1,
    shell:
        """
        export R_LIBS={params.R_LIBS}
        sstar threshold --score {input.score} --sim-data {input.summary} --quantile {wildcards.quantile} --output {output.quantiles}
        """


rule sstar_process_output:
    input:
        quantiles = rules.sstar_threshold.output.quantiles,
        true_tracts = output_dir + "/{seed}/" + output_prefix + ".truth.tracts.bed",
    output:
        inferred_tracts = "results/sstar/{params_set}/{demog}/nref_{nref}/ntgt_{ntgt}" + "/{seed}/{scenario}/sstar.1src.quantile.{quantile}.out.bed",
        accuracy = "results/sstar/{params_set}/{demog}/nref_{nref}/ntgt_{ntgt}" + "/{seed}/{scenario}/sstar.1src.quantile.{quantile}.out.accuracy",
    resources: time_min=3000, mem_mb=10000, cpus=1,
    threads: 1,
    run:
        process_sstar_1src_output(input.quantiles, output.inferred_tracts)
        precision, recall = cal_accuracy(input.true_tracts, output.inferred_tracts)
        with open(output.accuracy, 'w') as o:
            o.write(f'{wildcards.demog}\t{wildcards.scenario}\tnref_{wildcards.nref}_ntgt_{wildcards.ntgt}\t{wildcards.quantile}\t{precision}\t{recall}\n')


rule accuracy_summary:
    input:
        accuracy_files = expand("results/sstar/{params_set}/{demog}/nref_{nref}/ntgt_{ntgt}" + "/{seed}/{scenario}/sstar.1src.quantile.{quantile}.out.accuracy",
                                nref=nref, ntgt=ntgt, seed=seed_list, quantile=quantile_list, scenario=scenario_list, params_set=params_set, demog=demog_id),
    output:
        accuracy_table = sstar_output_dir + "/sstar_1src_accuracy.txt",
    resources: time_min=60, mem_mb=2000, cpus=1,
    threads: 1,
    shell:
        """
        cat {input.accuracy_files} | sed '1idemography\\tscenario\\tsample\\tcutoff\\tprecision\\trecall' > {output.accuracy_table}
        """

