import os
import demes
import numpy as np
import pandas as pd

from sstar_additional_functions import *


'''
rule all:
    input:
        sstar_output_dir + "/sstar_1src_accuracy.txt",
'''

rule mut_rec_combination:
    output:
        rates = sstar_output_dir_simulation + "/{scenario}/snps/rates.combination",
    params:
        seq_len = 50000,
        mut_rate = config["mut_rate"],
        rec_rate = config["rec_rate"],
    resources: time_min=60, mem_mb=5000, cpus=1,
    threads: 1,
    conda:
        "../../../envs/sstar-env.yaml",
    log:
        "logs/sstar/mutrec.{scenario}.log",
    benchmark:
        "benchmarks/sstar/mutrec.{scenario}.benchmark.txt",
    shell:
        """    
        python pipelines/rules/methods/sstar/mut_rec.py \
         --seqlen {params.seq_len}  --mutrate {params.mut_rate} \
         --recrate {params.rec_rate} --outfile '{output.rates}'
        """

rule simulate_glm_data:
    input:
        rates = rules.mut_rec_combination.output.rates,
    output:
        ms = 'results/sstar/{params_set}/{demog}/nref_{nref}/ntgt_{ntgt}' +"/simulation/{scenario}/snps/{snp_num}/sim1src.ms",
    log:
        "logs/sstar/msglm.{scenario}.{params_set}.{demog}.{nref}.{ntgt}.{snp_num}.log",
    benchmark:
        "benchmarks/sstar/msglm.{scenario}.{params_set}.{demog}.{nref}.{ntgt}.{snp_num}.benchmark.txt",
    params:
        nsamp = lambda wildcards: 2*(int(wildcards.nref)+int(1)),
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


#-----------------------------------------------------------------------------------------------------------------------
# create ind file lists for ms simulations
rule ms_ind_lists:
    output:
        ss_ref = 'results/sstar/{params_set}/{demog}/nref_{nref}/ntgt_{ntgt}' +"/simulation/{scenario}/sstarsim.ref.ind.list",
        ss_tgt = 'results/sstar/{params_set}/{demog}/nref_{nref}/ntgt_{ntgt}' +"/simulation/{scenario}/sstarsim.tgt.ind.list",
    log:
        "logs/sstar/ms_ind_lists.{scenario}.{params_set}.{demog}.{nref}.{ntgt}.log",
    benchmark:
        "benchmarks/sstar/ms_ind_lists.{scenario}.{params_set}.{demog}.{nref}.{ntgt}.benchmark.txt",
    params:
        nsamp = lambda wildcards: 2*(int(wildcards.nref)+int(1)),
        ploidy = config["ploidy"]
    resources: time_min=60, mem_mb=5000, cpus=1,
    threads: 1,
    run:
        create_ind_lists(params.nsamp, output.ss_ref, output.ss_tgt, params.ploidy, ind_prefix="tsk_")

rule ms2vcf:
    input:
        ms = rules.simulate_glm_data.output.ms,
    output:
        vcf = 'results/sstar/{params_set}/{demog}/nref_{nref}/ntgt_{ntgt}' +"/simulation/{scenario}/snps/{snp_num}/sim1src.vcf",
    log:
        "logs/sstar/ms2vcf.{scenario}.{params_set}.{demog}.{nref}.{ntgt}.{snp_num}.log",
    benchmark:
        "benchmarks/sstar/ms2vcf.{scenario}.{params_set}.{demog}.{nref}.{ntgt}.{snp_num}.benchmark.txt",
    params:
        nsamp = lambda wildcards: 2*(int(wildcards.nref)+int(1)),
        seq_len = 50000,
        ploidy = config["ploidy"]
    resources: 
        time = lambda wildcards: 360 if (int(wildcards.snp_num) < 340) else 1000,
        mem_gb = lambda wildcards: 10 if (int(wildcards.snp_num) < 340) else 100,
        cpus = 1
    threads: 1,
    run:
        ms2vcf(input.ms, output.vcf, params.nsamp, params.seq_len, params.ploidy, ind_prefix="tsk_")
#-----------------------------------------------------------------------------------------------------------------------


rule cal_score:
    input:
        vcf = rules.ms2vcf.output.vcf,
	ref_list = rules.ms_ind_lists.output.ss_ref,
	tgt_list = rules.ms_ind_lists.output.ss_tgt,
    output:
        score = 'results/sstar/{params_set}/{demog}/nref_{nref}/ntgt_{ntgt}' +"/simulation/{scenario}/snps/{snp_num}/sim1src.sstar.scores",
    params:
        seq_len = 50000,
    conda:
        "../../../envs/sstar-env.yaml",
    benchmark:
        "benchmarks/sstar/sstarcalcscore.{scenario}.{params_set}.{demog}.{nref}.{ntgt}.{snp_num}.benchmark.txt",
    log:
        "logs/sstar/sstarcalcscore.{scenario}.{params_set}.{demog}.{nref}.{ntgt}.{snp_num}.log",
    resources:
        time = lambda wildcards: 360 if (int(wildcards.snp_num) < 340) else 1000,
        mem_gb = lambda wildcards: 10 if (int(wildcards.snp_num) < 340) else 100,
        cpus = 1
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
    log:
        "logs/sstar/sstarcalcquant.{scenario}.{params_set}.{demog}.{nref}.{ntgt}.{snp_num}.log",
    benchmark:
        "benchmarks/sstar/sstarcalcquant.{scenario}.{params_set}.{demog}.{nref}.{ntgt}.{snp_num}.benchmark.txt",
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
    log:
        "logs/sstar/sstarquantsum.{scenario}.{params_set}.{demog}.{nref}.{ntgt}.log",
    benchmark:
        "benchmarks/sstar/sstarquantsum.{scenario}.{params_set}.{demog}.{nref}.{ntgt}.benchmark.txt",
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
    conda:
        "../../../envs/sstar-env.yaml",
    log:
        "logs/sstar/sstarscore.{params_set}.{demog}.{nref}.{ntgt}.{seed}.log",
    benchmark:
        "benchmarks/sstar/sstarscore.{params_set}.{demog}.{nref}.{ntgt}.{seed}.benchmark.txt",
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
    conda:
        "../../../envs/sstar-env.yaml",
    log:
        "logs/sstar/sstarthresh.{params_set}.{demog}.{nref}.{ntgt}.{seed}.{scenario}.{quantile}.log",
    benchmark:
        "benchmarks/sstar/sstarthresh.{params_set}.{demog}.{nref}.{ntgt}.{seed}.{scenario}.{quantile}.benchmark.txt",
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
    log:
        "logs/sstar/sstarprocout.{params_set}.{demog}.{nref}.{ntgt}.{seed}.{scenario}.{quantile}.log",
    benchmark:
        "benchmarks/sstar/sstarprocout.{params_set}.{demog}.{nref}.{ntgt}.{seed}.{scenario}.{quantile}.benchmark.txt",
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
    log:
        "logs/sstar/accuracy.log",
    benchmark:
        "benchmarks/sstar/accuracy.benchmark.txt",
    shell:
        """
        cat {input.accuracy_files} | sed '1idemography\\tscenario\\tsample\\tcutoff\\tprecision\\trecall' > {output.accuracy_table}
        """

