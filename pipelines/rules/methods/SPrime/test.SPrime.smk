import numpy as np
import os

from sprime_additional_functions import *


'''
rule all:
    input:
        expand(sprime_output_dir + "/{demog}/nref_{nref}/ntgt_{ntgt}/{seed}/{threshold}/sprime.1src.out.{threshold}.bed", demog=demog_id, nref=nref, ntgt=ntgt, seed = seed_list, threshold = threshold_list)
'''

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
    log:
        "logs/sprime/sprimerun.{demog}.{nref}.{ntgt}.{seed}/{threshold}.log",
    benchmark:
        "benchmarks/sprime/sprimerun.{demog}.{nref}.{ntgt}.{seed}/{threshold}.benchmark.txt",
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
    log:
        "logs/sprime/sprimeprocout.{demog}.{nref}.{ntgt}.{seed}/{threshold}.log",
    benchmark:
        "benchmarks/sprime/sprimeprocout.{demog}.{nref}.{ntgt}.{seed}/{threshold}.benchmark.txt",
    resources: time_min=60, mem_mb=2000, cpus=1,
    threads: 1,
    run:
       process_sprime_output(input.scores, output.inferred_tracts)
       precision, recall = cal_accuracy(input.true_tracts, output.inferred_tracts)
       with open(output.accuracy, 'w') as o:
           o.write(f'{wildcards.demog}\tnref_{wildcards.nref}_ntgt_{wildcards.ntgt}\t{wildcards.threshold}\t{precision}\t{recall}\n')


rule sprime_accuracy:
    input:
        accuracy_files = expand(sprime_output_dir + "/{demog}/nref_{nref}/ntgt_{ntgt}/{seed}/{threshold}/sprime.1src.out.{threshold}.accuracy", demog=demog_id, nref=nref, ntgt=ntgt, seed=seed_list, threshold=threshold_list),
    output:
        accuracy_table = os.path.join(sprime_output_dir + "{demog}/nref_{nref}/ntgt_{ntgt}/{seed}/sprime_accuracy.txt"),
    log:
        "logs/sprime/sprime_accuracy.{demog}.{nref}.{ntgt}.{seed}.log",
    benchmark:
        "benchmarks/sprime/sprime_accuracy.{demog}.{nref}.{ntgt}.{seed}.benchmark.txt",
    resources: time_min=60, mem_mb=2000, cpus=1,
    threads: 1,
    shell:
        """
        cat {input.accuracy_files} | sed '1idemography\\tscenario\\tsample\\tcutoff\\tprecision\\trecall' > {output.accuracy_table}
        """

