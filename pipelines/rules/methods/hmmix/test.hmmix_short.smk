import os

#-----------------------------------------------------------------------------------------------------------------------
# run section of rule process_test_data in pipelines/rules/methods/hmmix/test.hmmix.smk sent to process_hmmix.py
#-----------------------------------------------------------------------------------------------------------------------

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
    conda:
        "../../../envs/hmmix-env.yaml",
    log:
        "logs/hmmix/proctestdata.{seed}.log",
    benchmark:
        "benchmarks/hmmix/proctestdata.{seed}.benchmark.txt",
    params:
        skov_output_dir = skov_output_dir,
        ref_id = ref_set,
        ref_id = config["ref_id"],
        src_id = config["src_id"]
    shell:
        """    
        python pipelines/rules/methods/hmmix/process_hmmix.py \
            --inref_list '{input.ref_list}' \
            --intgt_list '{input.tgt_list}' \
            --in_vcf '{input.vcf}' \
            --in_true_tracts '{input.true_tracts}' \
            --skov_output_dir '{params.skov_output_dir}' \
            --seed {wildcards.seed} \
            --out_outgroup_file '{output.outgroup_file}' \
            --out_prob_file '{output.prob_file}' \
            --ref_set {params.ref_set} \
            --ref_id '{params.ref_id}' \
            --src_id '{params.src_id}' 
        """