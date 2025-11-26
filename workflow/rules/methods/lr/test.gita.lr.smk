rule process_test_data:
    input:
        vcf = rules.compress_vcf.output.vcf,
        ref = rules.simulate_test_data.output.ref,
        tgt = rules.simulate_test_data.output.tgt,
        feature_config = feature_config,
    output:
        features = "results/data/test/{test_demog}/nref_{nref}/ntgt_{ntgt}/{geno_state}/{test_seed}/{output_prefix}.features.gz",
    params:
        seq_len = seq_len["train"],
        win_step = win_step,
        is_phased = lambda wildcards: '--phased' if wildcards.geno_state == 'phased' else '',
        output_dir = os.path.join(output_dir["test"], "{geno_state}/{test_seed}"),
        features = "results/data/test/{test_demog}/nref_{nref}/ntgt_{ntgt}/{geno_state}/{test_seed}/{output_prefix}.features",
    log:
        "logs/process_test_data/{test_demog}/nref_{nref}/ntgt_{ntgt}/{geno_state}/{test_seed}/{output_prefix}.log",
    resources:
        partition = "himem,gpu",
        mem_gb = 32,
        cpus = 16,
    shell:
        """
        gita lr preprocess --vcf {input.vcf} --ref {input.ref} --tgt {input.tgt} --features {input.feature_config} \
                           --win-len {params.seq_len} --win-step {params.win_step} --output-prefix {wildcards.output_prefix} --output-dir {params.output_dir} \
                           --worker {resources.cpus} {params.is_phased}
        gzip -c {params.features} > {output.features} 2> {log}
        rm {params.features} 2>> {log}
        """


rule test_logistic_regression_model:
    input:
        features = rules.process_test_data.output.features,
        model_file = rules.train_logistic_regression_model.output.model_file,
    output:
        pred = "results/performance/train_{train_demog}_test_{test_demog}/nref_{nref}/ntgt_{ntgt}/{geno_state}/train_{train_seed}_test_{test_seed}/{output_prefix}.lr.pred.gz",
    params:
        pred = "results/performance/train_{train_demog}_test_{test_demog}/nref_{nref}/ntgt_{ntgt}/{geno_state}/train_{train_seed}_test_{test_seed}/{output_prefix}.lr.pred",
        output_dir = "results/performance/train_{train_demog}_test_{test_demog}/nref_{nref}/ntgt_{ntgt}/{geno_state}/train_{train_seed}_test_{test_seed}",
    log:
        "logs/test_logistic_regression_model/train_{train_demog}_test_{test_demog}/nref_{nref}/ntgt_{ntgt}/{geno_state}/train_{train_seed}_test_{test_seed}/{output_prefix}.log",
    resources:
        partition = "himem,gpu",
        mem_gb = 100,
    shell:
        """
        gita lr infer --features {input.features} --model-file {input.model_file} --output-prefix {output_prefix} --output-dir {params.output_dir}
        gzip -c {params.pred} > {output.pred} 2> {log}
        rm {params.pred} 2>> {log}
        """


rule evaluate_logistic_regression_model:
    input:
        pred = rules.test_logistic_regression_model.output.pred,
        true_tracts = lambda wildcards: rules.get_phased_true_tracts.output.bed if wildcards.geno_state == 'phased' else rules.get_unphased_true_tracts.output.bed,
    output:
        inferred_tracts = "results/performance/train_{train_demog}_test_{test_demog}/nref_{nref}/ntgt_{ntgt}/{geno_state}/train_{train_seed}_test_{test_seed}/{output_prefix}.lr.{cutoff}.inferred.tracts.bed",
        performance = "results/performance/train_{train_demog}_test_{test_demog}/nref_{nref}/ntgt_{ntgt}/{geno_state}/train_{train_seed}_test_{test_seed}/{output_prefix}.lr.{cutoff}.performance",
        summary = "results/performance/train_{train_demog}_test_{test_demog}/nref_{nref}/ntgt_{ntgt}/{geno_state}/train_{train_seed}_test_{test_seed}/{output_prefix}.lr.{cutoff}.performance.summary",
    params:
        feature_id = feature_id,
    log:
        "logs/evaluate_logistic_regression_model/train_{train_demog}_test_{test_demog}/nref_{nref}/ntgt_{ntgt}/{geno_state}/train_{train_seed}_test_{test_seed}/{output_prefix}.{cutoff}.log",
    resources:
        partition = "himem,gpu",
    shell:
        """
        zcat {input.pred} | sed '1d' | awk '$NF>{wildcards.cutoff}' | awk 'BEGIN{{OFS="\\t"}}{{print $1,$2,$3,$4}}' > {output.inferred_tracts} 2>> {log}
        gita eval --true-tracts {input.true_tracts} --inferred-tracts {output.inferred_tracts} --output {output.performance} 2>> {log}
        grep -w Summary {output.performance} | awk -v train_demog={wildcards.train_demog} -v test_demog={wildcards.test_demog} \
                                                  -v geno_state={wildcards.geno_state} -v nref={wildcards.nref} -v ntgt={wildcards.ntgt} \
                                                  -v train_seed={wildcards.train_seed} -v test_seed={wildcards.test_seed} \
                                                  -v feature={params.feature_id} -v cutoff={wildcards.cutoff} \
                                                  'BEGIN{{OFS="\\t"}}{{print train_demog,test_demog,geno_state,nref,ntgt,train_seed,test_seed,feature,cutoff,$2,$3,$4,$5,$6}}' > {output.summary} 2>> {log}
        """ 


rule summary:
    input:
        summaries = expand("results/performance/train_{train_demog}_test_{test_demog}/nref_{nref}/ntgt_{ntgt}/{geno_state}/train_{train_seed}_test_{test_seed}/{output_prefix}.lr.{cutoff}.performance.summary",
                           train_demog=train_demog, test_demog=test_demog, nref=nref, ntgt=ntgt, geno_state=geno_state_list, output_prefix=output_prefix, train_seed=seed_list["train"], test_seed=seed_list["test"], cutoff=cutoff_list),
    output:
        summary = "results/performance/{output_prefix}.performance.summary",
    log:
        "logs/summary/{output_prefix}.log",
    resources:
        partition = "himem,gpu",
        mem_gb = 16,
    shell:
        """
        cat {input.summaries} | sed '1iTraining Data Demography\\tTest Data Demography\\tGenotype\\tNref\\tNtgt\\tTraining Seed\\tTest Seed\\tFeature\\tCutoff\\tPrecision\\tRecall\\tTrue_tracts_length\\tInferred_tracts_length\\tOverlaps_length' > {output.summary} 2>> {log}
        """
