rule process_test_data:
    input:
        vcf = rules.compress_vcf.output.vcf,
        ref = rules.simulate_test_data.output.ref,
        tgt = rules.simulate_test_data.output.tgt,
        feature_config = feature_config,
    output:
        features = "results/data/test/{test_demog}/" +
                   "nref_{nref}/ntgt_{ntgt}/{test_seed}/" +
                   "{output_prefix}.{geno_state}.features.gz",
    params:
        seq_len = seq_len["train"],
        win_step = win_step,
        is_phased = lambda wildcards: '--phased' if wildcards.geno_state == 'phased' else '',
        output_dir = os.path.join(output_dir["test"], "{test_seed}"),
        features = "results/data/test/{test_demog}/" +
                   "nref_{nref}/ntgt_{ntgt}/{test_seed}/" +
                   "{output_prefix}.{geno_state}.features",
        chr_name = 1,
        nprocess = 1000,
    log:
        "logs/process_test_data/{test_demog}/" +
        "nref_{nref}/ntgt_{ntgt}/{test_seed}/" +
        "{output_prefix}.{geno_state}.log",
    resources:
        mem_gb = 128,
        cpus = 16,
        time = 360,
    shell:
        """
        gaia lr preprocess --vcf {input.vcf} --ref {input.ref} --tgt {input.tgt} --feature-config {input.feature_config} \
                           --win-len {params.seq_len} --win-step {params.win_step} --output-prefix {wildcards.output_prefix}.{wildcards.geno_state} \
                           --output-dir {params.output_dir} --nprocess {params.nprocess} {params.is_phased} --chr-name {params.chr_name}
        gzip -c {params.features} > {output.features} 2> {log}
        rm {params.features} 2>> {log}
        """


rule test_logistic_regression_model:
    input:
        features = rules.process_test_data.output.features,
        model_file = rules.train_logistic_regression_model.output.model_file,
    output:
        pred = "results/performance/train_{train_demog}_test_{test_demog}/" +
               "nref_{nref}/ntgt_{ntgt}/{geno_state}/{dataset_state}/" +
               "train_{train_seed}_test_{test_seed}/{output_prefix}.lr.pred.gz",
    params:
        pred = "results/performance/train_{train_demog}_test_{test_demog}/" +
               "nref_{nref}/ntgt_{ntgt}/{geno_state}/{dataset_state}/" +
               "train_{train_seed}_test_{test_seed}/{output_prefix}.lr.pred",
    log:
        "logs/test_logistic_regression_model/" +
        "train_{train_demog}_test_{test_demog}/" +
        "nref_{nref}/ntgt_{ntgt}/{geno_state}/" +
        "{dataset_state}/train_{train_seed}_test_{test_seed}/" +
        "{output_prefix}.log",
    resources:
        partition = "himem",
        mem_gb = 128,
    shell:
        """
        gaia lr infer --inference-data {input.features} --model-file {input.model_file} --output-file {params.pred}
        gzip -c {params.pred} > {output.pred} 2> {log}
        rm {params.pred} 2>> {log}
        """


rule evaluate_logistic_regression_model:
    input:
        pred = rules.test_logistic_regression_model.output.pred,
        true_tracts = lambda wildcards: (
            rules.get_phased_true_tracts.output.bed 
            if wildcards.geno_state == 'phased' 
            else rules.get_unphased_true_tracts.output.bed
        ),
    output:
        inferred_tracts = "results/performance/train_{train_demog}_test_{test_demog}/" +
                          "nref_{nref}/ntgt_{ntgt}/{geno_state}/{dataset_state}/" +
                          "train_{train_seed}_test_{test_seed}/" +
                          "{output_prefix}.lr.{cutoff}.inferred.tracts.bed",
        performance = "results/performance/train_{train_demog}_test_{test_demog}/" +
                      "nref_{nref}/ntgt_{ntgt}/{geno_state}/{dataset_state}/" +
                      "train_{train_seed}_test_{test_seed}/" +
                      "{output_prefix}.lr.{cutoff}.performance",
        summary = "results/performance/train_{train_demog}_test_{test_demog}/" +
                  "nref_{nref}/ntgt_{ntgt}/{geno_state}/{dataset_state}/" +
                  "train_{train_seed}_test_{test_seed}/" +
                  "{output_prefix}.lr.{cutoff}.performance.summary",
    params:
        feature_id = feature_id,
        ntgt = ntgt,
        seq_len = seq_len["test"],
        ploidy = ploidy,
    log:
        "logs/evaluate_logistic_regression_model/train_{train_demog}_test_{test_demog}/" +
        "nref_{nref}/ntgt_{ntgt}/{geno_state}/{dataset_state}/" +
        "train_{train_seed}_test_{test_seed}/{output_prefix}.{cutoff}.log",
    resources:
        partition = "himem,gpu",
    shell:
        """
        zcat {input.pred} | sed '1d' | awk '$NF>{wildcards.cutoff}' | awk 'BEGIN{{OFS="\\t"}}{{print $1,$2,$3,$4}}' > {output.inferred_tracts} 2>> {log}
        gaia eval --true-tracts {input.true_tracts} --inferred-tracts {output.inferred_tracts} --seq-len {params.seq_len} \
                  --sample-size {params.ntgt} --ploidy {params.ploidy} --output {output.performance} 2>> {log}
        grep -w Average {output.performance} | awk -v train_demog={wildcards.train_demog} -v test_demog={wildcards.test_demog} \
                                                  -v geno_state={wildcards.geno_state} -v dataset_state={wildcards.dataset_state} -v nref={wildcards.nref} -v ntgt={wildcards.ntgt} \
                                                  -v train_seed={wildcards.train_seed} -v test_seed={wildcards.test_seed} \
                                                  -v feature={params.feature_id} -v cutoff={wildcards.cutoff} \
                                                  'BEGIN{{OFS="\\t"}}{{print train_demog,test_demog,geno_state,dataset_state,nref,ntgt,train_seed,test_seed,feature,cutoff,$2,$3,$4,$5,$6,$7,$8}}' > {output.summary} 2>> {log}
        """ 


rule summary:
    input:
        summaries = expand(
            "results/performance/train_{train_demog}_test_{test_demog}/" +
            "nref_{nref}/ntgt_{ntgt}/{geno_state}/{dataset_state}/" +
            "train_{train_seed}_test_{test_seed}/" +
            "{output_prefix}.lr.{cutoff}.performance.summary",
            train_demog=train_demog, 
            test_demog=test_demog, 
            nref=nref, 
            ntgt=ntgt, 
            geno_state=geno_state_list, 
            dataset_state=dataset_state_list, 
            output_prefix=output_prefix, 
            train_seed=seed_list["train"], 
            test_seed=seed_list["test"], 
            cutoff=cutoff_list
        ),
    output:
        summary = "results/performance/{output_prefix}.performance.summary",
    log:
        "logs/summary/{output_prefix}.log",
    resources:
        partition = "himem,gpu",
        mem_gb = 32,
    shell:
        """
        cat {input.summaries} | sed '1iTraining Data Demography\\tTest Data Demography\\tGenotype\\tData balance\\tNref\\tNtgt\\tTraining Seed\\tTest Seed\\tFeature\\tCutoff\\tSequence_length\\tTrue_tracts_length\\tInferred_tracts_length\\tTrue_positives_length\\tFalse_positives_length\\tTrue_negatives_length\\tFalse_negatives_length' > {output.summary} 2>> {log}
        """
