rule simulate_training_data:
    input:
        demes_file = demes["train"],
        features_config = feature_config,
    output:
        features = "results/data/train/{train_demog}/" +
                   "nref_{nref}/ntgt_{ntgt}/{geno_state}/" +
                   "{dataset_state}/{train_seed}/{output_prefix}.features",
    params:
        nref = nref,
        ntgt = ntgt,
        nfeature = nfeature,
        ref_id = ref_id["train"],
        tgt_id = tgt_id["train"],
        src_id = src_id["train"],
        mut_rate = mut_rate["train"],
        rec_rate = rec_rate["train"],
        seq_len = seq_len["train"], 
        is_phased = lambda wildcards: '--phased' if wildcards.geno_state == 'phased' else '',
        is_balanced = lambda wildcards: '--force-balanced' if wildcards.dataset_state == 'balanced' else '',
        output_dir = os.path.join(output_dir["train"], "{geno_state}/{dataset_state}/{train_seed}"),
        nrep = 1000,
        nprocess = 1000,
    log:
        "logs/simulate_training_data/{train_demog}/" +
        "nref_{nref}/ntgt_{ntgt}/{geno_state}/" +
        "{dataset_state}/{train_seed}/{output_prefix}.log",
    resources:
        cpus = 16, 
        mem_gb = 512,
        time = 4320,
        partition = "himem",
    shell:
        """
        gaia lr simulate --demes {input.demes_file} --nref {params.nref} --ntgt {params.ntgt} \
                         --ref-id {params.ref_id} --tgt-id {params.tgt_id} --src-id {params.src_id} \
                         --mut-rate {params.mut_rate} --rec-rate {params.rec_rate} --seq-len {params.seq_len} \
                         --output-prefix {wildcards.output_prefix} --output-dir {params.output_dir} \
                         --seed {wildcards.train_seed} --replicate {params.nrep} --nprocess {params.nprocess} \
                         --feature-config {input.features_config} --nfeature {params.nfeature} \
                         {params.is_phased} {params.is_balanced} 2>> {log}
        """


rule train_logistic_regression_model:
    input:
        features = rules.simulate_training_data.output.features
    output:
        model_file = "results/data/train/{train_demog}/" +
                     "nref_{nref}/ntgt_{ntgt}/{geno_state}/" +
                     "{dataset_state}/{train_seed}/{output_prefix}.lr.model",
    log:
        "logs/train_logistic_regression_model/{train_demog}/" +
        "nref_{nref}/ntgt_{ntgt}/{geno_state}/" +
        "{dataset_state}/{train_seed}/{output_prefix}.log", 
    resources:
        partition = "himem,gpu", 
        mem_gb = 32,
    shell:
        """
        gaia lr train --training-data {input.features} --model-file {output.model_file} --seed {wildcards.train_seed} 2>> {log}
        """
