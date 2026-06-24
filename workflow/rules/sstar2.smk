# Copyright 2026 Xin Huang
#
# GNU General Public License v3.0
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, please see
#
#    https://www.gnu.org/licenses/gpl-3.0.en.html


rule render_sstar2_config:
    input:
        method_config=get_sstar2_method_config_path,
        vcf=rules.concat_analysis_vcf.output.vcf,
        ref=rules.create_analysis_sample_lists.output.ref,
        tgt=rules.create_analysis_sample_lists.output.tgt,
    output:
        config="results/sstar2/{analysis}/config/{analysis}.sstar2.config.yaml",
    log:
        "logs/sstar2/render_sstar2_config.{analysis}.log",
    conda:
        "../envs/introscape.yaml"
    script:
        "../scripts/render_sstar2_config.py"


rule train_sstar2_model:
    input:
        demes=get_analysis_demography,
        config=rules.render_sstar2_config.output.config,
    output:
        model="results/sstar2/{analysis}/model/{analysis}.sstar2.model.onnx",
    log:
        "logs/sstar2/train_sstar2_model.{analysis}.log",
    resources:
        mem_mb=64000,
        cpus=32,
        time=360,
    conda:
        "../envs/sstar2.yaml"
    shell:
        """
        mkdir -p results/sstar2/{wildcards.analysis}/model logs/sstar2
        sstar2 train --demes {input.demes} --config {input.config} --output {output.model} > {log} 2>&1
        """


rule infer_sstar2_tracts:
    input:
        model=rules.train_sstar2_model.output.model,
        config=rules.render_sstar2_config.output.config,
    output:
        features="results/sstar2/{analysis}/inference/{analysis}.sstar2.features.tsv",
        predictions="results/sstar2/{analysis}/inference/{analysis}.sstar2.pred.tsv",
        tracts="results/sstar2/{analysis}/inference/{analysis}.sstar2.inferred.tracts.bed",
    log:
        "logs/sstar2/infer_sstar2_tracts.{analysis}.log",
    resources:
        mem_mb=64000,
        cpus=32,
    conda:
        "../envs/sstar2.yaml"
    shell:
        """
        mkdir -p results/sstar2/{wildcards.analysis}/inference logs/sstar2
        sstar2 infer \
          --model {input.model} \
          --config {input.config} \
          --feat-file {output.features} \
          --pred-file {output.predictions} \
          --tract-file {output.tracts} > {log} 2>&1
        """


rule convert_sstar2_tracts_to_candidates:
    input:
        tracts=rules.infer_sstar2_tracts.output.tracts,
    output:
        candidates="results/sstar2/{analysis}/candidates/{analysis}.sstar2.candidate.regions.tsv",
    params:
        dataset=get_analysis_dataset_id,
    log:
        "logs/sstar2/convert_sstar2_tracts_to_candidates.{analysis}.log",
    conda:
        "../envs/introscape.yaml"
    script:
        "../scripts/sstar2_tracts_to_candidates.py"
