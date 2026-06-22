# Copyright 2025 Xin Huang
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


rule analyze_2src_samples:
    input:
        vcf = rules.add_2src_1KG.output.vcf,
        config = "config/analysis/1KG_w_{w}_y_{y}_z_{z}.yaml",
        anc_alleles = rules.fa2bed.output.bed,
    output:
        scores = "results/sai/2src/1KG/nea_den/w_{w}_y_{y}_z_{z}/{allele_type}/1KG.nea_den.{approach}.chr{i}.w_{w}_y_{y}_z_{z}.{allele_type}.scores.tsv",
        u_logs = "results/sai/2src/1KG/nea_den/w_{w}_y_{y}_z_{z}/{allele_type}/1KG.nea_den.{approach}.chr{i}.w_{w}_y_{y}_z_{z}.{allele_type}.scores.U.log",
        q_logs = "results/sai/2src/1KG/nea_den/w_{w}_y_{y}_z_{z}/{allele_type}/1KG.nea_den.{approach}.chr{i}.w_{w}_y_{y}_z_{z}.{allele_type}.scores.Q.log",
    params:
        chr_name = "{i}",
        win_len = 40000,
        win_step = 40000,
        anc_alleles = lambda wildcards: "" if wildcards.allele_type != "derived" else f"--anc-alleles results/polarized_data/anc_info/hg19.chr{wildcards.i}.anc.alleles.bed",
    resources:
        mem_gb = 128,
    shell:
        """
        sai score --vcf {input.vcf} --chr-name {params.chr_name} --output {output.scores} --win-len {params.win_len} --win-step {params.win_step} --config {input.config} {params.anc_alleles}
        """


rule get_2src_samples_outliers:
    input:
        scores = expand("results/sai/2src/1KG/nea_den/w_{w}_y_{y}_z_{z}/{allele_type}/1KG.nea_den.{approach}.chr{i}.w_{w}_y_{y}_z_{z}.{allele_type}.scores.tsv", i=list(range(1,23)), allow_missing=True),
    output:
        all_scores = "results/sai/2src/1KG/nea_den/w_{w}_y_{y}_z_{z}/{allele_type}/1KG.nea_den.{approach}.w_{w}_y_{y}_z_{z}.{allele_type}.scores.tsv",
        u_outliers = "results/sai/2src/1KG/nea_den/w_{w}_y_{y}_z_{z}/{allele_type}/1KG.nea_den.{approach}.w_{w}_y_{y}_z_{z}.{allele_type}.scores.U.0.999.outliers.tsv",
        q_outliers = "results/sai/2src/1KG/nea_den/w_{w}_y_{y}_z_{z}/{allele_type}/1KG.nea_den.{approach}.w_{w}_y_{y}_z_{z}.{allele_type}.scores.Q.0.999.outliers.tsv",
    params:
        outlier_quantile = 0.999,
        output_prefix = "results/sai/2src/1KG/nea_den/w_{w}_y_{y}_z_{z}/{allele_type}/1KG.nea_den.{approach}.w_{w}_y_{y}_z_{z}.{allele_type}.scores",
    shell:
        """
        set +o pipefail
        cat {input.scores} | head -1 > {output.all_scores}
        cat {input.scores} | grep -v "Chrom" >> {output.all_scores}
        sai outlier --score {output.all_scores} --output-prefix {params.output_prefix} --quantile {params.outlier_quantile}
        """
