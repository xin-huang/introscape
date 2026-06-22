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


import numpy as np


rule annotate_pan_candidates:
    input:
        outliers = "results/sai/1src/pan/PPA/w_{w}_y_{y}/{allele_type}/pan.PPA.w_{w}_y_{y}.{allele_type}.scores.{stat}.0.999.outliers.tsv",
        annotation = expand("results/annotated_data/pan/pan.chr{i}.biallelic.snps.hg38_multianno.txt", i=list(range(1,23)), allow_missing=True),
    output:
        annotated_candidates = "results/sai/1src/pan/PPA/w_{w}_y_{y}/{allele_type}/pan.PPA.w_{w}_y_{y}.{allele_type}.scores.{stat}.0.999.outliers.annotated.tsv",
    resources:
        mem_gb = 128,
    script:
       "../scripts/get_annotated_candidates.py" 


rule get_1KG_overlapping_candidate_windows:
    input:
        u_outliers = "results/sai/2src/1KG/nea_den/w_{w}_y_{y}_z_{z}/{allele_type}/1KG.nea_den.{approach}.w_{w}_y_{y}_z_{z}.{allele_type}.scores.U.0.999.outliers.tsv",
        q_outliers = "results/sai/2src/1KG/nea_den/w_{w}_y_{y}_z_{z}/{allele_type}/1KG.nea_den.{approach}.w_{w}_y_{y}_z_{z}.{allele_type}.scores.Q.0.999.outliers.tsv",
    output:
        overlapping = "results/sai/2src/1KG/nea_den/w_{w}_y_{y}_z_{z}/{allele_type}/1KG.nea_den.{approach}.w_{w}_y_{y}_z_{z}.{allele_type}.scores.U.Q.0.999.overlapping.outliers.tsv",
    script:
        "../scripts/get_overlapping_windows.py"


rule get_1KG_overlapping_candidate_SNPs:
    input:
        windows = rules.get_1KG_overlapping_candidate_windows.output.overlapping,
        U_snps = expand(
            "results/sai/2src/1KG/nea_den/w_{w}_y_{y}_z_{z}/{allele_type}/1KG.nea_den.{approach}.chr{i}.w_{w}_y_{y}_z_{z}.{allele_type}.scores.U.log",
            i=np.arange(1,23),
            allow_missing=True, 
        ),
        Q_snps = expand(
            "results/sai/2src/1KG/nea_den/w_{w}_y_{y}_z_{z}/{allele_type}/1KG.nea_den.{approach}.chr{i}.w_{w}_y_{y}_z_{z}.{allele_type}.scores.Q.log",
            i=np.arange(1,23),
            allow_missing=True,
        ),
    output:
        overlapping = "results/sai/2src/1KG/nea_den/w_{w}_y_{y}_z_{z}/{allele_type}/1KG.nea_den.{approach}.w_{w}_y_{y}_z_{z}.{allele_type}.scores.U.Q.0.999.overlapping.outliers.snps.tsv",
    script:
        "../scripts/get_overlapping_snps.py"


rule get_1KG_overlapping_candidate_SNPs_genes:
    input:
        candidate_snps = rules.get_1KG_overlapping_candidate_SNPs.output.overlapping,
        annotation = expand("results/annotated_data/1KG/1KG.nea_den.{approach}.chr{i}.hg19_multianno.txt", i=list(range(1,23)), allow_missing=True),
    output:
        candidate_genes = "results/sai/2src/1KG/nea_den/w_{w}_y_{y}_z_{z}/{allele_type}/1KG.nea_den.{approach}.w_{w}_y_{y}_z_{z}.{allele_type}.scores.{stat}.0.999.overlapping.outliers.snps.genes.tsv",
    shell:
        """
        sed '1d' {input.candidate_snps} | awk '{{print $4}}' | awk -F ":" '{{print $1"\\t"$2}}' | while IFS=$'\\t' read -r c p; do awk -v pos=${{p}} '$2==pos' results/annotated_data/1KG/1KG.nea_den.{wildcards.approach}.chr${{c}}.hg19_multianno.txt; done > {output.candidate_genes}
        """


rule get_pan_candidate_genes:
    input:
        annotated_candidates = rules.annotate_pan_candidates.output.annotated_candidates,
    output:
        candidate_genes = "results/sai/1src/pan/PPA/w_{w}_y_{y}/{allele_type}/pan.PPA.w_{w}_y_{y}.{allele_type}.scores.{stat}.0.999.outliers.candidate.genes.tsv",
    shell:
        """
        sed '1d' {input.annotated_candidates} | awk '{{print $7}}' | grep -v ";" | sort | uniq > {output.candidate_genes}
        """


rule get_pan_polarized_data_candidate_genes:
    input:
        candidate_genes = expand(
            "results/sai/1src/pan/PPA/w_{w}_y_{y}/derived/pan.PPA.w_{w}_y_{y}.derived.scores.{stat}.0.999.outliers.candidate.genes.tsv",
            stat=["DD", "Danc", "Dplus", "df", "fd", "U", "Q"],
            allow_missing=True,
        ),
    output:
        overlapping_genes = "results/sai/1src/pan/PPA/w_{w}_y_{y}/derived/pan.PPA.w_{w}_y_{y}.derived.scores.0.999.outliers.overlapping.candidate.genes.tsv",
    shell:
        """
        cat {input.candidate_genes} | sort | uniq -c | sort -rnk 1,1 | awk '$1>1' > {output.overlapping_genes}
        """


rule get_pan_unpolarized_data_candidate_genes:
    input:
        candidate_genes = expand(
            "results/sai/1src/pan/PPA/w_{w}_y_{y}/any/pan.PPA.w_{w}_y_{y}.any.scores.{stat}.0.999.outliers.candidate.genes.tsv",
            stat=["DD", "U", "Q"],
            allow_missing=True,
        ),
    output:
        overlapping_genes = "results/sai/1src/pan/PPA/w_{w}_y_{y}/any/pan.PPA.w_{w}_y_{y}.any.scores.0.999.outliers.overlapping.candidate.genes.tsv",
    shell:
        """
        cat {input.candidate_genes} | sort | uniq -c | sort -rnk 1,1 | awk '$1>1' > {output.overlapping_genes}
        """
