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


rule analyze_pan_samples:
    input:
        vcf = rules.extract_pan_biallelic_snps.output.vcf,
        config = "config/analysis/pan_w_{w}_y_{y}.yaml",
        anc_alleles = rules.extract_pan_anc_info.output.anc_alleles,
    output:
        scores = "results/sai/1src/pan/PPA/w_{w}_y_{y}/{allele_type}/pan.PPA.chr{i}.w_{w}_y_{y}.{allele_type}.scores.tsv",
    params:
        chr_name = "chr{i}",
        win_len = 50000,
        win_step = 50000,
        anc_alleles = lambda wildcards: "" if wildcards.allele_type != "derived" else f"--anc-alleles results/polarized_data/anc_info/hg38.chr{wildcards.i}.anc.alleles.bed",
    resources:
        mem_gb = 128,
    shell:
        """
        sai score --vcf {input.vcf} --chr-name {params.chr_name} --output {output.scores} --win-len {params.win_len} --win-step {params.win_step} --config {input.config} {params.anc_alleles}
        """


rule get_pan_outliers:
    input:
        scores = expand("results/sai/1src/pan/PPA/w_{w}_y_{y}/{allele_type}/pan.PPA.chr{i}.w_{w}_y_{y}.{allele_type}.scores.tsv", i=list(range(1,23)), allow_missing=True),
    output:
        all_scores = "results/sai/1src/pan/PPA/w_{w}_y_{y}/{allele_type}/pan.PPA.w_{w}_y_{y}.{allele_type}.scores.tsv",
        dd_outliers = "results/sai/1src/pan/PPA/w_{w}_y_{y}/{allele_type}/pan.PPA.w_{w}_y_{y}.{allele_type}.scores.DD.0.999.outliers.tsv",
        danc_outliers = "results/sai/1src/pan/PPA/w_{w}_y_{y}/{allele_type}/pan.PPA.w_{w}_y_{y}.{allele_type}.scores.Danc.0.999.outliers.tsv",
        dplus_outliers = "results/sai/1src/pan/PPA/w_{w}_y_{y}/{allele_type}/pan.PPA.w_{w}_y_{y}.{allele_type}.scores.Dplus.0.999.outliers.tsv",
        df_outliers = "results/sai/1src/pan/PPA/w_{w}_y_{y}/{allele_type}/pan.PPA.w_{w}_y_{y}.{allele_type}.scores.df.0.999.outliers.tsv",
        fd_outliers = "results/sai/1src/pan/PPA/w_{w}_y_{y}/{allele_type}/pan.PPA.w_{w}_y_{y}.{allele_type}.scores.fd.0.999.outliers.tsv",
        u_outliers = "results/sai/1src/pan/PPA/w_{w}_y_{y}/{allele_type}/pan.PPA.w_{w}_y_{y}.{allele_type}.scores.U.0.999.outliers.tsv",
        q_outliers = "results/sai/1src/pan/PPA/w_{w}_y_{y}/{allele_type}/pan.PPA.w_{w}_y_{y}.{allele_type}.scores.Q.0.999.outliers.tsv",
    params:
        outlier_quantile = 0.999,
        output_prefix = "results/sai/1src/pan/PPA/w_{w}_y_{y}/{allele_type}/pan.PPA.w_{w}_y_{y}.{allele_type}.scores",
    shell:
        """
        set +o pipefail
        cat {input.scores} | head -1 > {output.all_scores}
        cat {input.scores} | grep -v "Chrom" >> {output.all_scores}
        sai outlier --score {output.all_scores} --output-prefix {params.output_prefix} --quantile {params.outlier_quantile}
        """
