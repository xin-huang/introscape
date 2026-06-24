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


rule extract_biallelic_snps:
    input:
        vcf=get_vcf_input_path,
    output:
        vcf="results/processed_data/{dataset}/all/{chrom}.biallelic.snps.vcf.gz",
        idx="results/processed_data/{dataset}/all/{chrom}.biallelic.snps.vcf.gz.tbi",
    log:
        "logs/preprocess/extract_biallelic_snps.{dataset}.{chrom}.log",
    conda:
        "../envs/introscape.yaml"
    shell:
        """
        mkdir -p results/processed_data/{wildcards.dataset}/all logs/preprocess
        ( bcftools view {input.vcf} -v snps -m 2 -M 2 -g ^miss | bgzip -c > {output.vcf} ) 2> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """


rule create_analysis_sample_lists:
    input:
        metadata=get_metadata,
    output:
        ref="results/samples/{analysis}/ref.samples.list",
        tgt="results/samples/{analysis}/tgt.samples.list",
        src="results/samples/{analysis}/src.samples.list",
    params:
        ref_population=get_ref_population,
        tgt_population=get_tgt_population,
        src_population=get_src_population,
    log:
        "logs/preprocess/create_analysis_sample_lists.{analysis}.log",
    conda:
        "../envs/introscape.yaml"
    script:
        "../scripts/make_sample_lists.py"


rule merge_analysis_samples:
    input:
        ref=rules.create_analysis_sample_lists.output.ref,
        tgt=rules.create_analysis_sample_lists.output.tgt,
    output:
        samples="results/samples/{analysis}/analysis.samples.list",
    log:
        "logs/preprocess/merge_analysis_samples.{analysis}.log",
    conda:
        "../envs/introscape.yaml"
    shell:
        """
        mkdir -p results/samples/{wildcards.analysis} logs/preprocess
        cat {input.ref} {input.tgt} > {output.samples} 2> {log}
        """


rule extract_analysis_data:
    input:
        vcf=get_biallelic_vcf_for_analysis,
        idx=get_biallelic_index_for_analysis,
        samples=rules.merge_analysis_samples.output.samples,
    output:
        vcf="results/processed_data/{analysis}/chromosomes/{chrom}.biallelic.snps.vcf.gz",
        idx="results/processed_data/{analysis}/chromosomes/{chrom}.biallelic.snps.vcf.gz.tbi",
    log:
        "logs/preprocess/extract_analysis_data.{analysis}.{chrom}.log",
    conda:
        "../envs/introscape.yaml"
    shell:
        """
        mkdir -p results/processed_data/{wildcards.analysis}/chromosomes logs/preprocess
        ( bcftools view {input.vcf} -S {input.samples} --force-samples | \
          bcftools view -i "AC>0 && AC<AN" | \
          bgzip -c > {output.vcf} ) 2> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """


rule concat_analysis_vcf:
    input:
        vcfs=get_analysis_chrom_vcfs,
        idxs=get_analysis_chrom_indexes,
    output:
        vcf="results/processed_data/{analysis}/{analysis}.biallelic.snps.vcf.gz",
        idx="results/processed_data/{analysis}/{analysis}.biallelic.snps.vcf.gz.tbi",
    log:
        "logs/preprocess/concat_analysis_vcf.{analysis}.log",
    conda:
        "../envs/introscape.yaml"
    shell:
        """
        mkdir -p results/processed_data/{wildcards.analysis} logs/preprocess
        ( bcftools concat {input.vcfs} | bgzip -c > {output.vcf} ) 2> {log}
        tabix -p vcf {output.vcf} 2>> {log}
        """
