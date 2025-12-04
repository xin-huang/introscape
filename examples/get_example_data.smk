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


ruleorder: extract_archaic_biallelic_snps > extract_1KG_biallelic_snps

example_pops = ["YRI", "CHS"]

# Separate rule for just creating examples
rule create_examples:
    input:
        expand("examples/data/Human/{modern_pop}/chr21.biallelic.snps.vcf.gz", modern_pop=example_pops),
        expand("examples/data/Human/{archaic_pop}/chr21.biallelic.snps.vcf.gz", archaic_pop=["NeaAltai", "Denisova"]),
        expand("examples/data/Human/ancestral_alleles/homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_chr21.bed.gz"),
        "examples/data/Human/merged/chr21.biallelic.snps.vcf.gz",
        "examples/data/Human/metadata/example_metadata.txt",


rule download_1KG_genomes:
    output:
        vcf="examples/data/Human/1KG/full_chr21.vcf.gz",
        idx="examples/data/Human/1KG/full_chr21.vcf.gz.tbi",
    shell:
        """
        wget -c https://ftp.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -O {output.vcf}
        wget -c https://ftp.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi -O {output.idx}
        """


rule download_1KG_info:
    output:
        samples="examples/data/Human/1KG/integrated_call_samples_v3.20130502.ALL.panel",
    shell:
        """
        wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel -O {output.samples}
        """


rule download_nea_genome:
    input:
    output:
        vcf = "examples/data/Human/NeaAltai/chr21_mq25_mapab100.vcf.gz",
        idx = "examples/data/Human/NeaAltai/chr21_mq25_mapab100.vcf.gz.tbi",
    shell:
        """
        wget -c http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Altai/chr21_mq25_mapab100.vcf.gz -O {output.vcf}
        wget -c http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Altai/chr21_mq25_mapab100.vcf.gz.tbi -O {output.idx}
        """


rule download_den_genome:
    input:
    output:
        vcf = "examples/data/Human/Denisova/chr21_mq25_mapab100.vcf.gz",
        idx = "examples/data/Human/Denisova/chr21_mq25_mapab100.vcf.gz.tbi",
    shell:
        """
        wget -c http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Denisova/chr21_mq25_mapab100.vcf.gz -O {output.vcf}
        wget -c http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Denisova/chr21_mq25_mapab100.vcf.gz.tbi -O {output.idx}
        """


rule create_example_metadata:
    input:
        samples=rules.download_1KG_info.output.samples,
        nea=rules.download_nea_genome.output.vcf,
        den=rules.download_den_genome.output.vcf,
    output:
        metadata="examples/data/Human/metadata/example_metadata.txt",
    params:
        pop1=f"{example_pops[0]}",
        pop2=f"{example_pops[1]}",
    shell:
        """
        grep -w {params.pop1} {input.samples} | awk '{{print $1}}' | head -5 | awk -v pop={params.pop1} '{{print pop"\\t"$1}}' > {output.metadata}
        grep -w {params.pop2} {input.samples} | awk '{{print $1}}' | head -5 | awk -v pop={params.pop2} '{{print pop"\\t"$1}}' >> {output.metadata}
        bcftools query -l {input.nea} | awk '{{print "NeaAltai\\t"$1}}' >> {output.metadata}
        bcftools query -l {input.den} | awk '{{print "Denisova\\t"$1}}' >> {output.metadata}
        sed -i '1iPopulation\\tSample' {output.metadata}
        """


rule download_ensembl_ancestral_alleles:
    output:
        anc_alleles="examples/data/Human/ancestral_alleles/homo_sapiens_ancestor_GRCh37_e71.tar.bz2",
    shell:
        """
        wget -c https://ftp.ensembl.org/pub/release-75/fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh37_e71.tar.bz2 -O {output.anc_alleles}
        tar -xvjf examples/data/Human/ancestral_alleles/homo_sapiens_ancestor_GRCh37_e71.tar.bz2 -C examples/data/Human/ancestral_alleles
        """


rule extract_anc_info:
    input:
        anc_alleles=rules.download_ensembl_ancestral_alleles.output.anc_alleles,
    output:
        bed=temp("examples/data/Human/ancestral_alleles/homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_chr21.bed"),
    params:
        fasta="examples/data/Human/ancestral_alleles/homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_21.fa",
    run:
        import pysam
        import re
   
        fasta = pysam.FastaFile(params.fasta)

        with open(output.bed, "wt") as out:
            for raw_chrom in fasta.references:
                match = re.search(
                    r"GRCh\d+:(\d+|X|Y)", raw_chrom
                )  # Extract chromosome number
                if not match:
                    print(f"Skipping unrecognized chromosome name: {raw_chrom}")
                    continue
                chrom = match.group(1)  # Keep only the number (no "chr" prefix)

                seq = fasta.fetch(raw_chrom).upper()
                for pos, base in enumerate(seq):
                    if base in "ACGT":  # Filter out 'N'
                        out.write(f"chr{chrom}\t{pos}\t{pos+1}\t{base}\n")

        fasta.close()


rule compress_anc_info:
    input:
        bed=rules.extract_anc_info.output.bed,
    output:
        bed="examples/data/Human/ancestral_alleles/homo_sapiens_ancestor_GRCh37_e71/homo_sapiens_ancestor_chr21.bed.gz",
    shell:
        """
        bgzip -c {input.bed} > {output.bed}
        tabix -p bed {output.bed}
        """


rule extract_1KG_biallelic_snps:
    input:
        vcf=rules.download_1KG_genomes.output.vcf,
        metadata=rules.create_example_metadata.output.metadata,
    output:
        vcf="examples/data/Human/{modern_pop}/chr21.biallelic.snps.vcf.gz",
        idx="examples/data/Human/{modern_pop}/chr21.biallelic.snps.vcf.gz.tbi",
    shell:
        """
        bcftools view {input.vcf} -S <(grep -w {wildcards.modern_pop} {input.metadata} | awk '{{print $2}}') | \
        bcftools view -v snps -m 2 -M 2 -g ^miss | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule extract_archaic_biallelic_snps:
    input:
        vcf="examples/data/Human/{archaic_pop}/chr21_mq25_mapab100.vcf.gz",
        metadata=rules.create_example_metadata.output.metadata,
    output:
        vcf="examples/data/Human/{archaic_pop}/chr21.biallelic.snps.vcf.gz",
        idx="examples/data/Human/{archaic_pop}/chr21.biallelic.snps.vcf.gz.tbi",
    shell:
        """
        bcftools view {input.vcf} -S <(grep -w {wildcards.archaic_pop} {input.metadata} | awk '{{print $2}}') | \
        bcftools view -v snps -m 2 -M 2 -g ^miss | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule merge_vcf:
    input:
        yri="examples/data/Human/YRI/chr21.biallelic.snps.vcf.gz",
        chs="examples/data/Human/CHS/chr21.biallelic.snps.vcf.gz",
        nea="examples/data/Human/NeaAltai/chr21.biallelic.snps.vcf.gz",
        den="examples/data/Human/Denisova/chr21.biallelic.snps.vcf.gz",
    output:
        vcf="examples/data/Human/merged/chr21.biallelic.snps.vcf.gz",
        idx="examples/data/Human/merged/chr21.biallelic.snps.vcf.gz.tbi",
    shell:
        """
        bcftools merge {input.yri} {input.chs} {input.nea} {input.den} | bcftools view -v snps -m 2 -M 2 -g ^miss | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """
