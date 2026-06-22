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


rule download_nea_genome:
    input:
    output:
        vcf = "resources/data/NeaAltai/not_filtered/chr{i}_mq25_mapab100.vcf.gz",
        index = "resources/data/NeaAltai/not_filtered/chr{i}_mq25_mapab100.vcf.gz.tbi",
    shell:
        """
        wget -c http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Altai/chr{wildcards.i}_mq25_mapab100.vcf.gz -O {output.vcf}
        wget -c http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Altai/chr{wildcards.i}_mq25_mapab100.vcf.gz.tbi -O {output.index}
        """


rule download_den_genome:
    input:
    output:
        vcf = "resources/data/Denisova/not_filtered/chr{i}_mq25_mapab100.vcf.gz",
        index = "resources/data/Denisova/not_filtered/chr{i}_mq25_mapab100.vcf.gz.tbi",
    shell:
        """
        wget -c http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Denisova/chr{wildcards.i}_mq25_mapab100.vcf.gz -O {output.vcf}
        wget -c http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Denisova/chr{wildcards.i}_mq25_mapab100.vcf.gz.tbi -O {output.index}
        """


rule download_nea_filter:
    input:
    output:
        bed = "resources/data/NeaAltai/FilterBed/chr{i}_mask.bed.gz",
    shell:
        """
        wget -c http://cdna.eva.mpg.de/neandertal/Vindija/FilterBed/Altai/chr{wildcards.i}_mask.bed.gz -O {output.bed}
        """


rule download_den_filter:
    input:
    output:
        bed = "resources/data/Denisova/FilterBed/chr{i}_mask.bed.gz",
    shell:
        """
        wget -c http://cdna.eva.mpg.de/neandertal/Vindija/FilterBed/Denisova/chr{wildcards.i}_mask.bed.gz -O {output.bed}
        """


rule filter_nea_genome:
    input:
        vcf = rules.download_nea_genome.output.vcf,
        bed = rules.download_nea_filter.output.bed,
    output:
        vcf = "resources/data/NeaAltai/filtered/chr{i}_mq25_mapab100.vcf.gz",
        index = "resources/data/NeaAltai/filtered/chr{i}_mq25_mapab100.vcf.gz.tbi",
    shell:
        """
        bcftools view {input.vcf} -R {input.bed} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule filter_den_genome:
    input:
        vcf = rules.download_den_genome.output.vcf,
        bed = rules.download_den_filter.output.bed,
    output:
        vcf = "resources/data/Denisova/filtered/chr{i}_mq25_mapab100.vcf.gz",
        index = "resources/data/Denisova/filtered/chr{i}_mq25_mapab100.vcf.gz.tbi",
    shell:
        """
        bcftools view {input.vcf} -R {input.bed} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule download_1KG_genomes:
    input:
    output:
        vcf = "resources/data/1KG/ALL.chr{i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz",
        index = "resources/data/1KG/ALL.chr{i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi",
    shell:
        """
        wget -c https://ftp.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr{wildcards.i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -O {output.vcf}
        wget -c https://ftp.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr{wildcards.i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi -O {output.index}
        """


rule download_1KG_info:
    input:
    output:
        ref = "resources/data/1KG/human_g1k_v37.fasta",
        samples = "resources/data/1KG/integrated_call_samples_v3.20130502.ALL.panel",
    shell:
        """
        wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel -O {output.samples}
        wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
        gzip -d human_g1k_v37.fasta.gz || true
        mv human_g1k_v37.fasta {output.ref}
        """


rule download_pan_genomes:
    input:
    output:
        vcf = "resources/data/pan/pan.chr{i}.filteranno.vcf.gz",
        index = "resources/data/pan/pan.chr{i}.filteranno.vcf.gz.csi",
    shell:
        """
        wget -c https://phaidra.univie.ac.at/pfsa/o_2066302/merged_segregating/Pan/Pan_wild_filtered/chr{wildcards.i}.filteranno.vcf.gz -O {output.vcf}
        wget -c https://phaidra.univie.ac.at/pfsa/o_2066302/merged_segregating/Pan/Pan_wild_filtered/chr{wildcards.i}.filteranno.vcf.gz.csi -O {output.index}
        """


rule download_pan_metadata:
    input:
    output:
        txt = "resources/data/pan/metadata_full.txt",
    shell:
        """
        wget -c https://phaidra.univie.ac.at/pfsa/o_2066302/metadata_full.txt -O {output.txt}
        """


rule download_hg19_anc_info:
    input:
    output:
        expand("resources/data/anc_info/hg19/homo_sapiens_ancestor_{i}.fa", i=np.arange(1,23)),
    params:
        hg19_dir = "resources/data/anc_info/hg19",
    shell:
        """
        wget -c https://ftp.ensembl.org/pub/release-75/fasta/ancestral_alleles/homo_sapiens_ancestor_GRCh37_e71.tar.bz2
        mkdir -p {params.hg19_dir}
        tar -xjvf homo_sapiens_ancestor_GRCh37_e71.tar.bz2 -C {params.hg19_dir} --strip-components=1
        rm homo_sapiens_ancestor_GRCh37_e71.tar.bz2
        """


rule download_refgenomes:
    input:
    output:
        hg38 = "resources/data/refgenomes/hg38.fa",
        rheMac10 = "resources/data/refgenomes/rheMac10.fa",
        chain = "resources/data/refgenomes/hg38ToRheMac10.over.chain",
    params:
        dir = "resources/data/refgenomes",
    shell:
        """
        wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -O {params.dir}/hg38.fa.gz
        gzip -d {params.dir}/hg38.fa.gz
        wget -c https://hgdownload.soe.ucsc.edu/goldenPath/rheMac10/bigZips/rheMac10.fa.gz -O {params.dir}/rheMac10.fa.gz
        gzip -d {params.dir}/rheMac10.fa.gz
        wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToRheMac10.over.chain.gz -O {params.dir}/hg38ToRheMac10.over.chain.gz
        gzip -d {params.dir}/hg38ToRheMac10.over.chain.gz
        """


rule download_beagle:
    input:
    output:
        beagle = "resources/tools/beagle/beagle.06Aug24.a91.jar",
    shell:
        """
        wget -c https://faculty.washington.edu/browning/beagle/beagle.06Aug24.a91.jar
        mv beagle.06Aug24.a91.jar {output.beagle}
        """


rule download_annovar_db:
    input:
    output:
        hg19_avsnp150 = "resources/tools/annovar/humandb/hg19_avsnp150.txt",
        hg19_dbnsfp42c = "resources/tools/annovar/humandb/hg19_dbnsfp42c.txt",
        hg38_refgene = "resources/tools/annovar/humandb/hg38_refGene.txt",
        hg38_avsnp150 = "resources/tools/annovar/humandb/hg38_avsnp150.txt",
        hg38_dbnsfp42c = "resources/tools/annovar/humandb/hg38_dbnsfp42c.txt",
    shell:
        """
        resources/tools/annovar/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar avsnp150 resources/tools/annovar/humandb/
        resources/tools/annovar/annotate_variation.pl -downdb -buildver hg19 -webfrom annovar dbnsfp42c resources/tools/annovar/humandb/
        resources/tools/annovar/annotate_variation.pl -downdb -buildver hg38 -webfrom annovar refGene resources/tools/annovar/humandb/
        resources/tools/annovar/annotate_variation.pl -downdb -buildver hg38 -webfrom annovar avsnp150 resources/tools/annovar/humandb/
        resources/tools/annovar/annotate_variation.pl -downdb -buildver hg38 -webfrom annovar dbnsfp42c resources/tools/annovar/humandb/
        """


rule fa2bed:
    input:
        fa = "resources/data/anc_info/hg19/homo_sapiens_ancestor_{i}.fa",
    output:
        bed = "results/polarized_data/anc_info/hg19.chr{i}.anc.alleles.bed",
    script:
        "../scripts/fa2bed.py"
