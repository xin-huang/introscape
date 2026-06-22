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


rule create_1KG_samples:
    input:
        samples = rules.download_1KG_info.output.samples,
        nea_vcf = "resources/data/NeaAltai/not_filtered/chr21_mq25_mapab100.vcf.gz",
        den_vcf = "resources/data/Denisova/not_filtered/chr21_mq25_mapab100.vcf.gz",
    output:
        ref = "results/processed_data/1KG/samples/1KG.ref.samples.txt",
        tgt = "results/processed_data/1KG/samples/1KG.tgt.samples.txt",
        nea = "results/processed_data/1KG/samples/1KG.nea.samples.txt",
        den = "results/processed_data/1KG/samples/1KG.den.samples.txt",
        nea_den = "results/processed_data/1KG/samples/1KG.nea_den.samples.txt",
    shell:
        """
        grep -w 'AFR' {input.samples} | grep -v ACB | grep -v ASW | awk '{{print "AFR+EAS\\t"$1}}' > {output.ref}
        grep -w 'EAS' {input.samples} | awk '{{print "AFR+EAS\\t"$1}}' >> {output.ref}
        grep -w 'EUR' {input.samples} | awk '{{print $3"\\t"$1}}' > {output.tgt}
        bcftools query -l {input.nea_vcf} | awk '{{print "NEA\\t"$0}}' > {output.nea}
        bcftools query -l {input.den_vcf} | awk '{{print "DEN\\t"$0}}' > {output.den}
        cat {output.nea} > {output.nea_den}
        cat {output.den} >> {output.nea_den}
        """


rule extract_1KG_biallelic_snps:
    input:
        vcf = rules.download_1KG_genomes.output.vcf,
        ref = rules.create_1KG_samples.output.ref,
        tgt = rules.create_1KG_samples.output.tgt,
    output:
        vcf = "results/processed_data/1KG/snps/1KG.chr{i}.biallelic.snps.vcf.gz",
        index = "results/processed_data/1KG/snps/1KG.chr{i}.biallelic.snps.vcf.gz.tbi",
    shell:
        """
        bcftools view {input.vcf} -S <(cat {input.ref} {input.tgt} | awk '{{print $2}}') | bcftools view -m 2 -M 2 -v snps | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule add_1src_1KG:
    input:
        vcf = "results/processed_data/1KG/snps/1KG.chr{i}.biallelic.snps.vcf.gz", 
        nea_vcf = "resources/data/NeaAltai/{approach}/chr{i}_mq25_mapab100.vcf.gz",
    output:
        vcf = "results/processed_data/1KG/1src/1KG.nea.{approach}.chr{i}.vcf.gz",
        index = "results/processed_data/1KG/1src/1KG.nea.{approach}.chr{i}.vcf.gz.tbi",
    resources:
        time = 2160,
    shell:
        """
        bcftools merge {input.vcf} {input.nea_vcf} | bcftools view -v snps -m 2 -M 2 -g ^miss | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule add_2src_1KG:
    input:
        vcf = rules.add_1src_1KG.output.vcf,
        den_vcf = "resources/data/Denisova/{approach}/chr{i}_mq25_mapab100.vcf.gz",
    output:
        vcf = "results/processed_data/1KG/2src/1KG.nea_den.{approach}.chr{i}.vcf.gz",
    resources:
        time = 2160,
    shell:
        """
        bcftools merge {input.vcf} {input.den_vcf} | bcftools view -v snps -m 2 -M 2 -g ^miss | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule annotate_2src_1KG:
    input:
        vcf = rules.add_2src_1KG.output.vcf,
        avsnp150 = rules.download_annovar_db.output.hg19_avsnp150,
        dbnsfp42c = rules.download_annovar_db.output.hg19_dbnsfp42c,
    output:
        avinput = "results/annotated_data/1KG/1KG.nea_den.{approach}.chr{i}.avinput",
        txt = "results/annotated_data/1KG/1KG.nea_den.{approach}.chr{i}.hg19_multianno.txt",
    resources:
        cpus=8, mem_gb=32,
    params:
        output_prefix = "results/annotated_data/1KG/1KG.nea_den.{approach}.chr{i}",
    shell:
        """
        bcftools query -f "%CHROM\t%POS\t%POS\t%REF\t%ALT\n" {input.vcf} > {output.avinput}
        resources/tools/annovar/table_annovar.pl {output.avinput} resources/tools/annovar/humandb/ -buildver hg19 -out {params.output_prefix} -remove -protocol refGene,avsnp150,dbnsfp42c -operation g,f,f -nastring . --thread {resources.cpus}
        """


rule polarize_2src_1KG:
    input:
        vcf = rules.add_2src_1KG.output.vcf,
        anc_alleles = rules.fa2bed.output.bed,
    output:
        vcf = "results/polarized_data/1KG/2src/1KG.nea_den.{approach}.chr{i}.vcf.gz",
        index = "results/polarized_data/1KG/2src/1KG.nea_den.{approach}.chr{i}.vcf.gz.tbi",
    resources:
        mem_gb=32,
    script:
        "../scripts/polarize_snps.py"
