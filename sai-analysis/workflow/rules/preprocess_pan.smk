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


rule create_pan_samples:
    input:
        metadata = rules.download_pan_metadata.output.txt,
    output:
        ref = "results/processed_data/pan/samples/pan.ref.samples.txt",
        tgt = "results/processed_data/pan/samples/pan.tgt.samples.txt",
        src = "results/processed_data/pan/samples/pan.src.samples.txt",
    shell:
        """
        grep -v captive {input.metadata} | awk '$2~/^PTV$/{{print $2"\\t"$4}}' > {output.ref}
        grep -v captive {input.metadata} | awk '$2~/^PTT$/{{print $2"\\t"$4}}' > {output.tgt}
        grep -v captive {input.metadata} | awk '$2~/^PPA$/{{print $2"\\t"$4}}' > {output.src}
        """


rule extract_pan_biallelic_snps:
    input:
        vcf = rules.download_pan_genomes.output.vcf,
        ref = rules.create_pan_samples.output.ref,
        tgt = rules.create_pan_samples.output.tgt,
        src = rules.create_pan_samples.output.src,
    output:
        vcf = "results/processed_data/pan/1src/pan.chr{i}.biallelic.snps.vcf.gz",
        index = "results/processed_data/pan/1src/pan.chr{i}.biallelic.snps.vcf.gz.tbi",
    shell:
        """
        bcftools view {input.vcf} -S <(cat {input.ref} {input.tgt} {input.src} | awk '{{print $2}}') | bcftools view -v snps -m 2 -M 2 -g ^miss | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule annotate_pan:
    input:
        vcf = rules.extract_pan_biallelic_snps.output.vcf,
        refgene = rules.download_annovar_db.output.hg38_refgene,
        avsnp150 = rules.download_annovar_db.output.hg38_avsnp150,
        dbnsfp42c = rules.download_annovar_db.output.hg38_dbnsfp42c,
    output:
        avinput = "results/annotated_data/pan/pan.chr{i}.biallelic.snps.avinput",
        txt = "results/annotated_data/pan/pan.chr{i}.biallelic.snps.hg38_multianno.txt",
    resources:
        cpus=8, mem_gb=32,
    params:
        output_prefix = "results/annotated_data/pan/pan.chr{i}.biallelic.snps",
    shell:
        """
        bcftools query -f "%CHROM\\t%POS\\t%POS\\t%REF\\t%ALT\\n" {input.vcf} > {output.avinput}
        resources/tools/annovar/table_annovar.pl {output.avinput} resources/tools/annovar/humandb/ -buildver hg38 -out {params.output_prefix} -remove -protocol refGene,avsnp150,dbnsfp42c -operation g,f,f -nastring . --thread {resources.cpus}
        """


rule extract_pan_anc_info:
    input:
        hg38 = rules.download_refgenomes.output.hg38,
        rheMac10 = rules.download_refgenomes.output.rheMac10,
        chain = rules.download_refgenomes.output.chain,
    output:
        anc_alleles = "results/polarized_data/anc_info/hg38.chr{i}.anc.alleles.bed",
    resources:
        mem_gb = 128,
    shell:
        """
        python workflow/scripts/get_ancestral_info.py --src-fasta {input.hg38} --tgt-fasta {input.rheMac10} --liftover-chain {input.chain} --chr-name chr{wildcards.i} --output {output.anc_alleles}
        """


rule polarize:
    input:
        vcf = rules.extract_pan_biallelic_snps.output.vcf,
        anc_alleles = rules.extract_pan_anc_info.output.anc_alleles,
    output:
        vcf = "results/polarized_data/pan/pan.chr{i}.biallelic.snps.vcf.gz",
        index = "results/polarized_data/pan/pan.chr{i}.biallelic.snps.vcf.gz.tbi",
    resources:
        mem_gb = 32,
    script:
        "../scripts/polarize_snps.py"
