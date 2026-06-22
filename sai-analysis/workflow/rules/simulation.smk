# Copyright 2025 Xin Huang and Simon Chen
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


ruleorder: calc_stats_adaptive_introgression_mispecify_anc_alleles > calc_stats_adaptive_introgression

import numpy as np

np.random.seed(4836)
seed_list = { 
    "wo": np.random.randint(1, 2**31, num_runs),
    "neutral": np.random.randint(1, 2**31, num_runs),
    "adaptive": np.random.randint(1, 2**31, num_runs),
}

rule simulate_neutral_introgression:
    input:
        slim_script = "workflow/scripts/racimo2017_{scenario}_introgression.slim",
    output:
        final_vcf = "results/simulated_data/run{run}/{scenario}_introgression/simulation.run{run}.vcf.gz",
        anc_alleles = "results/simulated_data/run{run}/{scenario}_introgression/simulation.run{run}.anc.alleles.bed",
    params:
        src_vcf = "results/simulated_data/run{run}/{scenario}_introgression/src.run{run}.vcf",
        tgt_vcf = "results/simulated_data/run{run}/{scenario}_introgression/tgt.run{run}.vcf",
        ref_vcf = "results/simulated_data/run{run}/{scenario}_introgression/ref.run{run}.vcf",
        seed = lambda wildcards: seed_list[wildcards.scenario][int(wildcards.run)],
    shell:
        """
        slim -s {params.seed} {input.slim_script} | awk -F "#OUT:" 'BEGIN{{RS=""}}{{print $2 > "{params.src_vcf}"; print $3 > "{params.tgt_vcf}"; print $4 > "{params.ref_vcf}"}}'
        
        sed -i '/^$/d;/SV/d' {params.src_vcf}
        sed -i '/^$/d;/SV/d' {params.tgt_vcf}
        sed -i '/^$/d;/SV/d' {params.ref_vcf}
        
        bcftools reheader -s <(echo "src_i1") {params.src_vcf} | bgzip > {params.src_vcf}.gz
        bcftools reheader -s <(seq -f "tgt_i%g" 1 503) {params.tgt_vcf} | bgzip > {params.tgt_vcf}.gz
        bcftools reheader -s <(seq -f "ref_i%g" 1 1008) {params.ref_vcf} | bgzip > {params.ref_vcf}.gz
        
        tabix -p vcf {params.src_vcf}.gz
        tabix -p vcf {params.tgt_vcf}.gz
        tabix -p vcf {params.ref_vcf}.gz
        
        bcftools merge {params.src_vcf}.gz {params.tgt_vcf}.gz {params.ref_vcf}.gz --missing-to-ref | bcftools view -e 'INFO/MULTIALLELIC=1' | bgzip -c > {output.final_vcf}
        tabix -p vcf {output.final_vcf}
        bcftools query -f "%CHROM\\t%POS\\t%POS\\t%REF\\n" {output.final_vcf} | awk '{{print $1"\\t"$2-1"\\t"$3"\\t"$4}}' > {output.anc_alleles}
        rm {params.src_vcf}* {params.tgt_vcf}* {params.ref_vcf}*
        """


rule simulate_adaptive_introgression:
    input:
        slim_script = "workflow/scripts/racimo2017_adaptive_introgression.slim",
    output:
        final_vcf = "results/simulated_data/run{run}/adaptive_introgression/simulation.run{run}.s{s}.vcf.gz",
        anc_alleles = "results/simulated_data/run{run}/adaptive_introgression/simulation.run{run}.s{s}.anc.alleles.bed",
    params:
        src_vcf = "results/simulated_data/run{run}/adaptive_introgression/src.run{run}.s{s}.vcf",
        tgt_vcf = "results/simulated_data/run{run}/adaptive_introgression/tgt.run{run}.s{s}.vcf",
        ref_vcf = "results/simulated_data/run{run}/adaptive_introgression/ref.run{run}.s{s}.vcf",
        seed = lambda wildcards: seed_list["adaptive"][int(wildcards.run)],
    shell:
        """
        slim -s {params.seed} -d s={wildcards.s} {input.slim_script} | awk -F "#OUT:" 'BEGIN{{RS=""}}{{print $2 > "{params.src_vcf}"; print $3 > "{params.tgt_vcf}"; print $4 > "{params.ref_vcf}"}}'

        sed -i '/^$/d;/SV/d' {params.src_vcf}
        sed -i '/^$/d;/SV/d' {params.tgt_vcf}
        sed -i '/^$/d;/SV/d' {params.ref_vcf}

        bcftools reheader -s <(echo "src_i1") {params.src_vcf} | bgzip > {params.src_vcf}.gz
        bcftools reheader -s <(seq -f "tgt_i%g" 1 503) {params.tgt_vcf} | bgzip > {params.tgt_vcf}.gz
        bcftools reheader -s <(seq -f "ref_i%g" 1 1008) {params.ref_vcf} | bgzip > {params.ref_vcf}.gz

        tabix -p vcf {params.src_vcf}.gz
        tabix -p vcf {params.tgt_vcf}.gz
        tabix -p vcf {params.ref_vcf}.gz

        bcftools merge {params.src_vcf}.gz {params.tgt_vcf}.gz {params.ref_vcf}.gz --missing-to-ref | bcftools view -e 'INFO/MULTIALLELIC=1' | bgzip -c > {output.final_vcf}
        tabix -p vcf {output.final_vcf}
        bcftools query -f "%CHROM\\t%POS\\t%POS\\t%REF\\n" {output.final_vcf} | awk '{{print $1"\\t"$2-1"\\t"$3"\\t"$4}}' > {output.anc_alleles}
        rm {params.src_vcf}* {params.tgt_vcf}* {params.ref_vcf}*
        """


rule calc_stats_neutral_introgression:
    input:
        vcf = rules.simulate_neutral_introgression.output.final_vcf,
        anc_alleles = rules.simulate_neutral_introgression.output.anc_alleles,
    output:
        tsv = "results/simulated_data/run{run}/{scenario}_introgression/simulation.run{run}.{scenario}_introgression.tsv",
    params:
        config = "config/simulation/simulated_data.yaml",
        chr_name = 1,
        win_len = 40000,
        win_step = 40000,
    shell:
        """
        sai score --vcf {input.vcf} --chr-name {params.chr_name} --config {params.config} --output {output.tsv} --win-len {params.win_len} --win-step {params.win_step} --anc-alleles {input.anc_alleles}
        """


rule calc_stats_adaptive_introgression:
    input:
        vcf = rules.simulate_adaptive_introgression.output.final_vcf,
        anc_alleles = rules.simulate_adaptive_introgression.output.anc_alleles,
    output:
        tsv = "results/simulated_data/run{run}/adaptive_introgression/simulation.run{run}.s{s}.tsv",
    params:
        config = "config/simulation/simulated_data.yaml",
        chr_name = 1,
        win_len = 40000,
        win_step = 40000,
    shell:
        """
        sai score --vcf {input.vcf} --chr-name {params.chr_name} --config {params.config} --output {output.tsv} --win-len {params.win_len} --win-step {params.win_step} --anc-alleles {input.anc_alleles}
        """


rule mispecify_anc_alleles_neutral_introgression:
    input:
        anc_alleles = rules.simulate_neutral_introgression.output.anc_alleles,
    output:
        anc_alleles = "results/simulated_data/run{run}/{scenario}_introgression/simulation.run{run}.anc.alleles.mis.{prop}.bed",
    params:
        seed = lambda wildcards: seed_list[wildcards.scenario][int(wildcards.run)],
    script:
       "../scripts/mispecify_anc_alleles.py"


rule mispecify_anc_alleles_adaptive_introgression:
    input:
        anc_alleles = rules.simulate_adaptive_introgression.output.anc_alleles,
    output:
        anc_alleles = "results/simulated_data/run{run}/adaptive_introgression/simulation.run{run}.s{s}.anc.alleles.mis.{prop}.bed",
    params:
        seed = lambda wildcards: seed_list["adaptive"][int(wildcards.run)],
    script:
       "../scripts/mispecify_anc_alleles.py"


rule calc_stats_neutral_introgression_mispecify_anc_alleles:
    input:
        vcf = rules.simulate_neutral_introgression.output.final_vcf,
        anc_alleles = rules.mispecify_anc_alleles_neutral_introgression.output.anc_alleles,
    output:
        tsv = "results/simulated_data/run{run}/{scenario}_introgression/simulation.run{run}.{scenario}_introgression.mis.{prop}.tsv",
    params:
        config = "config/simulation/simulated_data.yaml",
        chr_name = 1,
        win_len = 40000,
        win_step = 40000,
    shell:
        """
        sai score --vcf {input.vcf} --chr-name {params.chr_name} --config {params.config} --output {output.tsv} --win-len {params.win_len} --win-step {params.win_step} --anc-alleles {input.anc_alleles}
        """


rule calc_stats_adaptive_introgression_mispecify_anc_alleles:
    input:
        vcf = rules.simulate_adaptive_introgression.output.final_vcf,
        anc_alleles = rules.mispecify_anc_alleles_adaptive_introgression.output.anc_alleles,
    output:
        tsv = "results/simulated_data/run{run}/adaptive_introgression/simulation.run{run}.s{s}.mis.{prop}.tsv",
    params:
        config = "config/simulation/simulated_data.yaml",
        chr_name = 1,
        win_len = 40000,
        win_step = 40000,
    shell:
        """
        sai score --vcf {input.vcf} --chr-name {params.chr_name} --config {params.config} --output {output.tsv} --win-len {params.win_len} --win-step {params.win_step} --anc-alleles {input.anc_alleles}
        """
