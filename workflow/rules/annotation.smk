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


rule annotate_candidate_regions:
    input:
        candidates=rules.convert_sstar2_tracts_to_candidates.output.candidates,
        annotation=get_genome_annotation,
    output:
        annotated="results/sstar2/{analysis}/candidates/{analysis}.sstar2.annotated.candidate.regions.tsv",
        genes="results/sstar2/{analysis}/candidates/{analysis}.sstar2.candidate.genes.tsv",
    log:
        "logs/annotation/annotate_candidate_regions.{analysis}.log",
    conda:
        "../envs/introscape.yaml"
    script:
        "../scripts/annotate_candidate_regions.py"
