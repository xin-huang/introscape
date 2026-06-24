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


rule candidate_genes_html:
    input:
        tsv=rules.annotate_candidate_regions.output.genes,
    output:
        html="results/sstar2/{analysis}/candidates/{analysis}.sstar2.candidate.genes.html",
    params:
        title=get_method_config_title,
    log:
        "logs/report/candidate_genes_html.{analysis}.log",
    conda:
        "../envs/introscape.yaml"
    script:
        "../scripts/tsv_to_html.py"
