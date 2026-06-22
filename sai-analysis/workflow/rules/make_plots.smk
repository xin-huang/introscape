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


ruleorder: plot_pr_curve > plot_roc_curve
ruleorder: plot_pr_curve_mispecify_anc_alleles > plot_pr_curve


rule plot_roc_curve:
    input:
        neu = expand(
            "results/simulated_data/run{run}/{scenario}_introgression/simulation.run{run}.{scenario}_introgression.tsv",
            run=list(range(0, num_runs)),
            allow_missing=True,
        ),
        sel = expand(
            "results/simulated_data/run{run}/adaptive_introgression/simulation.run{run}.s{s}.tsv",
            run=list(range(0, num_runs)),
            allow_missing=True,
        ),
    output:
        plot = "results/plots/roc/{scenario}_introgression_vs_adaptive_introgression_s{s}.svg",
    params:
        seed = 4836,
        replicates = 1000,
    script:
        "../scripts/plot_roc_curve.py"


rule plot_pr_curve:
    input:
        wo = expand(
            "results/simulated_data/run{run}/wo_introgression/simulation.run{run}.wo_introgression.tsv",
            run=list(range(0, num_runs)),
            allow_missing=True,
        ),
        neu = expand(
           "results/simulated_data/run{run}/neutral_introgression/simulation.run{run}.neutral_introgression.tsv",
           run=list(range(0, num_runs)),
           allow_missing=True,
        ),
        sel = expand(
            "results/simulated_data/run{run}/adaptive_introgression/simulation.run{run}.s{s}.tsv",
            run=list(range(0, num_runs)),
            allow_missing=True,
        ),
    output:
        plot = "results/plots/pr/non_adaptive_introgression_vs_adaptive_introgression_s{s}.svg",
    params:
        seed = 4836,
        replicates = 1000,
    script:
        "../scripts/plot_pr_curve.py"


rule plot_pr_curve_mispecify_anc_alleles:
    input:
        wo = expand(
            "results/simulated_data/run{run}/wo_introgression/simulation.run{run}.wo_introgression.mis.{prop}.tsv",
            run=list(range(0, num_runs)),
            allow_missing=True,
        ),
        neu = expand(
           "results/simulated_data/run{run}/neutral_introgression/simulation.run{run}.neutral_introgression.mis.{prop}.tsv",
           run=list(range(0, num_runs)),
           allow_missing=True,
        ),
        sel = expand(
            "results/simulated_data/run{run}/adaptive_introgression/simulation.run{run}.s{s}.mis.{prop}.tsv",
            run=list(range(0, num_runs)),
            allow_missing=True,
        ),
    output:
        plot = "results/plots/pr/non_adaptive_introgression_vs_adaptive_introgression_s{s}.mis.{prop}.svg",
    params:
        seed = 4836,
        replicates = 1000,
    script:
        "../scripts/plot_pr_curve.py"
