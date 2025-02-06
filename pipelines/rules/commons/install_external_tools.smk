# GNU General Public License v3.0
# Copyright 2024 Xin Huang
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

rule all:
    input:
        "resources/msdir/ms",
        "resources/SPrime/sprime.jar",
        "resources/flags/.skovscripts.downloaded"


rule install_ms_for_sstar:
    output:
        "resources/msdir/ms",
    log:
        "logs/install_ms_for_sstar/install_ms_for_sstar.log",
    shell:
        """
        cd resources
        tar -xvf ms.tar.gz
        cd msdir
        gcc -o ms ms.c streec.c rand1.c -lm
        """


rule install_SPrime:
    output:
        "resources/SPrime/sprime.jar",
        directory("resources/SPrime/sprimepipeline"),
    log:
        "logs/install_SPrime/install_SPrime.log",
    params:
        outdir = "resources/SPrime"
    shell:
        """
	mkdir -p {params.outdir}
	cd {params.outdir}
        wget https://faculty.washington.edu/browning/sprime.jar
        git clone https://github.com/YingZhou001/sprimepipeline
        chmod a+x sprimepipeline/pub.pipeline.pbs/tools/map_arch_genome/map_arch
        sed 's/out<-c()/out<-data.frame()/' sprimepipeline/pub.pipeline.pbs/tools/score_summary.r > tmp
        mv tmp sprimepipeline/pub.pipeline.pbs/tools/score_summary.r
        """

rule get_hmmix_scripts:
    output:
        download_flag = "resources/flags/.skovscripts.downloaded"
    resources: nodes=1, ntasks=1, time_min=60, mem_gb=30, cpus=1
    log:
        "logs/gitclone_skovscripts.log"
    params:
        outdir = "resources/skov_scripts"
    shell:
        """
        mkdir -p {params.outdir}
        cd {params.outdir}
        git clone https://github.com/LauritsSkov/Introgression-detection
	cd ../../
	touch {output.download_flag}
        """
