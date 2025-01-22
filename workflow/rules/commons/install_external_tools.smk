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


rule install_ms_for_sstar:
    input:
    output:
        "resources/msdir/ms",
    log:
        "logs/install_ms_for_sstar/install_ms_for_sstar.log",
    shell:
        """
        cd resources 2>> {log}
        tar -xvf ms.tar.gz 2>> {log}
        cd msdir 2>> {log}
        gcc -o ms ms.c streec.c rand1.c -lm 2>> {log}
        """


rule install_SPrime:
    input:
    output:
        "resources/SPrime/sprime.jar",
        directory("resources/SPrime/sprimepipeline"),
    log:
        "logs/install_SPrime/install_SPrime.log",
    shell:
        """
        cd resources 2>> {log}
        mkdir SPrime && cd SPrime 2>> {log}
        wget https://faculty.washington.edu/browning/sprime.jar 2>> {log}
        git clone https://github.com/YingZhou001/sprimepipeline 2>> {log}
        chmod a+x sprimepipeline/pub.pipeline.pbs/tools/map_arch_genome/map_arch 2>> {log}
        sed 's/out<-c()/out<-data.frame()/' sprimepipeline/pub.pipeline.pbs/tools/score_summary.r > tmp 2>> {log}
        mv tmp sprimepipeline/pub.pipeline.pbs/tools/score_summary.r 2>> {log}
        """
