rule install_smk:
    shell:
        """
        mkdir SPrime && cd SPrime
        wget https://faculty.washington.edu/browning/sprime.jar
        git clone https://github.com/YingZhou001/sprimepipeline
        chmod a+x sprimepipeline/pub.pipeline.pbs/tools/map_arch_genome/map_arch
        sed 's/out<-c()/out<-data.frame()/' sprimepipeline/pub.pipeline.pbs/tools/score_summary.r > tmp
        mv tmp sprimepipeline/pub.pipeline.pbs/tools/score_summary.r
        """
