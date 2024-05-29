rule install_ms_for_sstar:
    shell:
        """
        echo $CONDA_PREFIX
        mkdir ms
        cd ms
        tar -xvf ../pipelines/envs/ms.tar.gz
        cd msdir
        gcc -o ms ms.c streec.c rand1.c -lm
        """
