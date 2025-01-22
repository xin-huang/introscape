rule all:
    input:
        "resources/tools/volcanofinder_v1.0/VolcanoFinder",


rule download_volcanofinder:
    input:
    output:
        file = "resources/tools/volcanofinder_v1.0.tar.gz",
    log:
        "logs/volcanofinder/download.log"
    shell:
        """
        wget -c https://degiorgiogroup.fau.edu/volcanofinder_v1.0.tar.gz -O {output.file} > {log} 2>&1
        """


rule decompress_volcanofinder:
    input:
        file = rules.download_volcanofinder.output.file,
    output:
        dir = directory("resources/tools/volcanofinder_v1.0/"),
    log:
        "logs/volcanofinder/decompress.log",
    shell:
        """
        tar -xvf {input.file} > {log} 2>&1
        """


rule compile_volcanofinder:
    input:
        dir = rules.decompress_volcanofinder.output.dir,
    output:
        file = "resources/tools/volcanofinder_v1.0/VolcanoFinder",
    log:
        "logs/volcanofinder/compile.log",
    shell:
        """
        cd {input.dir} && \
        sed -i -e 's/ \*data;/;/' \
               -e '22iextern struct datatype *data;' \
               -e 's/ \*data_rec;/;/' \
               -e '33iextern struct datatype_rec *data_rec;' \
               -e 's/ \*data_bvalue;/;/' \
               -e '42iextern struct datatype_bvalue *data_bvalue;' VolcanoFinder.h && \
        make > {log} 2>&1
        """
