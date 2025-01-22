# annotate variants

#-----------------------------------------------------------------------------------------------------------------------

# 1) snpeff

# annt = annotation type (here annt= snpeff)

# download snpeff database for given species 
# eg snpeff_reference=GRCh38.99

# list prebuilt databases: 
# java -jar snpEff.jar databases 
# java -jar $SNPEFF/snpEff.jar databases | grep Homo_sapiens

rule snpeff_download:
    output:
        directory("resources/snpeff/{snpeff_reference}")
    log:
        "logs/snpeff/download/{snpeff_reference}.log"
    params:
        reference="{snpeff_reference}"
    resources:
        mem_mb=1024
    wrapper:
        "v5.5.0/bio/snpeff/download"


# java -Xmx8g -jar path/to/snpEff/snpEff.jar  -c path/to/snpEff/snpEff.config -v GRCh38.99 in.vcf.gz > in.ann.vcf.gz

# apply snpeff 
rule snpeff:
    input:
        calls="{downst_dir}/{samp}_{chrom}.vcf.gz", # (vcf, bcf, or vcf.gz)
        db="resources/snpeff/{snpeff_reference}" # path to reference db downloaded with the snpeff download wrapper
    output:
        calls="{annt}/{samp}_{chrom}_ann.vcf.gz"   # annotated calls (vcf, bcf, or vcf.gz)
    log:
        "logs/snpeff/{samp}_{chrom}_ann.log"
    resources:
        java_opts="-XX:ParallelGCThreads=10",
        mem_mb=4096
    wrapper:
        "v5.5.0/bio/snpeff/annotate"

#-----------------------------------------------------------------------------------------------------------------------
