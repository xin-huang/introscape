# annotate variants

#-----------------------------------------------------------------------------------------------------------------------

# download indexed vep cache

# eg snpeff_reference=GRCh38.99
# https://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/homo_sapiens_vep_110_GRCh38.tar.gz

# vep_species = homo_sapiens
# vep_build = GRCh38
# {vep_release} = 110


rule get_indexed_vep_cache:
    output:
        directory("resources/vep/indexed_cache"),
    params:
        species="{vep_species}",
        build="{vep_build}",
        release="{vep_release}",
        indexed=True,
    log:
        "logs/vep/indexed_cache.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    wrapper:
        "v5.5.0/bio/vep/cache"

# download vep plugin

rule download_vep_plugins:
    output:
        directory("resources/vep/plugins")
    params:
        release=100
    log:
        "logs/vep/download_vep_plugins.log",
    wrapper:
        "v5.5.0/bio/vep/plugins"


# annotate variants with vep

rule annotate_variants:
    input:
        calls="{downst_dir}/{samp}_{chrom}.vcf.gz",  # .vcf, .vcf.gz or .bcf
        cache="resources/vep/indexed_cache",  # can be omitted if fasta and gff are specified
        plugins="resources/vep/plugins"
    output:
        calls="{annt}/{samp}_{chrom}_ann.vcf.gz"  # .vcf, .vcf.gz or .bcf
    params:
        # Pass a list of plugins to use, see https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
            # include polyphen & sift info 
        plugins=["PolyPhen_SIFT"],
        extra="--everything"  # optional: extra arguments
    log:
        "logs/vep/annotate_{samp}_{chrom}.log",
    threads: 4
    wrapper:
        "v5.5.0/bio/vep/annotate"

