import os

configfile: "vf.gape.config.yaml"

chrom = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

b = config["block"]
nb = config["nblock"]
s = config["samp"]
downst_dir = config["downst_dir"]
#-----------------------------------------------------------------------------------------------------------------------

rule all:
	input:
		expand("{species}/vf_proc_vcf/{samp}_{chrom}.filt.vcf.gz", samp=s, chrom=chrom,species=config["species"]),
		expand("{species}/vf_proc_vcf/allelefreq/{samp}_{chrom}_SF2_polym.input", samp=s, chrom=chrom,species=config["species"]),
		expand("{species}/vf_proc_vcf/sfs/{samp}_autosomes_sfs.input.txt", samp=s,species=config["species"]),
		expand("{species}/vf_proc_vcf/{samp}_{chrom}_{block}_{nblock}", samp=s, chrom=chrom, block=b, nblock=nb,species=config["species"]),
		expand("{species}/vf_proc_vcf/{samp}_{chrom}_merged", samp=s, chrom=chrom,species=config["species"]),
		expand("{species}/vf_out/{samp}_vfoutliers.bed", samp=s, species=config["species"]),
		expand("{annt}/{samp}_{chrom}_ann.vcf.gz", samp=s, chrom=chrom, annt=config["annt"]))

#-----------------------------------------------------------------------------------------------------------------------

# preprocessing dependent on data type (diff naming conventions of the vcfs)
if config["type"] == "1000g":
    include: "rules/proc_1000g.smk"
elif config["type"] == "gape":
    include: "rules/proc_gape.smk"

# generate input files for volcanofinder & apply volcanofinder
include: "rules/af.sfs.smk"
include: "rules/vf.blocks.smk"

# filter volcanofinder candidates & prioritise vcf for annotation
# if introg regions from s* &/or hmmix are present - overlap these with vf candidates & then prioritise vcf for annotation

if config["introg_regions_present"] == "yes":
    include: "rules/vf.downst2.smk"
elif config["introg_regions_present"] == "no":
    include: "rules/vf.downst1.smk"

# annotate candidate variants in downstream vcf with snpeff or vep
    # ie vcf after filtering by vf candidates
if config['annotate_downst'] == "yes" and config['annt'] == "snpeff":
    include: "rules/snpeff.smk"
if config['annotate_downst'] == "yes" and config['annt'] == "vep":
    include: "rules/vep.smk"
