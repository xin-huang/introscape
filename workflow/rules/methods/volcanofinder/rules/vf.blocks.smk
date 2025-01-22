#-----------------------------------------------------------------------------------------------------------------------
# run vf in arbitary number of blocks

# expanlation of params:
# test size every kilobase, then 𝑔 = 1000
# FreqFile = allele frequency input
# SpectFile = empirical unnormalized site frequency

# BLOCK = current block to compute test sites 
# NBLOCK =  number of contiguous blocks in which to break a chromosome.

rule run_vf_in_blocks:
	input:
		FreqFile = "{species}/vf_proc_vcf/allelefreq/{samp}_{chrom}_SF2_polym.input",
		SpectFile = "{species}/vf_proc_vcf/sfs/{samp}_autosomes_sfs.input.txt"
	log:
		"{species}/logs/run_vf_in_blocks.{samp}_{chrom}_{block}_{nblock}.log"
	output:
		"{species}/vf_out/{samp}_{chrom}_{block}_{nblock}"
	resources: nodes=1, ntasks=1, time_min=360, mem_gb=30, cpus=2
	params:
		scriptdir = {config['scriptdir']},
		outfile = "vf_out/{samp}_{chrom}",
		block = lambda wildcards: wildcards.block,
		nblock = lambda wildcards: wildcards.nblock
	shell:
		"""
		{params.scriptdir}/VolcanoFinder -big 1000 {input.FreqFile} {input.SpectFile} -1 1 1 {params.outfile} {params.block} {params.nblock}
		"""

#-----------------------------------------------------------------------------------------------------------------------
rule merge_vf_in_blocks:
	input:
		expand("{species}/vf_out/{samp}_{chrom}_{block}_{nblock}", samp=s, chrom=chrom, block=b, nblock=nb, species=config["species"])
	log:
		"{species}/logs/merge_vf_in_blocks.{samp}_{chrom}_{nblock}.log"
	output:
		"{species}/vf_out/{samp}_{chrom}_{nblock}_merged"
	resources: nodes=1, ntasks=1, time_min=360, mem_gb=30, cpus=2
	params:
		scriptdir = {config['scriptdir']},
		outfile = "{species}/vf_out/{samp}_{chrom}_{nblock}",
		nblock = lambda wildcards: wildcards.nblock
	shell:
		"""
		{params.scriptdir}/VolcanoFinder -m {params.outfile} {params.nblock}
		"""
#-----------------------------------------------------------------------------------------------------------------------
