#-----------------------------------------------------------------------------------------------------------------------
# generate allele freq file
	# count allele freq (1st awk command)
	# add header line (echo command)
	# filter to remove  0 sites (ref homozygotes) (2nd awk command)

rule make_allele_freq:
	input:
		fvcf = "{species}/vf_proc_vcf/{samp}_{chrom}.filt.vcf.gz"
	log:
		"{species}/logs/make.allelefreq.file.{samp}_{chrom}.log"
	output:
		"{species}/vf_proc_vcf/allelefreq/{samp}_{chrom}_SF2.input",
		"{species}/vf_proc_vcf/allelefreq/{samp}_{chrom}_tmp.txt",
		"{species}/vf_proc_vcf/allelefreq/{samp}_{chrom}_SF2_polym.input"
	resources: nodes=1, ntasks=1, time_min=360, mem_gb=60, cpus=2
	shell:
		"""
		vcftools --counts2 --derived --gzvcf {input.fvcf} --stdout | awk 'NR<=1 {{next}} {{print $2"\t"$6"\t"$4"\t0"}}' > {output[0]}
		echo -e "position\tx\tn\tfolded" | cat - {output[0]} > {output[1]}  && mv {output[1]} {output[0]}
		awk '{{ if ($2 != 0) {{ print }} }}' {output[0]} > {output[2]}
		"""

#-----------------------------------------------------------------------------------------------------------------------

# convert vcf -> genotype matrix ( genotypes coded as 0, 1, 2, -1 (missing data))

# 1) bcftools query: extract only genotypes
# 2) sed: Replace 0/0 by 0 # Replace 0/1 by 1 # Replace 1/1 by 2 # Replace ./. by -1 # Replace . by -1

rule vcf_to_gtmatrix:
	input:
		fvcf = "{species}/vf_proc_vcf/{samp}_{chrom}.filt.vcf.gz"
	log:
		"{species}/logs/vcf_to_gtmatrix.{samp}_{chrom}.log"
	output:
		"{species}/vf_proc_vcf/gtmatrix/{samp}_{chrom}.GT",

	resources: nodes=1, ntasks=1, time_min=360, mem_gb=30, cpus=2
	shell:
		"""
		bcftools query -f '[%GT\t]\n' {input.fvcf} > {output};
		sed -i "s/0|0/0/g;s/0|1/1/g;s/1|0/1/g;s/1|1/2/g;s/\.|./-1/g;s/\./-1/g" {output}
		"""

#-----------------------------------------------------------------------------------------------------------------------

# calc sfs per chrom in R - 1st R script

rule sfs_per_chrom:
	input:
		fvcf = "{species}/vf_proc_vcf/{samp}_{chrom}.filt.vcf.gz",
		gtinp = "{species}/vf_proc_vcf/gtmatrix/{samp}_{chrom}.GT",
	log:
		"{species}/logs/sfs_per_chrom.{samp}_{chrom}.log"
	output:
		"{species}/vf_proc_vcf/sfs/{samp}_{chrom}_sfscounts",
	resources: nodes=1, ntasks=1, time_min=480, mem_gb=60, cpus=2
	params:
		scriptdir = {config['scriptdir']}
	shell:
		"""
		Rscript --vanilla {params.scriptdir}/vf.sfs.1.R '{input.fvcf}' '{input.gtinp}' '{output}'
		"""

#-----------------------------------------------------------------------------------------------------------------------

# normalise sfs to format req for volcanofinder

rule sfs_for_vf:
	input:
		expand("{species}/vf_proc_vcf/sfs/{samp}_{chrom}_sfscounts", samp=config["samp"], chrom=chrom, species=config["species"])
	log:
		"{species}/logs/sfs_for_vf.{samp}.log"
	output:
		"{species}/vf_proc_vcf/sfs/{samp}_autosomes_sfscounts",
		"{species}/vf_proc_vcf/sfs/{samp}_autosomes_sfs.input.txt"
	resources: nodes=1, ntasks=1, time_min=360, mem_gb=60, cpus=2
	params:
		scriptdir = config["scriptdir"],
		sp = "{samp}"
	shell:
		"""
		Rscript --vanilla {params.scriptdir}/vf.sfs.2.R '{params.sp}'
		"""
#-----------------------------------------------------------------------------------------------------------------------
