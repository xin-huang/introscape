#----------------------------------------------------------------------------------------------------------------------

# vf: extract outliers & filter vcf to candidate regions

rule vf_outliers:
	input:
		expand("{species}/vf_out/{samp}_{chrom}_{nblock}_merged", samp=s, chrom=chrom, block=b, nblock=nb, species=config["species"])
    output:
        "{species}/vf_out/{samp}_vfoutliers.bed"
    resources: nodes=1, ntasks=1, time=45, mem_gb=80, cpus=4
    log:
        "{species}/logs/{samp}_vfoutliers.log"
    params:
        inpref = "{species}/vf_out/{samp}",
        nblock = "{nblock}",
        cutoff = config["cutoff"],
        outpref = inpref,
        R_LIBS = config["R_LIBS"],
        R_script = config["R_script"],
        scriptdir = config["scriptdir"],
        chrom = f"{{chrom}}"
    shell:
        """
        export R_LIBS={params.R_LIBS}
        {params.R_script}/Rscript --vanilla {params.scriptdir}/vf.outliers.R '{params.inpref}' {params.nblock} {params.cutoff} '{params.outpref}'
        """

#-----------------------------------------------------------------------------------------------------------------------

rule filt_vcf_to_region:
	input:
		vcf = "{species}/vf_proc_vcf/{samp}_{chrom}.filt.vcf.gz",
		bed = "{species}/vf_out/{samp}_vfoutliers.bed"
	log:
		"{species}/logs/vcftovfoutl.{samp}_{chrom}.log"
	output:
		"{downst_dir}/{samp}_{chrom}.vcf.gz"
	resources: nodes=1, ntasks=1, time_min=360, mem_gb=30, cpus=2
	shell:
		"""
		bcftools view -R  {input.bed} {input.vcf} -Oz -o {output}
		"""
#-----------------------------------------------------------------------------------------------------------------------

