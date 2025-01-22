#----------------------------------------------------------------------------------------------------------------------

# vf: extract outliers, overlap with s* &/or hmmix regions & filter vcf to candidate regions

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

rule overlap_vf_intro:
    input:
        bed = "{species}/vf_out/{samp}_vfoutliers.bed"
    output:
        "{species}/vf_out/{samp}_vf_introg_outliers.bed"
    resources: nodes=1, ntasks=1, time=45, mem_gb=80, cpus=4
    log:
        "{species}/logs/{samp}_vf_introg_outliers.log"
    params:
        inpref = "{species}/vf_out/{samp}_vfoutliers",
        intro1 = config.get("sstar_bed", ""),
        intro2 = config.get("hmmix_bed", ""),
        outpref = "{species}/vf_out/{samp}",
        R_LIBS = config["R_LIBS"],
        R_script = config["R_script"],
        scriptdir = config["scriptdir"]
    shell:
        """
        export R_LIBS={params.R_LIBS}
        {params.R_script}/Rscript --vanilla {params.scriptdir}/overlap.vf.introg.R '{params.inpref}' '{params.intro1}' '{params.intro2}' '{params.outpref}'
        """


#-----------------------------------------------------------------------------------------------------------------------

rule ovfilt_vcf_to_region:
	input:
		vcf = "{species}/vf_proc_vcf/{samp}_{chrom}.filt.vcf.gz",
		bed = "{species}/vf_out/{samp}_vf_introg_outliers.bed"
	log:
		"{species}/logs/ovvcftovfoutl.{samp}_{chrom}.log"
	output:
		"{downst_dir}/{samp}_{chrom}.vcf.gz"
	resources: nodes=1, ntasks=1, time_min=360, mem_gb=30, cpus=2
	shell:
		"""
		bcftools view -R  {input.bed} {input.vcf} -Oz -o {output}
		"""
#-----------------------------------------------------------------------------------------------------------------------