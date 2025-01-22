# 1) filter to ingroup ids, biallelic snps & no missing 

rule filt_ingroup_biallelic_snp_nomissing:
	input:
		vcf = "{config['inpath']}chr{chrom}.vcf.gz",
		sampl = "{config['samp_dir']}/{samp}.samples"
	log:
		"{species}/logs/vf_proc_vcf_filt1.{samp}_{chrom}.log"
	output:
		"{species}/vf_proc_vcf/{samp}_{chrom}.vcf.gz"
	resources: nodes=1, ntasks=1, time_min=480, mem_gb=60, cpus=2
	shell:
		"""
		bcftools view -S {input.sampl} {input.vcf} | bcftools view -m2 -M2 -v snps - | bcftools filter --exclude 'GT~"\."' -Oz -o {output}
		"""

#-----------------------------------------------------------------------------------------------------------------------
# 2) tabix vcf

rule tabix:
	input:
		vcf = "{species}/vf_proc_vcf/{samp}_{chrom}.vcf.gz"
	log:
		"{species}/logs/vf_proc_vcf_tabix.{samp}_{chrom}.log"
	output:
		"{species}/vf_proc_vcf/{samp}_{chrom}.vcf.gz.tbi"
	resources: nodes=1, ntasks=1, time_min=60, mem_gb=30, cpus=1
	shell:
		"""
		tabix -p vcf {input.vcf}
		"""

#-----------------------------------------------------------------------------------------------------------------------

rule annotate_ancestral_alleles:
    input:
        vcf = "{species}/vf_proc_vcf/{samp}_{chrom}.vcf.gz",
        ancbed = "{config['ancbed']}_{chrom}.bed"
    output:
        anc_vcf = "{species}/vf_proc_vcf/{samp}_{chrom}.filt.vcf.gz"
    resources: nodes=1, ntasks=1, time=360, mem_gb=80, cpus=4
    log:
        "{species}/logs/reannot_anc_allels_{samp}_{chrom}.log"
    params:
        hd_dir = config["hd_dir"]
    shell:
        """
		bcftools annotate -a {input.ancbed} -h {params.hd_dir}/annots.hdr -c CHROM,FROM,TO,AA {input.vcf} -Oz -o {output.anc_vcf}
        """

#-----------------------------------------------------------------------------------------------------------------------
