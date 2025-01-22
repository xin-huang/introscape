# volcanofinder_to_smk

convert steps necessary to run volcanofinder to a smk pipeline
https://github.com/h-pawar/gor_ghost_introg/tree/main/7.Volcanofinder/7.1.processforvf/scripts

Separate rules per step
 convert to smk file (0) ancestral allele states (1) processing data (2) apply vf (3) downstream steps

* run.vf.smk
 * calls rules & config files

- preprocess vcf, incl. add ancestral allele info
- generate input files for volcanofinder
- run volcanofinder
- outlier regions of vf statistic -> candidate regions
- if bed files for introg regions (eg hmmix &/or s*) are present overlap with the vf outliers -> candidate regions
- filter vcf to candidate regions
- annotate candidate regions using snpeff or vep info

#-----------------------------------------------------------------------------------------------------------------------

how to calc ancestral allele states not yet included
