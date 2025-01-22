require(data.table)
options(stringsAsFactors=F)
args = commandArgs(trailingOnly=TRUE)

# input vcf (filtered, only easterns, biallelic, no multiallelic ancestral states)
inputvcf=args[1]
# input (GT file)
input=args[2]
# output file
output=args[3]

# calc sfs per chrom

samples.in.vcf<-system(paste("bcftools query -l ",inputvcf," ",sep=""),intern=T) 

# number of individuals in ingroup
nind <- length(samples.in.vcf)

# number of SNPs per GT file
nrows_a <- system(paste("cat ",input," | wc -l",sep=""))

# calc number of snps per gt file

# read in the GT matrix for this chr
md_geno <- matrix(scan(file=input), ncol=nind, nrow=nrows_a, byrow=T)

obs.sfs<-table(rowSums(md_geno))
obs.sfs # should be 0 -2n

# write out the individual r objects (sfs of each chr)
save(obs.sfs,file=output)

q()
#-----------------------------------------------------------------------------------------------------------------------
