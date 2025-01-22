#-----------------------------------------------------------------------------------------------------------------------
options(stringsAsFactors=F)
options(scipen=999)
library(GenomicRanges)
library(rtracklayer)
#-----------------------------------------------------------------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
in_pref=args[1] # eg: f"vf_out/{{samp}}
nblock=args[2] # {{nblock}}
cutoff=as.numeric(args[3]) # eg 0.95
outpref=args[4]

#-----------------------------------------------------------------------------------------------------------------------
# read in volcanofinder output & extract outliers
vfout<-list()
for (i in 1:22) {
test<-read.table(paste(in_pref,"_",i,"_",nblock,"_merged",sep=""),sep="\t",header=T)
test$chr<-paste('chr',i,sep='')
vfout[[i]]<-test
}
vfoutall<-do.call(rbind,vfout)
# extract top subsets of volcanofinder scores, eg  > 95% outliers
vf_95<-vfoutall[which(vfoutall[,2]>=quantile(vfoutall$LR, probs = cutoff)),]

#-----------------------------------------------------------------------------------------------------------------------
# convert to granges
vf_95_ranges<-GRanges(seqnames=vf_95[,5],ranges=IRanges(start=as.numeric(vf_95[,1]),end=as.numeric(vf_95[,1]+1),names=vf_95[,5]),strand=rep("*",length(vf_95[,1])),LR=as.numeric(vf_95[,2]))
save(vf_95_ranges, file=paste(out_pref,"_vfoutliers",sep=""))
export.bed(vf_95_ranges, con=paste(out_pref,"_vfoutliers.bed",sep=""))
#-----------------------------------------------------------------------------------------------------------------------

# if introg regions present should intersect vf regions with those
