#-----------------------------------------------------------------------------------------------------------------------
options(stringsAsFactors=F)
options(scipen=999)
library(GenomicRanges)
library(rtracklayer)
#-----------------------------------------------------------------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)
# args here should be prefixes to files of the outliers of -
vf_outl=args[1] # volcanofinder
intro1=args[2] # sstar 
intro2=args[3] # hmmix
out_pref=args[4]
#-----------------------------------------------------------------------------------------------------------------------
bed_to_ranges<-function(in_pref) {
in_bed<-read.table(paste(in_pref,".bed",sep=""),sep="\t",header=F)
# bed = 0-based, gr = 1-based
in_ranges<-GRanges(seqnames=in_bed[,1],ranges=IRanges(start=as.numeric(in_bed[,2]+1),end=as.numeric(in_bed[,3]),names=in_bed[,1]),strand=rep("*",length(in_bed[,1])))
return(in_ranges)}

# overlap vf outliers with introg regions
vf_w_introg_overlap_fun<-function(vf_obj,introg_obj){
out_ov<-subsetByOverlaps(vf_obj,introg_obj)
return(out_ov)}
#-----------------------------------------------------------------------------------------------------------------------

vf_outliers<-bed_to_ranges(vf_outl)

#-----------------------------------------------------------------------------------------------------------------------

# if outliers of both hmmix & s* are provided - first extract the unique introg regions & query by this
	# else query by the single outliers provided (either s* or hmmix)
if (intro1 != "" & intro2 != "") {
sstar_subs<-bed_to_ranges(intro1)
skov_subs<-bed_to_ranges(intro2)
introg_overlap<-reduce(subsetByOverlaps(skov_subs,sstar_subs)) # unique overlapping regions
out_vi<-vf_w_introg_overlap_fun(vf_outliers, introg_overlap)
} else if (intro1 != "" & intro2 == "") {
sstar_subs<-bed_to_ranges(intro1)
out_vi<-vf_w_introg_overlap_fun(vf_outliers, sstar_subs)
} else if (intro1 == "" & intro2 != "") {
skov_subs<-bed_to_ranges(intro2)
out_vi<-vf_w_introg_overlap_fun(vf_outliers, skov_subs)
} else {
print("need to provide introgressed regions")
}

export.bed(out_vi, con=paste(out_pref,"_vf_introg_outliers.bed",sep=""))

#-----------------------------------------------------------------------------------------------------------------------
