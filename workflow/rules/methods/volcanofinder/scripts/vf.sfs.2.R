require(data.table)
options(stringsAsFactors=F)
args = commandArgs(trailingOnly=TRUE)
species=args[1]

savelist<-list()
for (a in 1:22){
load(file=paste(species,"_",a,"_sfscounts",sep=""), verbose=TRUE)
savelist[[a]]<-obs.sfs 
}
save(savelist,file=paste(species,"_autosomes_sfscounts",sep=""))

# unlist & sum
all<-data.frame(matrix(unlist(savelist), ncol=length(savelist[[1]]), byrow=TRUE),stringsAsFactors=FALSE)

aut.obs.sfs <-colSums(all)

# then only normalise once have calculated site categories for all chromosomes
# normalise 
test<-aut.obs.sfs/sum(aut.obs.sfs)

# upper limit of seq needs to be (nind)*2

ndiploid<-nind*2 

names(test)<-c(0, seq(1:ndiploid))

# drop first category
test1<-test[-1]

x<-as.data.table(test1)

x$V1<-as.numeric(names(test1))


# columns switched around -> need to be  1"\t"sfs_val 
# amend this before outputting as txt file

# ie correctly formatted sfs is the following - 
#(base) harvi@Harvinders-MacBook-Pro all_candidates % head  ~/Downloads/new_4_autosomes_sfs.input.txt
#1   0.002095229649634
#2   0.0009180433101545

y<-x[,c(2,1)]


write.table(y, file=paste(species,"_autosomes_sfs.input.txt",sep=""),sep="\t",row.names=FALSE,col.names=FALSE)

#-----------------------------------------------------------------------------------------------------------------------
# Derek Setter: The problem is simply that you had the columns switched around in your SFS file.

#[hpawar@cnb2 ~]$ head /scratch/devel/hpawar/admix/volcanofinder/input/2feb22/tmp/4_autosomes_sfs.input.txt
#0.00209522964963409 1
#0.000918043310154512    2
#0.000794831420940574    3

# ie correctly formatted sfs is the following - 
#(base) harvi@Harvinders-MacBook-Pro all_candidates % head  ~/Downloads/new_4_autosomes_sfs.input.txt
#1   0.002095229649634
#2   0.0009180433101545
#3   0.0007948314209405

#-----------------------------------------------------------------------------------------------------------------------
