# In this note, we use DESeq2 to normalize the raw counts
# These codes can be run both on the HPC and on local machine. 
# Specification of home.dir decides which machine this is run.
####################################
# 0. load in the raw count matrix and generate the 
home.dir<-"/home/yanxiting/driver_Grace"
source(file.path(home.dir,"Rprogram/my_functions.R"))

rawcount.filepath<-file.path(home.dir,"scratch/GRADS/SARC_results/Results_summary_PBMC_hg38/baseline/data/rawcounts_baseline_276_clean.txt")
rawcount.matrix<-read.table(rawcount.filepath,sep="\t",header=T,check.names=F)

annot.matrix<-rawcount.matrix[,1:6]
data.matrix<-rawcount.matrix[,7:ncol(rawcount.matrix)]





