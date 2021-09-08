# In this note, we use DESeq2 to normalize the raw counts
# These codes can be run both on the HPC and on local machine. 
# Specification of home.dir decides which machine this is run.

# https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html

####################################
# 0. load in the raw count matrix and generate the 
library(DESeq2)
home.dir<-"/home/yanxiting/driver_Grace"
source(file.path(home.dir,"Rprogram/my_functions.R"))

rawcount.filepath<-file.path(home.dir,"scratch/GRADS/SARC_results/Results_summary_PBMC_hg38/baseline/data/rawcounts_baseline_276_clean.txt")
rawcount.matrix<-read.table(rawcount.filepath,sep="\t",header=T,check.names=F)

annot.matrix<-rawcount.matrix[,1:6]
data.matrix<-rawcount.matrix[,7:ncol(rawcount.matrix)]

# generate the meta data
meta.data<-readRDS(file.path(home.dir,"scratch/GRADS/SARC_results/Results_summary_PBMC_hg38/baseline/data/clinic_matrix_merged.RDS"),refhook = NULL)

# match the row names of meta data and the column names of data.matrix
meta.data<-meta.data[colnames(data.matrix),]

dds <- DESeqDataSetFromMatrix(countData = data.matrix, colData = meta.data, design = ~ 1)
#View(counts(dds))
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

normalized_counts <- counts(dds, normalized=TRUE)

cmd.out<-cbind(annot.matrix,normalized_counts)

# output the normalized counts
output.filepath<-file.path(home.dir,"scratch/GRADS/SARC_results/Results_summary_PBMC_hg38/baseline/data","DESeq2_normalized_276_clean.txt")
write.table(cmd.out,file=output.filepath,append=F,sep="\t",row.names=F,col.names=T,quote=F)

output.filepath<-file.path(home.dir,"scratch/GRADS/SARC_results/Results_summary_PBMC_hg38/baseline/data","DESeq2_normalized_276_clean.RDS")
saveRDS(cmd.out,file=output.filepath,refhook = NULL)
