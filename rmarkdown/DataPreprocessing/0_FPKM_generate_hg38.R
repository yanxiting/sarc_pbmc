################################################
#1. create the FPKM matrix for the tophat2+bowtie2 results
library(gdata)
source("/gpfs/home/fas/kaminski/xy48/Rprogram/my_functions.R")

cufflinks.dir<-"/home/fas/kaminski/xy48/scratch/GRADS/results_cufflinks2_hg38"
output.dir<-"/home/fas/kaminski/xy48/scratch/GRADS"

# get the names of all the cufflinks samples
filenames<-list.files(cufflinks.dir)
filenames<-filenames[grep("_out",filenames)]


anno.matrix<-character()
fpkm.matrix<-numeric()

for(i in 1:length(filenames)){

temp<-read.table(file.path(cufflinks.dir,filenames[i],"genes.fpkm_tracking"),check.names=F,header=TRUE,sep="\t")
rownames(temp)<-paste(as.matrix(temp)[,"tracking_id"],"_",as.matrix(temp)[,"tss_id"],sep="")

if(i==1){
temp.names<-rownames(temp)
anno.matrix<-as.matrix(temp)[,1:9]
}else{
temp<-temp[temp.names,]
}

fpkm.matrix<-cbind(fpkm.matrix,temp[,"FPKM"])
}

sample.names<-unname(sapply(filenames,my.element.extract,splitchar="_out",index=1))

colnames(fpkm.matrix)<-paste(sample.names,"_FPKM",sep="")
rownames(fpkm.matrix)<-temp.names

output.filepath<-file.path(output.dir,"FPKM_matrix_hg38_noCorrection_allSamples_allGenes.txt")
cmd.out<-cbind(anno.matrix,fpkm.matrix)
write.table(cmd.out,file=output.filepath,append=F,sep="\t",row.names=F,col.names=T,quote=F)

####
library(gdata)
source("/gpfs/home/fas/kaminski/xy48/Rprogram/my_functions.R")
output.dir<-"/home/fas/kaminski/xy48/scratch/GRADS"

cmd.out<-read.table("/home/fas/kaminski/xy48/scratch/GRADS/FPKM_matrix_hg38_noCorrection_allSamples_allGenes.txt",sep="\t",header=T,check.names=F,stringsAsFactors = F)

idcorrection.matrix<-matrix(c(
"06S7133_005","02S7045_011",
"03S7224_006","05S7101_012",
"03B9320_007","06B9128_013",
"03B9219_008","01B9033_014",
"06B9221_009","04B9250_015",
"03S7071_010","03B9078_016"),ncol=2,byrow=T)

data.matrix<-cmd.out
data.matrix<-data.matrix[,10:ncol(data.matrix)]
sample.names1<-unname(sapply(colnames(data.matrix),my.element.extract,splitchar="_",index=1))
sample.names2<-unname(sapply(colnames(data.matrix),my.element.extract,splitchar="_",index=3))
colnames(data.matrix)<-paste(sample.names1,"_",sample.names2,sep="")


# switch the IDs
for(i in 1:nrow(idcorrection.matrix)){
temp.matrix<-data.matrix[,idcorrection.matrix[i,]]
temp.matrix<-temp.matrix[,c(2,1)]
data.matrix[,idcorrection.matrix[i,]]<-temp.matrix
}

rownames(data.matrix)<-paste(cmd.out[,1],"_",cmd.out[,6],sep="")
output.filepath<-file.path(output.dir,"FPKM_matrix_hg38_Corrected_allSamples_allGenes.txt")
cat("tracking_id\t",file=output.filepath,append=F)
write.table(data.matrix,file=output.filepath,append=T,sep="\t",row.names=T,col.names=T,quote=F)

###subset for the cleaned BAL data set
fpkm.old<-read.table("/home/fas/kaminski/xy48/scratch/GRADS/SARC_results/Data/SARC_BALFPKM_matrix_209_clean.txt",sep="\t",header=T,check.names=F,stringsAsFactors = F)

BAL.fpkm<-data.matrix[,colnames(fpkm.old)[2:ncol(fpkm.old)]]
output.filepath<-file.path(output.dir,"BAL_FPKM_matrix_hg38_Corrected_allGenes.txt")
cat("tracking_id\t",file=output.filepath,append=F)
write.table(BAL.fpkm,file=output.filepath,append=T,sep="\t",row.names=T,col.names=T,quote=F)













################################################
#1. create the FPKM matrix for the tophat2+bowtie2 results
library(gdata)
source("/gpfs/home/fas/kaminski/xy48/Rprogram/my_functions.R")

cufflinks.dir<-"/home/fas/kaminski/xy48/scratch/GRADS/results_cufflinks2_hg38"
output.dir<-"/home/fas/kaminski/xy48/scratch/GRADS"

# get the names of all the cufflinks samples
filenames<-list.files(cufflinks.dir)
filenames<-filenames[grep("_out",filenames)]


anno.matrix<-character()
fpkm.matrix<-numeric()

for(i in 1:length(filenames)){

temp<-read.table(file.path(cufflinks.dir,filenames[i],"genes.fpkm_tracking"),check.names=F,header=TRUE,sep="\t")
rownames(temp)<-paste(as.matrix(temp)[,"tracking_id"],"_",as.matrix(temp)[,"tss_id"],sep="")

if(i==1){
temp.names<-rownames(temp)
anno.matrix<-as.matrix(temp)[,1:9]
}else{
temp<-temp[temp.names,]
}

fpkm.matrix<-cbind(fpkm.matrix,temp[,"FPKM"])
}

sample.names<-unname(sapply(filenames,my.element.extract,splitchar="_out",index=1))

colnames(fpkm.matrix)<-paste(sample.names,"_FPKM",sep="")
rownames(fpkm.matrix)<-temp.names

output.filepath<-file.path(output.dir,"FPKM_matrix_hg38_noCorrection_allSamples_allGenes.txt")
cmd.out<-cbind(anno.matrix,fpkm.matrix)
write.table(cmd.out,file=output.filepath,append=F,sep="\t",row.names=F,col.names=T,quote=F)

####
library(gdata)
source("/gpfs/home/fas/kaminski/xy48/Rprogram/my_functions.R")
output.dir<-"/home/fas/kaminski/xy48/scratch/GRADS"

cmd.out<-read.table("/home/fas/kaminski/xy48/scratch/GRADS/FPKM_matrix_hg38_noCorrection_allSamples_allGenes.txt",sep="\t",header=T,check.names=F,stringsAsFactors = F)

idcorrection.matrix<-matrix(c(
"06S7133_005","02S7045_011",
"03S7224_006","05S7101_012",
"03B9320_007","06B9128_013",
"03B9219_008","01B9033_014",
"06B9221_009","04B9250_015",
"03S7071_010","03B9078_016"),ncol=2,byrow=T)

data.matrix<-cmd.out
data.matrix<-data.matrix[,10:ncol(data.matrix)]
sample.names1<-unname(sapply(colnames(data.matrix),my.element.extract,splitchar="_",index=1))
sample.names2<-unname(sapply(colnames(data.matrix),my.element.extract,splitchar="_",index=3))
colnames(data.matrix)<-paste(sample.names1,"_",sample.names2,sep="")


# switch the IDs
for(i in 1:nrow(idcorrection.matrix)){
temp.matrix<-data.matrix[,idcorrection.matrix[i,]]
temp.matrix<-temp.matrix[,c(2,1)]
data.matrix[,idcorrection.matrix[i,]]<-temp.matrix
}

rownames(data.matrix)<-paste(cmd.out[,1],"_",cmd.out[,6],sep="")
output.filepath<-file.path(output.dir,"FPKM_matrix_hg38_Corrected_allSamples_allGenes.txt")
cat("tracking_id\t",file=output.filepath,append=F)
write.table(data.matrix,file=output.filepath,append=T,sep="\t",row.names=T,col.names=T,quote=F)

###subset for the cleaned BAL data set
fpkm.old<-read.table("/home/fas/kaminski/xy48/scratch/GRADS/SARC_results/Data/SARC_BALFPKM_matrix_209_clean.txt",sep="\t",header=T,check.names=F,stringsAsFactors = F)

BAL.fpkm<-data.matrix[,colnames(fpkm.old)[2:ncol(fpkm.old)]]
output.filepath<-file.path(output.dir,"BAL_FPKM_matrix_hg38_Corrected_allGenes.txt")
cat("tracking_id\t",file=output.filepath,append=F)
write.table(BAL.fpkm,file=output.filepath,append=T,sep="\t",row.names=T,col.names=T,quote=F)





