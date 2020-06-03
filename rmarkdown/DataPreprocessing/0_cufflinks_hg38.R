# 1_mapping.R step3
# 3. merge the bam files from the two stages and apply cufflinks to estimate the FPKMs

# we have to convert bam file into sam file, add the header ot the sam file, and then convert it back to bam file, sort them and merge them using samtools

source("/gpfs/home/fas/kaminski/xy48/Rprogram/my_functions.R")
work.dir<-"/home/fas/kaminski/xy48/scratch/GRADS"

bam.dir<-file.path(work.dir,"results_STAR_hg38")
#gtf.filepath<-"/gpfs/home/fas/kaminski/bh234/scratch_public/genomes/Homo_sapiens/UCSC/hg19/Annotation/Archives/archive-2013-03-06-11-23-03/Genes/genes.gtf"
gtf.filepath<-"/home/fas/kaminski/xy48/scratch_kaminski/public/genomes/Homo_sapiens/UCSC/hg38/Annotation/Archives/archive-2015-08-14-08-18-15/Genes/genes.gtf"
genomefa.filepath<-"/gpfs/home/fas/kaminski/xy48/scratch_kaminski/public/genomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"
bam.header.filepath<-"/gpfs/home/fas/kaminski/xy48/Rprogram/hg38_header.txt"

output.dir<-file.path(work.dir,"results_cufflinks2_hg38")
script.dir<-file.path(work.dir,"scripts_cufflinks2_hg38")

if(file.exists(output.dir)==F){
dir.create(output.dir)
}

if(file.exists(script.dir)==F){
dir.create(script.dir)
}

filenames<-list.files(bam.dir)
temp1<-unname(sapply(filenames,my.element.extract,splitchar="_",index=1))
temp2<-unname(sapply(filenames,my.element.extract,splitchar="_",index=3))

sample.names<-paste(temp1,"_",temp2,sep="")

jobsub.filepath<-file.path(script.dir,"jobsub.bat")
if(file.exists(jobsub.filepath)){
file.remove(jobsub.filepath)
}
file.create(jobsub.filepath)

for(i in 1:length(sample.names)){

script.filepath<-file.path(script.dir,paste(sample.names[i],".sh",sep=""))
output.subdir<-file.path(output.dir,filenames[i])

final.bam<-file.path(bam.dir,filenames[i],"final_mapping.bam")

# generate the header for slurm
cmd.out<-"#!/bin/bash\n#SBATCH --partition=day,pi_kaminski,week\n#SBATCH --ntasks=1 --nodes=1 --cpus-per-task=10\n#SBATCH --mem=49152\n#SBATCH --time=24:00:00\n#SBATCH --mail-type=NONE\n"
cmd.out<-paste(cmd.out,"#SBATCH --job-name=",sample.names[i],"\n",sep="")
cmd.out<-paste(cmd.out,"#SBATCH --error=",script.filepath,".e%J\n",sep="")
cmd.out<-paste(cmd.out,"#SBATCH --output=",script.filepath,".o%J\n",sep="")

#cmd.out<-paste(cmd.out,"cufflinks -o ",output.subdir," -p 10 -G ",gtf.filepath," -u ",final.bam,"\n",sep="")
cmd.out<-paste(cmd.out,"/home/fas/kaminski/xy48/scratch_kaminski/public/softwares/cufflinks-2.2.1.Linux_x86_64/cufflinks -o ",output.subdir," -p 10 -G ",gtf.filepath," -u ",final.bam,"\n",sep="")

cat(cmd.out,file=script.filepath,append=F)

#cat("bsub < ",script.filepath,"\n",sep="",file=jobsub.filepath,append=T)
#cat(paste("bsub < ",script.filepath,"\n",sep=""),file=jobsub.filepath,append=T)
cat("sbatch < ",script.filepath,"\n",sep="",file=jobsub.filepath,append=T)

system(paste("chmod 700 ",script.filepath,sep=""))


}
system(paste("chmod 700 ",jobsub.filepath,sep=""))





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























############################################################################################################
# 2. run cufflinks on the results by STAR+bowtie2
source("/gpfs/home/fas/kaminski/xy48/Rprogram/my_functions.R")
work.dir<-"/home/fas/kaminski/xy48/scratch/GRADS"

bam.dir<-file.path(work.dir,"results_STAR_hg38")
#gtf.filepath<-"/gpfs/home/fas/kaminski/bh234/scratch_public/genomes/Homo_sapiens/UCSC/hg19/Annotation/Archives/archive-2013-03-06-11-23-03/Genes/genes.gtf"
gtf.filepath<-"/home/fas/kaminski/xy48/scratch_kaminski/public/genomes/Homo_sapiens/UCSC/hg38/Annotation/Archives/archive-2015-08-14-08-18-15/Genes/genes.gtf"
genomefa.filepath<-"/gpfs/home/fas/kaminski/xy48/scratch_kaminski/public/genomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"

output.dir<-file.path(work.dir,"results_cufflinks2_STARBowtie2_hg38")
script.dir<-file.path(work.dir,"scripts_cufflinks2_STARBowtie2_hg38")

if(file.exists(output.dir)==F){
dir.create(output.dir)
}

if(file.exists(script.dir)==F){
dir.create(script.dir)
}

filenames<-list.files(bam.dir)
temp1<-unname(sapply(filenames,my.element.extract,splitchar="_",index=1))
temp2<-unname(sapply(filenames,my.element.extract,splitchar="_",index=3))

sample.names<-paste(temp1,"_",temp2,sep="")

jobsub.filepath<-file.path(script.dir,"jobsub.bat")
if(file.exists(jobsub.filepath)){
file.remove(jobsub.filepath)
}
file.create(jobsub.filepath)

for(i in 1:length(sample.names)){

script.filepath<-file.path(script.dir,paste(sample.names[i],".sh",sep=""))
output.subdir<-file.path(output.dir,filenames[i])

final.bam<-file.path(bam.dir,filenames[i],"Final.STARBowtie2.bam")

# generate the header for slurm
cmd.out<-"#!/bin/bash\n#SBATCH --partition=day,pi_kaminski,week\n#SBATCH --ntasks=1 --nodes=1 --cpus-per-task=10\n#SBATCH --mem=49152\n#SBATCH --time=24:00:00\n#SBATCH --mail-type=NONE\n"
cmd.out<-paste(cmd.out,"#SBATCH --job-name=",sample.names[i],"\n",sep="")
cmd.out<-paste(cmd.out,"#SBATCH --error=",script.filepath,".e%J\n",sep="")
cmd.out<-paste(cmd.out,"#SBATCH --output=",script.filepath,".o%J\n",sep="")

#cmd.out<-paste(cmd.out,"cufflinks -o ",output.subdir," -p 10 -G ",gtf.filepath," -u ",final.bam,"\n",sep="")
cmd.out<-paste(cmd.out,"/home/fas/kaminski/xy48/scratch_kaminski/public/softwares/cufflinks-2.2.1.Linux_x86_64/cufflinks -o ",output.subdir," -p 10 -G ",gtf.filepath," -u ",final.bam,"\n",sep="")

cat(cmd.out,file=script.filepath,append=F)

#cat("bsub < ",script.filepath,"\n",sep="",file=jobsub.filepath,append=T)
#cat(paste("bsub < ",script.filepath,"\n",sep=""),file=jobsub.filepath,append=T)
cat("sbatch < ",script.filepath,"\n",sep="",file=jobsub.filepath,append=T)

system(paste("chmod 700 ",script.filepath,sep=""))


}
system(paste("chmod 700 ",jobsub.filepath,sep=""))






################################################
#1. create the FPKM matrix for the STAR+bowtie2 results
library(gdata)
source("/gpfs/home/fas/kaminski/xy48/Rprogram/my_functions.R")

cufflinks.dir<-"/home/fas/kaminski/xy48/scratch/GRADS/results_cufflinks2_STARBowtie2_hg38"
output.dir<-"/home/fas/kaminski/xy48/scratch/GRADS/SARC_results/Results_summary_BAL_hg38"

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

colnames(fpkm.matrix)<-sample.names
rownames(fpkm.matrix)<-temp.names

output.filepath<-file.path(output.dir,"FPKM_matrix_hg38_noCorrection_allSamples_allGenes.txt")
cmd.out<-cbind(anno.matrix,fpkm.matrix)
write.table(cmd.out,file=output.filepath,append=F,sep="\t",row.names=F,col.names=T,quote=F)

####
library(gdata)
source("/gpfs/home/fas/kaminski/xy48/Rprogram/my_functions.R")
output.dir<-"/home/fas/kaminski/xy48/scratch/GRADS/SARC_results/Results_summary_BAL_hg38"

cmd.out<-read.table("/home/fas/kaminski/xy48/scratch/GRADS/SARC_results/Results_summary_BAL_hg38/FPKM_matrix_hg38_noCorrection_allSamples_allGenes.txt",sep="\t",header=T,check.names=F,stringsAsFactors = F)

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
output.filepath<-file.path(output.dir,"FPKM_matrix_hg38_Corrected_allGenes_BAL.txt")
cat("tracking_id\t",file=output.filepath,append=F)
write.table(BAL.fpkm,file=output.filepath,append=T,sep="\t",row.names=T,col.names=T,quote=F)

temp.out<-cbind(cmd.out[,1:9],BAL.fpkm)
output.filepath<-file.path(output.dir,"FPKM_matrix_hg38_Corrected_allGenes_BAL_withAnnot.txt")
write.table(temp.out,file=output.filepath,append=F,sep="\t",row.names=F,col.names=T,quote=F)







