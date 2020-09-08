####################################
# 0. for each sample, generate the TPM separately
# use miniconda to install TPMcalculator and BAMtools
star.dir<-"/home/xy48/scratch/GRADS/SARC_results/Results_summary_PBMC_hg38/baseline/results_STAR_hg38"
gtf.filepath<-"/home/xy48/scratch_kaminski/public/genomes/Homo_sapiens/UCSC/hg38/Annotation/Archives/archive-2015-08-14-08-18-15/Genes/genes.gtf"
folder.names<-list.files(star.dir)
work.dir<-"/home/xy48/scratch/GRADS/SARC_results/Results_summary_PBMC_hg38/baseline"
source("~/Rprogram/my_functions.R")


output.dir<-file.path(work.dir,"results_TPMcalculator_hg38")
script.dir<-file.path(work.dir,"scripts_TPMcalculator_hg38")

if(file.exists(output.dir)==F){
  dir.create(output.dir)
}

if(file.exists(script.dir)==F){
  dir.create(script.dir)
}

jobsub.filepath<-file.path(script.dir,"jobsub.bat")
if(file.exists(jobsub.filepath)){
  file.remove(jobsub.filepath)
}
file.create(jobsub.filepath)

for(i in 1:length(folder.names)){

  temp.name<-paste(my.element.extract(folder.names[i],splitchar="_",index=1),"_",my.element.extract(folder.names[i],splitchar="_",index=3),sep="")
  
  output.subdir<-file.path(output.dir,temp.name)
  if(file.exists(output.subdir)==F){
    dir.create(output.subdir)
  }
  
  script.filepath<-file.path(script.dir,paste(temp.name,".sh",sep=""))
  bam.filepath<-file.path(star.dir,folder.names[i],"Final.STARBowtie2.bam")
  cmd.out<-"#!/bin/bash\n#SBATCH --partition=day,pi_kaminski\n#SBATCH --ntasks=1 --nodes=1 --cpus-per-task=1\n#SBATCH --mem=49152\n#SBATCH --time=24:00:00\n#SBATCH --mail-type=NONE\n"
  cmd.out<-paste(cmd.out,"#SBATCH --job-name=",temp.name,"\n",sep="")
  cmd.out<-paste(cmd.out,"#SBATCH --error=",script.filepath,".e%J\n",sep="")
  cmd.out<-paste(cmd.out,"#SBATCH --output=",script.filepath,".o%J\n",sep="")
  cmd.out<-paste(cmd.out,"module resotre tmpcalculator\n",sep="")
  cmd.out<-paste(cmd.out,"conda activate rnaseq\n",sep="")
  cmd.out<-paste(cmd.out,"cd ",output.subdir,"\n",sep="")
  cmd.out<-paste(cmd.out,"TPMCalculator -g ",gtf.filepath," -b ",bam.filepath,"\n",sep="")
  cat(cmd.out,file=script.filepath,append=F)
  
  cat("sbatch < ",script.filepath,"\n",sep="",file=jobsub.filepath,append=T)
  system(paste("chmod 700 ",script.filepath,"\n",sep=""))
}
system(paste("chmod 700 ",jobsub.filepath,sep=""))











####################################
# 1. merge the results of TPMs into one big matrix from individual samples
source("~/Rprogram/my_functions.R")
data.dir<-"/home/xy48/scratch/GRADS/SARC_results/Results_summary_PBMC_hg38/baseline/results_TPMcalculator_hg38"

filenames<-list.files(data.dir)

###############################
# load in the IDs that need to be switched with each other
idcorrection.matrix<-matrix(c(
"06S7133_005","02S7045_011",
"03S7224_006","05S7101_012",
"03B9320_007","06B9128_013",
"03B9219_008","01B9033_014",
"06B9221_009","04B9250_015",
"03S7071_010","03B9078_016"),ncol=2,byrow=T)

idcorrection.kitid.matrix<-matrix(sapply(idcorrection.matrix,my.element.extract,splitchar="_",index=1),nrow=nrow(idcorrection.matrix),byrow=F)

#sample.names<-sapply(filenames,my.element.remove,splitchar="_",index=-1)
sample.names<-filenames
  
count.table<-numeric()

for(i in 1:length(filenames)){
cat("i=",i,"\n",sep="")
# check if this sample needs to be corrected or not.

if(sample.names[i]%in%as.vector(idcorrection.matrix)){

temp.id<-idcorrection.matrix[idcorrection.matrix[,1]==sample.names[i] | idcorrection.matrix[,2]==sample.names[i],]
temp.id<-setdiff(temp.id,sample.names[i])

data.filepath<-file.path(data.dir,temp.id,"Final.STARBowtie2_genes.out")

if(file.exists(data.filepath)==F){
cat("WARNING: ",data.filepath," does not exist!\n",sep="")
next
}

}else{
data.filepath<-file.path(data.dir,filenames[i],"Final.STARBowtie2_genes.out")
}

temp<-read.table(data.filepath,sep="\t",comment.char="",header=T,as.is=TRUE,check.names=F)

if(i==1){
anno.matrix<-temp[,1:6]
count.table<-cbind(count.table,temp[,7])
gene.names<-paste(anno.matrix[,2],":",anno.matrix[,1],sep="")
}else{
rownames(temp)<-paste(temp[,2],":",temp[,1],sep="")
temp<-temp[gene.names,]
count.table<-cbind(count.table,temp[,7])
}

}

colnames(count.table)<-sample.names
output.filepath<-"/home/xy48/scratch/GRADS/SARC_results/Results_summary_PBMC_hg38/baseline/data/TPM_corrected_all.txt"
cmd.out<-cbind(anno.matrix,count.table)
write.table(cmd.out,file=output.filepath,append=F,sep="\t",row.names=F,col.names=T,quote=F)














####################################
# 2. extract the samples wanted and change the column names into the subject ID instead of the KITID_Barcode
source("~/Rprogram/my_functions.R")
data.filepath<-"/home/xy48/scratch/GRADS/SARC_results/Results_summary_PBMC_hg38/baseline/data/TPM_corrected_all.txt"
previous.fpkm.filepath<-"/home/xy48/scratch/GRADS/SARC_results/Data/SARC_PBMCFPKM_matrix_baseline_276_clean.txt"

masterkey.filepath<-"/home/xy48/scratch/GRADS/SARC_results/Data/GRADS Master Progress Key 6-13-16.xlsx"
output.filepath<-"/home/xy48/scratch/GRADS/SARC_results/Results_summary_PBMC_hg38/baseline/data/TPM_baseline_276_clean_GRADSID.txt"
# load in the master key matrix
library(gdata)
mkey<-read.xls(masterkey.filepath,sheet=3,skip=2,header=T,check.names=F,stringsAsFactors = F)
mkey <- mkey[,1:5]

# load in the previous fpkm matrix to get the list of samples to extract
temp<-read.table(previous.fpkm.filepath,sep="\t",header=T,as.is=TRUE,check.names=F)
sample.list<-colnames(temp)[2:ncol(temp)]

# extract the raw count matrix
count.table<-read.table(data.filepath,sep="\t",comment.char="",header=T,as.is=TRUE,check.names=F)
anno.matrix<-count.table[,1:6]
count.table<-count.table[,7:ncol(count.table)]

my.count.table<-count.table[,sample.list]
my.sample.names<-colnames(my.count.table)
my.kitid<-unname(sapply(my.sample.names,my.element.extract,splitchar="_",index=1))

# convert the KITID into subject IDs
my.key<-mkey[mkey[,2]!="",]
rownames(my.key)<-my.key[,2]
my.subject.id<-my.key[my.kitid,1]

colnames(my.count.table)<-my.subject.id

cmd.out<-cbind(anno.matrix,my.count.table)
write.table(cmd.out,file=output.filepath,append=F,sep="\t",row.names=F,col.names=T,quote=F)






####################################
# 3. extract all samples and change the column names into the subject ID instead of the KITID_Barcode
source("~/Rprogram/my_functions.R")
data.filepath<-"/home/xy48/scratch/GRADS/SARC_results/Results_summary_PBMC_hg38/baseline/data/TPM_corrected_all.txt"
#previous.fpkm.filepath<-"/home/xy48/scratch/GRADS/SARC_results/Data/SARC_PBMCFPKM_matrix_baseline_276_clean.txt"

masterkey.filepath<-"/home/xy48/scratch/GRADS/SARC_results/Data/GRADS Master Progress Key 6-13-16.xlsx"
output.filepath<-"/home/xy48/scratch/GRADS/SARC_results/Results_summary_PBMC_hg38/baseline/data/TPM_baseline_allsamples_GRADSID.txt"
# load in the master key matrix
library(gdata)
mkey<-read.xls(masterkey.filepath,sheet=3,skip=2,header=T,check.names=F,stringsAsFactors = F)
mkey <- mkey[,1:5]

# extract the raw count matrix
count.table<-read.table(data.filepath,sep="\t",comment.char="",header=T,as.is=TRUE,check.names=F)
anno.matrix<-count.table[,1:6]
count.table<-count.table[,7:ncol(count.table)]
# remove the non PBMC samples
count.table<-count.table[,substr(colnames(count.table),3,3)=="S"]

my.count.table<-count.table
my.sample.names<-colnames(my.count.table)
my.kitid<-unname(sapply(my.sample.names,my.element.extract,splitchar="_",index=1))

# convert the KITID into subject IDs
my.key<-mkey[mkey[,2]!="",]
rownames(my.key)<-my.key[,2]
my.subject.id<-my.key[my.kitid,1]

colnames(my.count.table)<-my.subject.id

cmd.out<-cbind(anno.matrix,my.count.table)
write.table(cmd.out,file=output.filepath,append=F,sep="\t",row.names=F,col.names=T,quote=F)

