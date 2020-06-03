####################################
# 0. for each sample, generate the raw count separately
source("/gpfs/home/fas/kaminski/xy48/Rprogram/my_functions.R")

r.filepath<-"/home/fas/kaminski/xy48/Rprogram/GRADS_PBMC/DataPreprocessing/0_rawcounts_hg38_persample.R"
tophat.dir<-"/home/fas/kaminski/xy48/scratch/GRADS/results_STAR_hg38"
work.dir<-"/home/fas/kaminski/xy48/scratch/GRADS"

script.dir<-file.path(work.dir,"scripts_rawcounts_hg38")
output.dir<-file.path(work.dir,"results_rawcounts_hg38")

if(file.exists(script.dir)==F){
dir.create(script.dir)
}

if(file.exists(output.dir)==F){
dir.create(output.dir)
}


jobsub.filepath<-file.path(script.dir,"jobsub.bat")
if(file.exists(jobsub.filepath)){
file.remove(jobsub.filepath)
}
file.create(jobsub.filepath)

filenames<-list.files(tophat.dir)
temp1<-unname(sapply(filenames,my.element.extract,splitchar="_",index=1))
temp2<-unname(sapply(filenames,my.element.extract,splitchar="_",index=3))

sample.names<-paste(temp1,"_",temp2,sep="")

for(i in 1:length(sample.names)){

script.filepath<-file.path(script.dir,paste(sample.names[i],".sh",sep=""))

# generate the header for slurm
cmd.out<-"#!/bin/bash\n#SBATCH --partition=day,pi_kaminski,week\n#SBATCH --ntasks=1 --nodes=1 --cpus-per-task=1\n#SBATCH --mem=8152\n#SBATCH --time=24:00:00\n#SBATCH --mail-type=NONE\n"
cmd.out<-paste(cmd.out,"#SBATCH --job-name=",sample.names[i],"\n",sep="")
cmd.out<-paste(cmd.out,"#SBATCH --error=",script.filepath,".e%J\n",sep="")
cmd.out<-paste(cmd.out,"#SBATCH --output=",script.filepath,".o%J\n",sep="")

cmd.out<-paste(cmd.out,"R --vanilla<<EOF\n",sep="")
cmd.out<-paste(cmd.out,"tophat.dir<-\"",tophat.dir,"\"\n",sep="")
cmd.out<-paste(cmd.out,"output.dir<-\"",output.dir,"\"\n",sep="")
cmd.out<-paste(cmd.out,"j=",i,"\n",sep="")
cmd.out<-paste(cmd.out,"source(\"",r.filepath,"\")\n",sep="")
cmd.out<-paste(cmd.out,"EOF\n",sep="")

cat(cmd.out,file=script.filepath,append=F)
cat("sbatch < ",script.filepath,"\n",sep="",file=jobsub.filepath,append=T)

system(paste("chmod 700 ",script.filepath,sep=""))
}
system(paste("chmod 700 ",jobsub.filepath,sep=""))













####################################
# 1. merge the results of raw counts into one big matrix from individual samples
source("~/Rprogram/my_functions.R")
data.dir<-"/home/fas/kaminski/xy48/scratch/GRADS/results_rawcounts_hg38"

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

sample.names<-sapply(filenames,my.element.remove,splitchar="_",index=-1)

count.table<-numeric()

for(i in 1:length(filenames)){
cat("i=",i,"\n",sep="")
# check if this sample needs to be corrected or not.

if(sample.names[i]%in%as.vector(idcorrection.matrix)){

temp.id<-idcorrection.matrix[idcorrection.matrix[,1]==sample.names[i] | idcorrection.matrix[,2]==sample.names[i],]
temp.id<-setdiff(temp.id,sample.names[i])

data.filepath<-file.path(data.dir,paste(temp.id,"_rawcounts.txt",sep=""))

if(file.exists(data.filepath)==F){
cat("WARNING: ",data.filepath," does not exist!\n",sep="")
next
}

}else{
data.filepath<-file.path(data.dir,filenames[i])
}

temp<-read.table(data.filepath,sep="\t",header=T,as.is=TRUE,check.names=F)

if(i==1){
anno.matrix<-temp[,1:6]
count.table<-cbind(count.table,temp[,7])
gene.names<-anno.matrix[,1]
}else{
rownames(temp)<-temp[,1]
temp<-temp[gene.names,]
count.table<-cbind(count.table,temp[,7])
}

}

colnames(count.table)<-sample.names
output.filepath<-"/home/fas/kaminski/xy48/scratch/GRADS/SARC_results/Results_summary_PBMC_baseline_hg38/data/rawcounts_corrected_all.txt"
cmd.out<-cbind(anno.matrix,count.table)
write.table(cmd.out,file=output.filepath,append=F,sep="\t",row.names=F,col.names=T,quote=F)
















####################################
# 2. extract the samples wanted and change the column names into the subject ID instead of the KITID_Barcode
source("~/Rprogram/my_functions.R")
data.filepath<-"/home/fas/kaminski/xy48/scratch/GRADS/SARC_results/Results_summary_PBMC_baseline_hg38/data/rawcounts_corrected_all.txt"
previous.fpkm.filepath<-"/home/fas/kaminski/xy48/scratch/GRADS/SARC_results/Data/SARC_PBMCFPKM_matrix_baseline_276_clean.txt"

masterkey.filepath<-"/home/fas/kaminski/xy48/scratch/GRADS/SARC_results/Data/GRADS Master Progress Key 6-13-16.xlsx"
output.filepath<-"/home/fas/kaminski/xy48/scratch/GRADS/SARC_results/Results_summary_PBMC_baseline_hg38/data/rawcounts_baseline_276_clean.txt"
# load in the master key matrix
library(gdata)
mkey<-read.xls(masterkey.filepath,sheet=3,skip=2,header=T,check.names=F,stringsAsFactors = F)
mkey <- mkey[,1:5]

# load in the previous fpkm matrix to get the list of samples to extract
temp<-read.table(previous.fpkm.filepath,sep="\t",header=T,as.is=TRUE,check.names=F)
sample.list<-colnames(temp)[2:ncol(temp)]

# extract the raw count matrix
count.table<-read.table(data.filepath,sep="\t",header=T,as.is=TRUE,check.names=F)
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
# 1. these codes generate the raw count matrix for all SARC BAL samples or a given list of SARC BAL samples
library(Rsubread)
library(edgeR)
source("~/Rprogram/my_functions.R")

#masterkey.filepath<-"/home/fas/kaminski/xy48/scratch/GRADS/SARC_results/Data/GRADS Master Progress Key 6-13-16.xlsx"
tophat.dir<-"/home/fas/kaminski/xy48/scratch/GRADS/results_STAR_hg38"

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


###############################
# load in the IDs that need to be switched with each other
dir.list<-list.files(tophat.dir)
dir.id<-unname(sapply(dir.list,my.element.extract,splitchar="_",index=1))
dir.barcode<-unname(sapply(dir.list,my.element.extract,splitchar="_",index=3))
dir.samplenames<-paste(dir.id,"_",dir.barcode,sep="")
count.table<-numeric()

result.dirs<-dir.list
for(j in 1:length(result.dirs)){

kit.id<-my.element.extract(result.dirs[j],splitchar="_",index=1)
barcode.id<-my.element.extract(result.dirs[j],splitchar="_",index=3)
sample.name<-paste(kit.id,"_",barcode.id,sep="")

if(sample.name%in%as.vector(idcorrection.matrix)){
temp.index<-(1:nrow(idcorrection.matrix))[idcorrection.matrix[,1]==sample.name | idcorrection.matrix[,2]==sample.name]
temp.id<-idcorrection.matrix[temp.index,]
temp.id<-temp.id[temp.id!=sample.name]

if(length(temp.id)>1){
cat("WARNING: MORE THAN 1 ID MATCH FOR j=",j,"\n")
}

temp.dir<-result.dirs[dir.samplenames==temp.id]

if(length(temp.dirs)>1){
cat("WARNING: MORE THAN 1 DIR MATCH FOR j=",j,"\n")
}

bam.filepath<-file.path(tophat.dir,temp.dir,"Final.STARBowtie2.bam")

}else{
bam.filepath<-file.path(tophat.dir,result.dirs[j],"Final.STARBowtie2.bam")
}

temp<-featureCounts(bam.filepath,annot.inbuilt="hg38",annot.ext="/home/fas/kaminski/xy48/scratch_kaminski/public/genomes/Homo_sapiens/UCSC/hg38/Annotation/Archives/archive-2015-08-14-08-18-15/Genes/genes.gtf",isGTFAnnotationFile=T,GTF.featureType="exon",GTF.attrType="gene_id",useMetaFeatures=T,allowMultiOverlap=T,nthreads=4,countMultiMappingReads=T)



if(j==1){
anno.matrix<-temp[[2]]
}
count.table<-cbind(count.table,temp[[1]])
}

cell.names<-unname(sapply(result.dirs,my.element.extract,splitchar="_",index=1))
cell.names2<-unname(sapply(result.dirs,my.element.extract,splitchar="_",index=3))

colnames(count.table)<-paste(cell.names,"_",cell.names2,sep="")

# output the count.table
#output.filepath<-file.path("..","edgeR",paste(sample.name,"_counttable.txt",sep=""))
output.filepath<-"/home/fas/kaminski/xy48/scratch/GRADS/SARC_results/Results_summary_PBMC_baseline_hg38/data/rawcounts_corrected_all.txt"
write.table(count.table,file=output.filepath,append=F,sep="\t",row.names=T,col.names=T,quote=F)

# output the annotation matrix
#output.filepath<-file.path("..","edgeR",paste(sample.name,"_annotation_matrix.txt",sep=""))
output.filepath<-"/home/fas/kaminski/xy48/scratch/GRADS/SARC_results/Results_summary_PBMC_baseline_hg38/data/annotation_matrix_corrected_all.txt"
write.table(anno.matrix,file=output.filepath,append=F,sep="\t",row.names=F,col.names=T,quote=F)




















library(gdata)
source("/home/fas/kaminski/xl462/scratch/GRADS/scripts/my_functions.R")
masterkey.filepath<-"/home/fas/kaminski/xl462/scratch/GRADS/datatable/GRADS Master Progress Key 6-13-16.xlsx"
cuff.dir<-"/home/fas/kaminski/xl462/scratch_public/Backup/GRADS_BQ/results_cufflinks2_hg19"
output.dir<-"/home/fas/kaminski/xl462/scratch/GRADS/output"
FPKM.output.dir<-paste(output.dir,"FPKM.output/",sep="")


if(file.exists(output.dir)==F){
dir.create(output.dir)
}


###############################
# load in the IDs that need to be switched with each other
idcorrection.matrix<-matrix(c(
"06S7133_005","02S7045_011",
"03S7224_006","05S7101_012",
"03B9320_007","06B9128_013",
"03B9219_008","01B9033_014",
"06B9221_009","04B9250_015",
"03S7071_010","03B9078_016"),ncol=2,byrow=T)

###############################
# load in the list of GRADS IDs and KIT IDs from the master key file
mkey<-read.xls(masterkey.filepath,sheet=3,skip=2,header=T,check.names=F,stringsAsFactors = F)
mkey <- mkey[,1:5]

# Extract the sarcoidosis GRADS IDs (3rd letter S)
temp<-unname(sapply(mkey[,1],my.element.extract,splitchar="",index=3))
sarc<-mkey[temp=="S",]

# The corresponding BAL Kit ID
BAL<-sarc[,4]
BAL<-BAL[BAL !=""]


###############################
# load in the list of sequencing data from the cufflinks folder
cufflist<-list.files(cuff.dir)
cuffid <- unname(sapply(cufflist, my.element.extract,splitchar="_",index=1))
cuffbarcode<-unname(sapply(cufflist, my.element.extract,splitchar="_",index=3))
cuffnames<-paste(cuffid,"_",cuffbarcode,sep="")
new.cufflist <- cufflist[cuffid %in% BAL]
length(unique(cuffid[cuffid %in% BAL])) # Check the number of samples corresponding to the BAL IDs # 215
BALcuff.dir<-file.path(cuff.dir,new.cufflist)

# genereate the column names for the FPKM matrix
temp.names.1<-unname(sapply(new.cufflist, my.element.extract,splitchar="_",index=1))
temp.names.2<-unname(sapply(new.cufflist, my.element.extract,splitchar="_",index=3))
reaction.names<-paste(temp.names.1,"_",temp.names.2,sep="")


fpkm.matrix<-numeric()
fpkm.matrix.anno<-character()

for (i in 1:length(BALcuff.dir)){

if(reaction.names[i]%in%idcorrection.matrix){
temp.vect<-idcorrection.matrix[apply(idcorrection.matrix==reaction.names[i],1,sum)>0,]	
temp.truenames<-setdiff(temp.vect,reaction.names[i])
data.filepath<-file.path(cuff.dir,cufflist[cuffnames==temp.truenames],"genes.fpkm_tracking")
}else{
data.filepath<-file.path(BALcuff.dir[i],"genes.fpkm_tracking")	
}

temp.matrix<-read.table(data.filepath,sep="\t",header=T,check.names=F, stringsAsFactors=F)
rownames(temp.matrix)<-paste(temp.matrix[,1],"_",temp.matrix$tss_id,sep="")

if(i==1){
gene.names<-rownames(temp.matrix)
fpkm.matrix<-cbind(fpkm.matrix,temp.matrix[,"FPKM"])
fpkm.matrix.anno<-temp.matrix[,1:9]
}else{
# temp.matrix<-temp.matrix[gene.names,]
fpkm.matrix<-cbind(fpkm.matrix,temp.matrix[gene.names,"FPKM"])	
}
}
rownames(fpkm.matrix)<-gene.names
colnames(fpkm.matrix)<-reaction.names 

# output the fpkm matrix and its annotation matrix
output.filepath<-file.path(output.dir,"FPKM_matrix_240.txt")
cat("row_names\t",file=output.filepath,append=F)
write.table(fpkm.matrix,file=output.filepath,sep="\t",col.names=T,row.names=T,append=T,quote=F)

output.filepath<-file.path(output.dir,"FPKM_matrix_240_annotation.txt")
cat("row_names\t",file=output.filepath,append=F)
write.table(fpkm.matrix.anno,file=output.filepath,sep="\t",col.names=T,row.names=T,append=T,quote=F)





