
##########################################################################################
# 2. map the reads back to hg38 using STAR+bowtie2 for all the SARC PBMC samples

###########################
# 2.1 generate the files that contains the read length distribution for each sample so that we can get the maximum read length to generate the STAR index file
source("~/Rprogram/my_functions.R")
fastq.dir<-"/home/fas/kaminski/xy48/scratch_kaminski/public/Backup/GRADS_BQ/fastq_JH/GRADS_fastq/fastq"
work.dir<-"/home/fas/kaminski/xy48/scratch/GRADS"
masterkey.filepath<-"/home/fas/kaminski/xy48/scratch_kaminski/xl462/GRADS/datatable/GRADS Master Progress Key 6-13-16.xlsx"
target.dir<-"/home/fas/kaminski/xy48/scratch/GRADS/SARC_results/fastq"
star.filepath<-"/home/fas/kaminski/xy48/scratch_kaminski/public/softwares/STAR-2.6.0c/bin/Linux_x86_64_static/STAR"
bowtie2.filepath<-"/home/fas/kaminski/xy48/scratch_kaminski/public/softwares/bowtie2-2.3.4.1-linux-x86_64/bowtie2"
bowtie.index.dir<-"/home/fas/kaminski/xy48/scratch_kaminski/public/genomes/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome"
picard.filepath<-"/home/fas/kaminski/xy48/scratch_kaminski/public/softwares/picard-2.18.10/bin/picard"
# load in the key matrix file to generate all the sample IDs we want
library(gdata)
key.matrix<-read.xls(masterkey.filepath,sheet=3,check.names=F,as.is=TRUE,header=T,skip=2)

# for A1AT samples
#my.key.matrix<-key.matrix[key.matrix[,5]%in%c("PiMZ","PiZZ Not on Augmenation Therapy","PiZZ on Augmenation Therapy"),]
#my.id<-c(my.key.matrix[,2],my.key.matrix[,4])


# for SARC BAL samples
#my.key.matrix<-key.matrix[!key.matrix[,5]%in%c("PiMZ","PiZZ Not on Augmenation Therapy","PiZZ on Augmenation Therapy"),]
#my.id<-my.key.matrix[,4]

# for SARC PBMC baseline samples
my.key.matrix<-key.matrix[!key.matrix[,5]%in%c("PiMZ","PiZZ Not on Augmenation Therapy","PiZZ on Augmenation Therapy"),]
my.id<-my.key.matrix[,2] # baseline
#my.id<-my.key.matrix[,3]  # 6 month


my.id<-my.id[my.id!=""]

extra.id<-c("06S7133","02S7045",
"03S7224","05S7101",
"03B9320","06B9128",
"03B9219","01B9033",
"06B9221","04B9250",
"03S7071","03B9078")

my.id<-c(my.id,extra.id)
my.id<-unique(my.id)

# load in all the fastq names under the fastq.dir and filter out the ones not in my.id
fastq.filenames<-list.files(fastq.dir)
fastq.id<-unname(sapply(fastq.filenames,my.element.extract,splitchar="_",index=1))
fastq.filenames<-fastq.filenames[fastq.id%in%my.id]

# create the script directory and the output directory
script.dir<-file.path(work.dir,"scripts_readlength")
output.dir<-file.path(work.dir,"results_readlength")


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

for(i in 1:length(fastq.filenames)){
	
kit.id<-my.element.extract(fastq.filenames[i],splitchar="_",index=1)
barcode.id<-my.element.extract(fastq.filenames[i],splitchar="_",index=3)
temp.name<-my.element.extract(fastq.filenames[i],splitchar="\\.",index=1)

script.filepath<-file.path(script.dir,paste(temp.name,".sh",sep=""))
output.filepath<-file.path(output.dir,paste(temp.name,".txt",sep=""))

cmd.out<-"#!/bin/bash\n#SBATCH --partition=day,week,pi_kaminski\n#SBATCH --ntasks=1 --nodes=1 --cpus-per-task=1\n#SBATCH --mem=2152\n#SBATCH --time=24:00:00\n#SBATCH --mail-type=NONE\n"
cmd.out<-paste(cmd.out,"#SBATCH --job-name=",kit.id,"_",barcode.id,"\n",sep="")
cmd.out<-paste(cmd.out,"#SBATCH --error=",script.filepath,".e%J\n",sep="")
cmd.out<-paste(cmd.out,"#SBATCH --output=",script.filepath,".o%J\n",sep="")

# get the maximum read length
cmd.line<-paste("zcat ",file.path(fastq.dir,fastq.filenames[i])," | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c >",output.filepath,sep="")
#temp.text<-system(cmd.line,intern=TRUE)
#max.readlength<-as.numeric(my.element.extract(temp.text[length(temp.text)],splitchar=" ",index=-1))
cmd.out<-paste(cmd.out,cmd.line,"\n",sep="")

cat(cmd.out,file=script.filepath,append=F)

cat("sbatch < ",script.filepath,"\n",sep="",file=jobsub.filepath,append=T)

system(paste("chmod 700 ",script.filepath,sep=""))

}
system(paste("chmod 700 ",jobsub.filepath,sep=""))




###########################
# 1.1 generate the scripts to run STAR+bowtie2 on each sample
source("~/Rprogram/my_functions.R")
fastq.dir<-"/home/fas/kaminski/xy48/scratch_kaminski/public/Backup/GRADS_BQ/fastq_JH/GRADS_fastq/fastq"
work.dir<-"/home/fas/kaminski/xy48/scratch/GRADS"
masterkey.filepath<-"/home/fas/kaminski/xy48/scratch_kaminski/xl462/GRADS/datatable/GRADS Master Progress Key 6-13-16.xlsx"
star.filepath<-"/home/fas/kaminski/xy48/scratch_kaminski/public/softwares/STAR-2.6.0c/bin/Linux_x86_64_static/STAR"
bowtie2.filepath<-"/home/fas/kaminski/xy48/scratch_kaminski/public/softwares/bowtie2-2.3.4.1-linux-x86_64/bowtie2"
bowtie.index.dir<-"/home/fas/kaminski/xy48/scratch_kaminski/public/genomes/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome"
picard.filepath<-"/home/fas/kaminski/xy48/scratch_kaminski/public/softwares/picard-2.18.10/bin/picard"
readlength.dir<-"/home/fas/kaminski/xy48/scratch/GRADS/results_readlength"

# load in the key matrix file to generate all the sample IDs we want
library(gdata)
key.matrix<-read.xls(masterkey.filepath,sheet=3,check.names=F,as.is=TRUE,header=T,skip=2)

# for SARC PBMC baseline samples
my.key.matrix<-key.matrix[!key.matrix[,5]%in%c("PiMZ","PiZZ Not on Augmenation Therapy","PiZZ on Augmenation Therapy"),]
my.id<-my.key.matrix[,2] # baseline

#my.key.matrix<-key.matrix[!key.matrix[,5]%in%c("PiMZ","PiZZ Not on Augmenation Therapy","PiZZ on Augmenation Therapy"),]
#my.id<-my.key.matrix[,4]
my.id<-my.id[my.id!=""]

extra.id<-c("06S7133","02S7045",
"03S7224","05S7101",
"03B9320","06B9128",
"03B9219","01B9033",
"06B9221","04B9250",
"03S7071","03B9078")

my.id<-c(my.id,extra.id)
my.id<-unique(my.id)

# load in all the fastq names under the fastq.dir and filter out the ones not in my.id
fastq.filenames<-list.files(fastq.dir)
fastq.id<-unname(sapply(fastq.filenames,my.element.extract,splitchar="_",index=1))
fastq.filenames<-fastq.filenames[fastq.id%in%my.id]

# create the script directory and the output directory
script.dir<-file.path(work.dir,"scripts_STAR_hg38")
output.dir<-file.path(work.dir,"results_STAR_hg38")


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

for(i in 1:length(fastq.filenames)){
#cat(i,"/",length(fastq.filenames),"\r")
kit.id<-my.element.extract(fastq.filenames[i],splitchar="_",index=1)
barcode.id<-my.element.extract(fastq.filenames[i],splitchar="_",index=3)
temp.name<-my.element.extract(fastq.filenames[i],splitchar="\\.",index=1)

script.filepath<-file.path(script.dir,paste(temp.name,".sh",sep=""))
output.subdir<-file.path(output.dir,paste(temp.name,"_out",sep=""))

if(file.exists(output.subdir)==F){
dir.create(output.subdir)
}

# get the maximum read length
readlength.filepath<-file.path(readlength.dir,paste(temp.name,".txt",sep=""))
temp<-read.table(readlength.filepath)
max.readlength<-temp[nrow(temp),2]

starindex.dir<-file.path(output.subdir,"STARIndex")
if(file.exists(starindex.dir)==T){
temp.names<-list.files(starindex.dir)
if(length(temp.names)>0){
for(j in 1:length(temp.names)){
file.remove(file.path(starindex.dir,temp.names[j]))
}
}
}else{
dir.create(starindex.dir)
}


# generate the header for slurm
cmd.out<-"#!/bin/bash\n#SBATCH --partition=day,week,pi_kaminski\n#SBATCH --ntasks=1 --nodes=1 --cpus-per-task=10\n#SBATCH --mem=49152\n#SBATCH --time=24:00:00\n#SBATCH --mail-type=NONE\n"
cmd.out<-paste(cmd.out,"#SBATCH --job-name=",kit.id,"_",barcode.id,"\n",sep="")
cmd.out<-paste(cmd.out,"#SBATCH --error=",script.filepath,".e%J\n",sep="")
cmd.out<-paste(cmd.out,"#SBATCH --output=",script.filepath,".o%J\n",sep="")

# generate the STAR index file using the maximum read length
cmd.out<-paste(cmd.out,star.filepath," --runMode genomeGenerate --outFileNamePrefix ",output.subdir,"/ --runThreadN 10 --genomeDir ",starindex.dir,"/ --genomeFastaFiles /home/fas/kaminski/xy48/scratch_kaminski/public/genomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa --sjdbGTFfile /home/fas/kaminski/xy48/scratch_kaminski/public/genomes/Homo_sapiens/UCSC/hg38/Annotation/Archives/archive-2015-08-14-08-18-15/Genes/genes.gtf --sjdbOverhang ",max.readlength-1,"\n",sep="")
# map the reads using STAR
#cmd.out<-paste(cmd.out,star.filepath," --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --genomeDir ",starindex.dir," --runThreadN 10 --readFilesIn ",file.path(fastq.dir,fastq.filenames[i])," --outSAMunmapped ",file.path(output.subdir,"unmappedSTAR.sam")," --outReadsUnmapped Fastx  --outFileNamePrefix ",output.subdir,"/ --chimSegmentMin 18 --chimScoreMin 12\n",sep="")
cmd.out<-paste(cmd.out,star.filepath," --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --genomeDir ",starindex.dir," --runThreadN 10 --readFilesIn ",file.path(fastq.dir,fastq.filenames[i])," --outReadsUnmapped Fastx  --outFileNamePrefix ",output.subdir,"/ --chimSegmentMin 18 --chimScoreMin 12\n",sep="")

# map the  unmapped reads using bowtie2
cmd.out<-paste(cmd.out,bowtie2.filepath," --local --very-sensitive-local -p 10 -q --mm -x ",bowtie.index.dir," -U ",file.path(output.subdir,"Unmapped.out.mate1")," --un ",file.path(output.subdir,"sbt2_unmap.fq")," | samtools view -uhS - | samtools sort - ",output.subdir,"/unmapped_remapBowtie2\n",sep="")

# merge the two mapping results
cmd.out<-paste(cmd.out,picard.filepath," MergeSamFiles USE_THREADING=true MSD=true AS=true I=",file.path(output.subdir,"Aligned.sortedByCoord.out.bam")," I=",file.path(output.subdir,"unmapped_remapBowtie2.bam")," O=",file.path(output.subdir,"Final.STARBowtie2.bam"),"\n",sep="")

# remove the STARIndex directory
cmd.out<-paste(cmd.out,"rm ",starindex.dir," -fr\n",sep="")

 
cat(cmd.out,file=script.filepath,append=F)

cat("sbatch < ",script.filepath,"\n",sep="",file=jobsub.filepath,append=T)

system(paste("chmod 700 ",script.filepath,sep=""))

}

system(paste("chmod 700 ",jobsub.filepath,sep=""))






