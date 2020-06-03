library(Rsubread)
library(edgeR)
source("~/Rprogram/my_functions.R")

#output.dir<-"/home/fas/kaminski/xy48/scratch/GRADS/SARC_results/Results_summary_PBMC_baseline_hg38/results_rawcounts_hg38"
#tophat.dir<-"/home/fas/kaminski/xy48/scratch/GRADS/results_STAR_hg38"
#j=1

###############################
# load in the IDs that need to be switched with each other
dir.list<-list.files(tophat.dir)
dir.id<-unname(sapply(dir.list,my.element.extract,splitchar="_",index=1))
dir.barcode<-unname(sapply(dir.list,my.element.extract,splitchar="_",index=3))
dir.samplenames<-paste(dir.id,"_",dir.barcode,sep="")
count.table<-numeric()

result.dirs<-dir.list
#for(j in 1:length(result.dirs)){

kit.id<-my.element.extract(result.dirs[j],splitchar="_",index=1)
barcode.id<-my.element.extract(result.dirs[j],splitchar="_",index=3)
sample.name<-paste(kit.id,"_",barcode.id,sep="")

bam.filepath<-file.path(tophat.dir,result.dirs[j],"Final.STARBowtie2.bam")

temp<-featureCounts(bam.filepath,annot.inbuilt="hg38",annot.ext="/home/ysm/xiting_yan/xy48/scratch_kaminski/public/genomes/Homo_sapiens/UCSC/hg38/Annotation/Archives/archive-2015-08-14-08-18-15/Genes/genes.gtf",isGTFAnnotationFile=T,GTF.featureType="exon",GTF.attrType="gene_id",useMetaFeatures=T,allowMultiOverlap=T,nthreads=1,countMultiMappingReads=T)

anno.matrix<-temp[[2]]
count.table<-cbind(count.table,temp[[1]])

#cell.names<-unname(sapply(result.dirs,my.element.extract,splitchar="_",index=1))
#cell.names2<-unname(sapply(result.dirs,my.element.extract,splitchar="_",index=3))

colnames(count.table)<-sample.name

# output the count.table
#output.filepath<-file.path("..","edgeR",paste(sample.name,"_counttable.txt",sep=""))
#output.filepath<-"/home/fas/kaminski/xy48/scratch/GRADS/SARC_results/Results_summary_PBMC_baseline_hg38/data/rawcounts_corrected_all.txt"
cmd.out<-cbind(anno.matrix,count.table)
output.filepath<-file.path(output.dir,paste(sample.name,"_rawcounts.txt",sep=""))
write.table(cmd.out,file=output.filepath,append=F,sep="\t",row.names=F,col.names=T,quote=F)

# output the annotation matrix
#output.filepath<-file.path("..","edgeR",paste(sample.name,"_annotation_matrix.txt",sep=""))
#output.filepath<-"/home/fas/kaminski/xy48/scratch/GRADS/SARC_results/Results_summary_PBMC_baseline_hg38/data/annotation_matrix_corrected_all.txt"
#write.table(anno.matrix,file=output.filepath,append=F,sep="\t",row.names=F,col.names=T,quote=F)

