rm(list=ls())
library(data.table)
####miRNA
miRNA_ENSEMBLE<-read.table("F:/HIRI/Gencode/miRNA_ENSEMBLE.txt",
                           header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)

###lncRNA
lncRNA_ENSEMBLE<-read.table("F:/HIRI/Gencode/lncRNA_ENSEMBLE.txt",
                            header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)

###mRNA
mRNA_ENSEMBLE<-read.table("F:/HIRI/Gencode/mRNA_ENSEMBLE.txt",
                          header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)


#####Count expression profile
Count_exp<-read.table("F:/HIRI/original data/GSE151648_liver-iri-counts.txt",
                     header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
####mRNA
mRNA_EXP<-Count_exp
a<-match(mRNA_EXP$Geneid,mRNA_ENSEMBLE$gene_id)
mRNA_EXP$gene<-mRNA_ENSEMBLE$gene_name[a]
mRNA_EXP<-na.omit(mRNA_EXP)

####Take the mean of repeated values
mRNA_EXP<-aggregate(mRNA_EXP,by=list(mRNA_EXP$gene),FUN=mean)
mRNA_EXP<-mRNA_EXP[,-c(2,83)]
#rownames(mRNA_EXP)<-NULL
library(tibble)
mRNA_EXP <- column_to_rownames(mRNA_EXP,"Group.1")
write.table(mRNA_EXP,"F:/HIRI/DEA/mRNA/mRNA_count.txt",sep="\t",quote=FALSE,row.names=T)

####lncRNA
lncRNA_EXP<-Count_exp
a<-match(lncRNA_EXP$Geneid,lncRNA_ENSEMBLE$gene_id)
lncRNA_EXP$gene<-lncRNA_ENSEMBLE$gene_name[a]
lncRNA_EXP<-na.omit(lncRNA_EXP)

####Take the mean of repeated values
lncRNA_EXP<-aggregate(lncRNA_EXP,by=list(lncRNA_EXP$gene),FUN=mean)
lncRNA_EXP<-lncRNA_EXP[,-c(2,83)]
#rownames(lncRNA_EXP)<-NULL
library(tibble)
lncRNA_EXP <- column_to_rownames(lncRNA_EXP,"Group.1")
write.table(lncRNA_EXP,"F:/HIRI/DEA/lncRNA/lncRNA_count.txt",sep="\t",quote=FALSE,row.names=T)

####miRNA
miRNA_EXP<-Count_exp
a<-match(miRNA_EXP$Geneid,miRNA_ENSEMBLE$gene_id)
miRNA_EXP$gene<-miRNA_ENSEMBLE$gene_name[a]
miRNA_EXP<-na.omit(miRNA_EXP)

####Take the mean of repeated values
miRNA_EXP<-aggregate(miRNA_EXP,by=list(miRNA_EXP$gene),FUN=mean)
miRNA_EXP<-miRNA_EXP[,-c(2,83)]
#rownames(miRNA_EXP)<-NULL
library(tibble)
miRNA_EXP <- column_to_rownames(miRNA_EXP,"Group.1")
write.table(miRNA_EXP,"F:/HIRI/DEA/miRNA/miRNA_count.txt",sep="\t",quote=FALSE,row.names=T)
