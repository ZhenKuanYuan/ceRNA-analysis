
options(stringsAsFactors = FALSE)

rm(list=ls())
FPKM_ID<-read.table("F:/HIRI/immune/count2FPKM.txt",
                  header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
Gencode<-read.table("F:/HIRI/Gencode/Gencode_result.txt",
                    header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
index<-match(rownames(FPKM_ID),Gencode$gene_id)
FPKM_ID$gene_name<-Gencode$gene_name[index]
FPKM_ID<-na.omit(FPKM_ID)
attach(FPKM_ID)
FPKM_ID<-aggregate(FPKM_ID,by=list(gene_name),FUN = mean)
rownames(FPKM_ID)<-FPKM_ID[,1]
FPKM_ID<-FPKM_ID[,-c(1,82)]

GSE151648_sample<-read.table("F:/HIRI/GSE151648_sample.txt",
                             header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)

Post_nomal_sample<-GSE151648_sample$samle[which(GSE151648_sample$timepoint %in% " Post-transplant"  & GSE151648_sample$iri%in%" IRI-")]
Post_IRI_sample<-GSE151648_sample$samle[which(GSE151648_sample$timepoint %in% " Post-transplant"  & GSE151648_sample$iri%in%" IRI+")]

Post_nomal_sample<-gsub("-","\\.",Post_nomal_sample)
Post_nomal_sample<-paste("X",Post_nomal_sample,sep = "")
Post_IRI_sample<-gsub("-","\\.",Post_IRI_sample)
Post_IRI_sample<-paste("X",Post_IRI_sample,sep = "")
Post_nomal_exp<-FPKM_ID[,Post_nomal_sample]
colnames(Post_nomal_exp)<-paste("N",1:ncol(Post_nomal_exp),sep = "")
Post_IRI_exp<-FPKM_ID[,Post_IRI_sample]
colnames(Post_IRI_exp)<-paste("IRI",1:ncol(Post_IRI_exp),sep = "")
gene_exp_result<-data.frame(Post_nomal_exp,Post_IRI_exp)
write.table(gene_exp_result,"F:/HIRI/immune/Cibersort/FPKM_exp.txt",sep="\t",quote=FALSE,row.names=T)
