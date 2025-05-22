####Differential expression analysis
#####mRNA
options(stringsAsFactors = FALSE)
###Clear the environment
rm(list=ls())
####Read the sample label
GSE151648_sample<-read.table("F:/HIRI/original data/GSE151648_sample.txt",
                             header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
####Read the expression profile
mRNA_EXP<-read.table("F:/HIRI/DEA/mRNA/mRNA_count.txt",
                      header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
######Retain patients with IRI- and IRI+ after transplantation
Post_nomal_sample<-GSE151648_sample$samle[which(GSE151648_sample$timepoint %in% " Post-transplant"  & GSE151648_sample$iri%in%" IRI-")]
Post_IRI_sample<-GSE151648_sample$samle[which(GSE151648_sample$timepoint %in% " Post-transplant"  & GSE151648_sample$iri%in%" IRI+")]
########pre/post sample
Post_nomal_sample<-gsub("-","\\.",Post_nomal_sample)
Post_nomal_sample<-paste("X",Post_nomal_sample,sep = "")
Post_IRI_sample<-gsub("-","\\.",Post_IRI_sample)
Post_IRI_sample<-paste("X",Post_IRI_sample,sep = "")
Post_nomal_exp<-mRNA_EXP[,Post_nomal_sample]
colnames(Post_nomal_exp)<-paste("N",1:ncol(Post_nomal_exp),sep = "")
Post_IRI_exp<-mRNA_EXP[,Post_IRI_sample]
colnames(Post_IRI_exp)<-paste("IRI",1:ncol(Post_IRI_exp),sep = "")

gene_exp_result<-data.frame(Post_nomal_exp,Post_IRI_exp)



#BiocManager::install("edgeR")
#install.packages("statmod")
library(edgeR)
library(statmod)
y<-DGEList(counts = gene_exp_result,genes = rownames(gene_exp_result))

#####Filtration and standardization
###Screen out the transcript with the highest count
o <- order(rowSums(y$counts), decreasing=TRUE)
y <- y[o,]
y<-DGEList(counts = y$counts,genes =y$genes)
####TMM standardization
y <- calcNormFactors(y)
###Draw MDS plot
###Check the distance between the samples and the coefficient of biological difference
plotMDS(y)
####Construct the design matrix
group<- factor(c(rep("control",length(Post_nomal_sample)), rep("treat",length(Post_IRI_sample))))
design <- model.matrix(~group)
#Estimate the data dispersion
y <- estimateDisp(y, design, robust=TRUE)
plotBCV(y)
###DEA
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
topTags(lrt)

result<-cbind(y$genes,lrt$table,FDR)
result = result[order(result$PValue),]
write.table(result,"F:/HIRI/DEA/mRNA/mRNA_result.txt",sep="\t",quote=FALSE,row.names=F)
sig_result<-result[result$PValue<0.05&abs(result$logFC)>1,]
index<-which(rownames(gene_exp_result)%in%sig_result$genes)
sig_result<-cbind(sig_result,gene_exp_result[index,])
write.table(sig_result,"F:/HIRI/DEA/mRNA/Sig_mRNA_result.txt",sep="\t",quote=FALSE,row.names=F)


#####lncRNA
options(stringsAsFactors = FALSE)

rm(list=ls())

GSE151648_sample<-read.table("F:/HIRI/GSE151648_sample.txt",
                             header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)

lncRNA_EXP<-read.table("F:/HIRI/DEA/lncRNA/lncRNA_count.txt",
                     header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)

Post_nomal_sample<-GSE151648_sample$samle[which(GSE151648_sample$timepoint %in% " Post-transplant"  & GSE151648_sample$iri%in%" IRI-")]
Post_IRI_sample<-GSE151648_sample$samle[which(GSE151648_sample$timepoint %in% " Post-transplant"  & GSE151648_sample$iri%in%" IRI+")]

Post_nomal_sample<-gsub("-","\\.",Post_nomal_sample)
Post_nomal_sample<-paste("X",Post_nomal_sample,sep = "")
Post_IRI_sample<-gsub("-","\\.",Post_IRI_sample)
Post_IRI_sample<-paste("X",Post_IRI_sample,sep = "")
Post_nomal_exp<-lncRNA_EXP[,Post_nomal_sample]
colnames(Post_nomal_exp)<-paste("N",1:ncol(Post_nomal_exp),sep = "")
Post_IRI_exp<-lncRNA_EXP[,Post_IRI_sample]
colnames(Post_IRI_exp)<-paste("IRI",1:ncol(Post_IRI_exp),sep = "")

gene_exp_result<-data.frame(Post_nomal_exp,Post_IRI_exp)

y<-DGEList(counts = gene_exp_result,genes = rownames(gene_exp_result))


o <- order(rowSums(y$counts), decreasing=TRUE)
y <- y[o,]
y<-DGEList(counts = y$counts,genes =y$genes)

y <- calcNormFactors(y)

group<- factor(c(rep("control",length(Post_nomal_sample)), rep("treat",length(Post_IRI_sample))))
design <- model.matrix(~group)

y <- estimateDisp(y, design, robust=TRUE)

fit <- glmFit(y, design)
lrt <- glmLRT(fit)

result<-cbind(y$genes,lrt$table,FDR)
result = result[order(result$PValue),]
write.table(result,"F:/HIRI/DEA/lncRNA/lncRNA_result.txt",sep="\t",quote=FALSE,row.names=F)
sig_result<-result[result$PValue<0.05&abs(result$logFC)>1,]
index<-which(rownames(gene_exp_result)%in%sig_result$genes)
sig_result<-cbind(sig_result,gene_exp_result[index,])
write.table(sig_result,"F:/HIRI/DEA/lncRNA/Sig_lncRNA_result.txt",sep="\t",quote=FALSE,row.names=F)


#####miRNA
options(stringsAsFactors = FALSE)

rm(list=ls())

GSE151648_sample<-read.table("F:/HIRI/GSE151648_sample.txt",
                             header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)

miRNA_EXP<-read.table("F:/HIRI/DEA/miRNA/miRNA_count.txt",
                     header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)

Post_nomal_sample<-GSE151648_sample$samle[which(GSE151648_sample$timepoint %in% " Post-transplant"  & GSE151648_sample$iri%in%" IRI-")]
Post_IRI_sample<-GSE151648_sample$samle[which(GSE151648_sample$timepoint %in% " Post-transplant"  & GSE151648_sample$iri%in%" IRI+")]

Post_nomal_sample<-gsub("-","\\.",Post_nomal_sample)
Post_nomal_sample<-paste("X",Post_nomal_sample,sep = "")
Post_IRI_sample<-gsub("-","\\.",Post_IRI_sample)
Post_IRI_sample<-paste("X",Post_IRI_sample,sep = "")
Post_nomal_exp<-miRNA_EXP[,Post_nomal_sample]
colnames(Post_nomal_exp)<-paste("N",1:ncol(Post_nomal_exp),sep = "")
Post_IRI_exp<-miRNA_EXP[,Post_IRI_sample]
colnames(Post_IRI_exp)<-paste("IRI",1:ncol(Post_IRI_exp),sep = "")

gene_exp_result<-data.frame(Post_nomal_exp,Post_IRI_exp)

y<-DGEList(counts = gene_exp_result,genes = rownames(gene_exp_result))


o <- order(rowSums(y$counts), decreasing=TRUE)
y <- y[o,]
y<-DGEList(counts = y$counts,genes =y$genes)

y <- calcNormFactors(y)

group<- factor(c(rep("control",length(Post_nomal_sample)), rep("treat",length(Post_IRI_sample))))
design <- model.matrix(~group)

y <- estimateDisp(y, design, robust=TRUE)

fit <- glmFit(y, design)
lrt <- glmLRT(fit)


result<-cbind(y$genes,lrt$table,FDR)
result = result[order(result$PValue),]
write.table(result,"F:/HIRI/DEA/miRNA/miRNA_result.txt",sep="\t",quote=FALSE,row.names=F)
sig_result<-result[result$PValue<0.05&abs(result$logFC)>1,]
index<-which(rownames(gene_exp_result)%in%sig_result$genes)
sig_result<-cbind(sig_result,gene_exp_result[index,])
write.table(sig_result,"F:/HIRI/DEA/miRNA/Sig_miRNA_result.txt",sep="\t",quote=FALSE,row.names=F)

