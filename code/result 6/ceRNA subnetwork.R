
options(stringsAsFactors = FALSE)

rm(list=ls())

ceRNA<-read.table("F:/HIRI/ceRNA/result/UpmRNA_DownmiRNA_UplncRNA_sanggi.txt",
                  header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
ceRNA<-unique(ceRNA)
Imm_gene<-c("IL24","PLAU","CCR5","FGF5")
ceRNA_IMM<-ceRNA[ceRNA$mRNA%in%Imm_gene,]
a<-c(unique(ceRNA_IMM$miRNA),unique(ceRNA_IMM$lncRNA))
ceRNA_IMM<-ceRNA_IMM[,2:3]
ceRNA_IMM<-rbind(ceRNA_IMM,a)

write.table(ceRNA_IMM,"F:/HIRI/ceRNA subnetwork/ceRNA_IMM.txt",sep="\t",quote=FALSE,row.names=F)
