rm(list = ls())
options(stringsAsFactors = F)
######upregulated mRNAs-downregulated miRNAs-upregulated lncRNAs
######upregulated mRNAs and upregulated miRNAs
up_mRNA<-read.table("F:/HIRI/DEA/ceRNA/Up_mRNA.txt",
                       header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)

#####
setwd("F:/HIRI/ceRNA/starbase_DIFFmiRNA_mRNA/DOWN_miRNA/")
Down_miRNA_files<-list.files()
j=1
####
DownmiRNA_UpmRNA<-NULL
j=1
for (j in 1:length(Down_miRNA_files)) {
  miRNA_mRNA<-read.table(Down_miRNA_files[j],header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
  DIFF_targetmRNA<-miRNA_mRNA[which(miRNA_mRNA$geneName%in%up_mRNA$x),]
  DownmiRNA_UpmRNA<-rbind(DownmiRNA_UpmRNA,DIFF_targetmRNA)
}
DownmiRNA_UpmRNA<-unique(DownmiRNA_UpmRNA[,1:4])
write.table(DownmiRNA_UpmRNA,"F:/HIRI/ceRNA/result/DownmiRNA_UpmRNA.txt",sep="\t",quote=FALSE,row.names = F)


######upregulated lncRNAs and downregulated miRNAs
#rm(list = ls())
options(stringsAsFactors = F)
######upregulated lncRNAs and downregulated miRNAs

up_lncRNA<-read.table("F:/HIRI/DEA/ceRNA/Up_lncRNA.txt",
                      header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)

setwd("F:/HIRI/ceRNA/starbase_DIFFmiRNA_lncRNA/DOWN_miRNA/")
Down_miRNA_files<-list.files()
j=1
DownmiRNA_UplncRNA<-NULL
j=1
for (j in 1:length(Down_miRNA_files)) {
  miRNA_lncRNA<-read.table(Down_miRNA_files[j],header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE,fill = T)
  DIFF_targetlncRNA<-miRNA_lncRNA[which(miRNA_lncRNA$geneName%in%up_lncRNA$x),]
  DownmiRNA_UplncRNA<-rbind(DownmiRNA_UplncRNA,DIFF_targetlncRNA)
}
DownmiRNA_UplncRNA<-unique(DownmiRNA_UplncRNA[,1:4])

write.table(DownmiRNA_UplncRNA,"F:/HIRI/ceRNA/result/DownmiRNA_UplncRNA.txt",sep="\t",quote=FALSE,row.names = F)


lncRNA_type<-data.frame(unique(DownmiRNA_UplncRNA$geneName),rep(3,length(unique(DownmiRNA_UplncRNA$geneName))))
colnames(lncRNA_type)<-c("gene","type")

index<-which(DownmiRNA_UpmRNA$miRNAname%in%DownmiRNA_UplncRNA$miRNAname)
DownmiRNA_UpmRNA<-DownmiRNA_UpmRNA[index,]
mRNA_type<-data.frame(unique(DownmiRNA_UpmRNA$geneName),rep(1,length(unique(DownmiRNA_UpmRNA$geneName))))
colnames(mRNA_type)<-c("gene","type")
UpmRNA_DownmiRNA_UplncRNA<-rbind(DownmiRNA_UpmRNA,DownmiRNA_UplncRNA)
UpmRNA_DownmiRNA_UplncRNA<-unique(UpmRNA_DownmiRNA_UplncRNA)
write.table(UpmRNA_DownmiRNA_UplncRNA,"F:/HIRI/ceRNA/result/UpmRNA_DownmiRNA_UplncRNA.txt",sep="\t",quote=FALSE,row.names = F)



