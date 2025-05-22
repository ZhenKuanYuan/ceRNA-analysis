
options(stringsAsFactors = FALSE)

rm(list=ls())
#####Riaz
hgnc_complete_set<-read.table("F:/HIRI/DEA/miRNA/ID conversion/hgnc_complete_set_result.txt",
                               header=T,sep = "\t", stringsAsFactors=FALSE,quote = "")
####miRNA
miRNA_result<-read.table("F:/HIRI/DEA/miRNA/Sig_miRNA_result.txt",
                              header=T,sep = "\t", stringsAsFactors=FALSE,quote = "")
index<-match(miRNA_result$genes,hgnc_complete_set$symbol)
sigID_miRNA<-data.frame(miRNA_result$genes,hgnc_complete_set$alias_symbol[index])
colnames(sigID_miRNA)<-c("ID","miRNA")
miRNA_result$genes<-hgnc_complete_set$alias_symbol[index]

#miRNA_result$genes<-gsub("hsa-mir-26a-1","hsa-mir-26",miRNA_result$genes)
#miRNA_result$genes<-gsub("\\|hsa-mir-33a","",miRNA_result$genes)
#miRNA_result$genes<-gsub("hsa-mir-125b-1","hsa-mir-125",miRNA_result$genes)
write.table(sigID_miRNA,"F:/HIRI/DEA/miRNA/sigID_miRNA.txt",col.names = T,row.names = F,sep = "\t",quote=F)
write.table(miRNA_result,"F:/HIRI/DEA/miRNA/Sig_miRNA_result_hsa.txt",col.names = T,row.names = F,sep = "\t",quote=F)
