##classify the differentially expressed RNAs as up-regulated and down-regulated
##mRNA
rm(list = ls())
options(stringsAsFactors = F)
######Read the differentially expressed mRNA
Diff_mRNA<-read.table("F:/HIRI/DEA/mRNA/Sig_mRNA_result.txt",
                      header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
Up_mRNA<-Diff_mRNA$genes[which(Diff_mRNA$logFC>1)]
down_mRNA<-Diff_mRNA$genes[which(Diff_mRNA$logFC<=-1)]
write.table(Up_mRNA,"F:/HIRI/ceRNA/Up_mRNA.txt",sep="\t",quote=FALSE,row.names = F)
write.table(down_mRNA,"F:/HIRI/ceRNA/down_mRNA.txt",sep="\t",quote=FALSE,row.names = F)

######LncRNA
rm(list = ls())
options(stringsAsFactors = F)

Diff_lncRNA<-read.table("F:/HIRI/DEA/lncRNA/Sig_lncRNA_result.txt",
                      header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
Up_lncRNA<-Diff_lncRNA$genes[which(Diff_lncRNA$logFC>1)]
down_lncRNA<-Diff_lncRNA$genes[which(Diff_lncRNA$logFC<=-1)]
write.table(Up_lncRNA,"F:/HIRI/ceRNA/Up_lncRNA.txt",sep="\t",quote=FALSE,row.names = F)
write.table(down_lncRNA,"F:/HIRI/ceRNA/down_lncRNA.txt",sep="\t",quote=FALSE,row.names = F)

######miRNA
rm(list = ls())
options(stringsAsFactors = F)

Diff_miRNA<-read.table("F:/HIRI/DEA/miRNA/Sig_miRNA_result_hsa.txt",
                        header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
Up_miRNA<-Diff_miRNA$genes[which(Diff_miRNA$logFC>1)]
down_miRNA<-Diff_miRNA$genes[which(Diff_miRNA$logFC<=-1)]
write.table(Up_miRNA,"F:/HIRI/ceRNA/Up_miRNA.txt",sep="\t",quote=FALSE,row.names = F)
write.table(down_miRNA,"F:/HIRI/ceRNA/down_miRNA.txt",sep="\t",quote=FALSE,row.names = F)


