#############GO enrichment analysis
#####Ä£¿é1
rm(list=ls())
library(ggplot2)
library(stringr)
library(clusterProfiler)
#dat <- dat %>% dplyr::group_by(group_type) %>% dplyr::do(head(., n = 13)) 
######cytoscape
sce.markers<-read.table("F:/HIRI/ceRNA/result/UpmRNA_DownmiRNA_UplncRNA_sanggi.txt",
                        header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
sce.markers<-sce.markers$mRNA
######Conversion of gene names and ids
ids = bitr(sce.markers, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
#sce.markers = merge(sce.markers, ids, by.x='V1', by.y='SYMBOL')
#####GO
go<-enrichGO(ids$ENTREZID, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
             qvalueCutoff = 0.2,keyType = 'ENTREZID')
go_result<-as.data.frame(go)
write.table(go_result,"F:/HIRI/UpmRNA function enrichment/GO_result.txt",sep="\t",quote=FALSE,row.names = F)

