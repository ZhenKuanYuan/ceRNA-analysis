rm(list=ls())
PPI<-read.table("F:/HIRI/PPI/PPI_result.txt",
                header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
mRNA<-read.table("F:/HIRI/ceRNA/result/UpmRNA_DownmiRNA_UplncRNA_sanggi.txt",
                     header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
mRNA<-mRNA$mRNA
gene_PPI<-PPI[which(PPI$hgnc_symbol.x %in% mRNA & PPI$hgnc_symbol.y %in% mRNA),]

gene_PPI$identical<-apply(gene_PPI, 1, function(x){
  x<-sort(x[4:5]);
  return(paste0(x[1],",",x[2]))})
gene_PPI<-gene_PPI[!duplicated(gene_PPI$identical),]
a<-unique(c(gene_PPI$hgnc_symbol.x,gene_PPI$hgnc_symbol.y))
length(which(mRNA%in%a))
write.table(gene_PPI,"F:/HIRI/PPI/gene_PPI.txt",
            sep="\t",quote=FALSE,row.names = F)
