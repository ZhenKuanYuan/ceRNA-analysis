
options(stringsAsFactors = FALSE)

rm(list=ls())
library(data.table)

Lnc_Pathways<-fread("F:/HIRI/ceRNA subnetwork/lncRNA/Lnc_Pathways_Sig.txt",
                             header=T,sep = "\t")

ceRNA<-read.table("F:/HIRI/ceRNA/result/UpmRNA_DownmiRNA_UplncRNA_sanggi.txt",
                         header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
lncRNA<-unique(ceRNA$lncRNA)
lncRNA
#######PARD6G-AS1
lncRNA_Path1<-Lnc_Pathways[Lnc_Pathways$lncRNA_symbol%in%"PARD6G-AS1",]
lncRNA_Path1$log10P=-log10(lncRNA_Path1$p_value)

p1<-ggplot(lncRNA_Path1, aes(x=immune_pathway, y=log10P)) + 
  geom_bar(stat = "identity", fill="#CDC0B0") + scale_fill_brewer(palette = "Blues")+
  coord_flip()+ theme_bw() +
  xlab(NULL)+ylab("-log10P")+
  theme(panel.grid=element_blank(), panel.border = element_blank())
ggsave("F:/HIRI/ceRNA subnetwork/lncRNA/PARD6G-AS1_pathway.pdf",p1,width = 4,height = 3)


#######SOX2-OT
lncRNA_Path2<-Lnc_Pathways[Lnc_Pathways$lncRNA_symbol%in%"SOX2-OT",]
lncRNA_Path2$log10P=-log10(lncRNA_Path2$p_value)

p2<-ggplot(lncRNA_Path2, aes(x=immune_pathway, y=log10P)) + 
  geom_bar(stat = "identity", fill="#6495ed") + scale_fill_brewer(palette = "Blues")+
  coord_flip()+ theme_bw() +
  xlab(NULL)+ylab("-log10P")+
  theme(panel.grid=element_blank(), panel.border = element_blank())
p2
ggsave("F:/HIRI/ceRNA subnetwork/lncRNA/SOX2-OT_pathway.pdf",p2,width = 4,height = 3)


#######PKP4-AS1
lncRNA_Path3<-Lnc_Pathways[Lnc_Pathways$lncRNA_symbol%in%"PKP4-AS1",]
lncRNA_Path3$log10P=-log10(lncRNA_Path3$p_value)

p3<-ggplot(lncRNA_Path3, aes(x=immune_pathway, y=log10P)) + 
  geom_bar(stat = "identity", fill="#9dcb9d") + scale_fill_brewer(palette = "Blues")+
  coord_flip()+ theme_bw() +
  xlab(NULL)+ylab("-log10P")+
  theme(panel.grid=element_blank(), panel.border = element_blank())
p3
ggsave("F:/HIRI/ceRNA subnetwork/lncRNA/PKP4-AS1_pathway.pdf",p3,width = 4,height = 3)
