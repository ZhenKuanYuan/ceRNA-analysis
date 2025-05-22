
options(stringsAsFactors = FALSE)

rm(list=ls())
library(data.table)

Mir_Pathways<-fread("F:/HIRI/ceRNA subnetwork/miRNA/miR_pathway_sig.txt",
                    header=T,sep = "\t")
Mir_Pathways<-unique(Mir_Pathways[,3:5])
colnames(Mir_Pathways)

colnames(Mir_Pathways)<-gsub(" ","_",colnames(Mir_Pathways))
attach(Mir_Pathways)
Mir_Pathways<-aggregate(Mir_Pathways,by=list(miRNA_Symbol,Immune_Pathway),FUN=mean)
Mir_Pathways<-Mir_Pathways[,c(1,2,5)]
colnames(Mir_Pathways)[1:2]<-c("miRNA_Symbol","Immune_Pathway")
###ceRNAÖÐµÄDEG
ceRNA<-read.table("F:/HIRI/ceRNA/result/UpmRNA_DownmiRNA_UplncRNA_sanggi.txt",
                  header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
miRNA<-unique(ceRNA$miRNA)
miRNA
#######hsa-miR-125b-5p
miRNA_Path1<-Mir_Pathways[Mir_Pathways$miRNA_Symbol%in%"hsa-miR-125b-5p",]
miRNA_Path1$log10P=-log10(miRNA_Path1$P_Value)


library(ggplot2)
p1=ggplot(miRNA_Path1, aes(x=Immune_Pathway, y=log10P)) +
  geom_segment( aes(x=Immune_Pathway, xend=Immune_Pathway, y=0, yend=log10P), color="#4a708b",size=1) +
  geom_point( color="#4a708b", size=4) +
  theme_light() +
  coord_flip() +
  xlab(NULL)+ylab("-log10P")+
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )+theme(panel.grid=element_blank(), panel.border = element_blank())
ggsave("F:/HIRI/ceRNA subnetwork/miRNA/hsa-miR-125b-5p_pathway.pdf",p1,width = 4,height = 4)


#######hsa-miR-2355-3p
miRNA_Path2<-Mir_Pathways[Mir_Pathways$miRNA_Symbol%in%"hsa-miR-2355-3p",]
miRNA_Path2$log10P=-log10(miRNA_Path2$P_Value)

library(ggplot2)
p2=ggplot(miRNA_Path2, aes(x=Immune_Pathway, y=log10P)) +
  geom_segment( aes(x=Immune_Pathway, xend=Immune_Pathway, y=0, yend=log10P), color="#cd6839",size=1) +
  geom_point( color="#cd6839", size=4) +
  theme_light() +
  coord_flip() +
  xlab(NULL)+ylab("-log10P")+
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )+theme(panel.grid=element_blank(), panel.border = element_blank())
p2
ggsave("F:/HIRI/ceRNA subnetwork/miRNA/hsa-miR-2355-3p_pathway.pdf",p2,width = 4,height = 4)


#######hsa-miR-2355-5p
miRNA_Path3<-Mir_Pathways[Mir_Pathways$miRNA_Symbol%in%"hsa-miR-2355-5p",]
miRNA_Path3$log10P=-log10(miRNA_Path3$P_Value)


library(ggplot2)
p3=ggplot(miRNA_Path3, aes(x=Immune_Pathway, y=log10P)) +
  geom_segment( aes(x=Immune_Pathway, xend=Immune_Pathway, y=0, yend=log10P), color="#cd3278",size=1) +
  geom_point( color="#cd3278", size=4) +
  theme_light() +
  coord_flip() +
  xlab(NULL)+ylab("-log10P")+
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )+theme(panel.grid=element_blank(), panel.border = element_blank())
p3
ggsave("F:/HIRI/ceRNA subnetwork/miRNA/hsa-miR-2355-5p_pathway.pdf",p3,width = 4,height = 4)
