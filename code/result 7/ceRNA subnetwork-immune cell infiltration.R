rm(list = ls())

ceRNA_IMM<-read.table("F:/HIRI/ceRNA subnetwork/ceRNA_IMM.txt",
                         header = T,sep="\t",stringsAsFactors=FALSE,quote = "")


FPKM_exp<-read.table("F:/HIRI/immune/Cibersort/FPKM_exp.txt",
                           header = T,sep="\t",stringsAsFactors=FALSE,quote = "")
ce_FPKM<-FPKM_exp[ceRNA_IMM$mRNA,]

mi_FPKM_exp<-read.table("F:/HIRI/immune/miRNA_FPKM_result.txt",
                     header = T,sep="\t",stringsAsFactors=FALSE,quote = "")
mi_FPKM<-mi_FPKM_exp[rownames(mi_FPKM_exp)%in%"MIR125B1",]
mi_FPKM<-mi_FPKM[,-41]
ce_FPKM_result<-rbind(ce_FPKM,mi_FPKM)
ce_FPKM_result<-t(ce_FPKM_result)


CIBERSORT<-read.table("F:/HIRI/immune/Cibersort/CIBERSORT-Results.txt",
                      header = T,sep="\t",stringsAsFactors=FALSE,quote = "")
CIBERSORT_result<-CIBERSORT[,2:23]
rownames(CIBERSORT_result)<-CIBERSORT$Mixture

a<-apply(CIBERSORT_result,2,sum)
index<-which(a==0)
CIBERSORT_result<-CIBERSORT_result[,-index]
CIBERSORT_result<-as.matrix(CIBERSORT_result)

cor_data<-cbind(ce_FPKM_result,CIBERSORT_result)
library(Hmisc)
cortest <- rcorr(as.matrix(cor_data), type = "spearman")




#library(psych)
#cor_result <- corr.test(ce_FPKM_result,CIBERSORT_result, method = "spearman",adjust="fdr")
cor_r<-cortest$r
cor_p<-cortest$P
cor_p[is.na(cor_p)]<-1
library(corrplot)
color_2<-colorRampPalette(c("CornflowerBlue","White","tomato"))  

corrplot(cor_r,col = color_2(10),method = "color",
         tl.col="black",tl.cex = 0.5,cl.pos = "r",cl.ratio = 0.05,
         p.mat =cor_p, sig.level = c(.001, .01, .05),outline="white",
         insig = "label_sig",pch.cex = 0.5, pch.col = "black")
#####4,5

#####IL24
cor_r<-as.data.frame(cor_r)
cor_p<-as.data.frame(cor_p)

IL24_result<-data.frame(rownames(cor_r),cor_r$IL24,cor_p$IL24)
colnames(IL24_result)<-c("Imm_cell","r","p")
IL24_result<-IL24_result[IL24_result$p<0.05,]
IL24_result$log10P=-log10(IL24_result$p)

library(ggplot2)
library(forcats)

library(RColorBrewer) 
cols<-c("blue","red")
#cols <- rev(brewer.pal(7, 'RdYlBu'))
p<-ggplot(IL24_result, showCategory = 20, 
          aes(r, fct_reorder(Imm_cell, r))) + 
  geom_segment(aes(xend=0, yend = Imm_cell)) +
  geom_point(aes(color=log10P, size = abs(r))) +
  scale_colour_gradientn(colours = cols)+ 
  # scale_color_viridis_c(guide=guide_colorbar(reverse=F),option = "C") +
  # scale_size_continuous(range=c(2, 10)) +
  theme_minimal() + 
  ylab(NULL) 
p
myfile<-"F:/HIRI/immune/IL24_corr.pdf"
ggsave(myfile,p,width = 5,height = 5)

#####PLAU
PLAU_result<-data.frame(rownames(cor_r),cor_r$PLAU,cor_p$PLAU)
colnames(PLAU_result)<-c("Imm_cell","r","p")
PLAU_result<-PLAU_result[PLAU_result$p<0.05,]

PLAU_result$log10P=-log10(PLAU_result$p)
library(ggplot2)
library(forcats)

library(RColorBrewer) 
cols<-c("blue","red")
#cols <- rev(brewer.pal(7, 'RdYlBu'))
p<-ggplot(PLAU_result, showCategory = 20, 
          aes(r, fct_reorder(Imm_cell, r))) + 
  geom_segment(aes(xend=0, yend = Imm_cell)) +
  geom_point(aes(color=log10P, size = abs(r))) +
  scale_colour_gradientn(colours = cols)+ 
  # scale_color_viridis_c(guide=guide_colorbar(reverse=F),option = "C") +
  # scale_size_continuous(range=c(2, 10)) +
  theme_minimal() + 
  ylab(NULL) 
p
myfile<-"F:/HIRI/immune/PLAU_corr.pdf"
ggsave(myfile,p,width = 5,height = 5)

#####CCR5
CCR5_result<-data.frame(rownames(cor_r),cor_r$CCR5,cor_p$CCR5)
colnames(CCR5_result)<-c("Imm_cell","r","p")
CCR5_result<-CCR5_result[CCR5_result$p<0.05,]
CCR5_result$log10P=-log10(CCR5_result$p)
library(ggplot2)
library(forcats)

library(RColorBrewer) 
cols<-c("blue","red")
#cols <- rev(brewer.pal(7, 'RdYlBu'))
p<-ggplot(CCR5_result, showCategory = 20, 
          aes(r, fct_reorder(Imm_cell, r))) + 
  geom_segment(aes(xend=0, yend = Imm_cell)) +
  geom_point(aes(color=log10P, size = abs(r))) +
  scale_colour_gradientn(colours = cols)+ 
  # scale_color_viridis_c(guide=guide_colorbar(reverse=F),option = "C") +
  # scale_size_continuous(range=c(2, 10)) +
  theme_minimal() + 
  ylab(NULL) 
p
myfile<-"F:/HIRI/immune/CCR5_corr.pdf"
ggsave(myfile,p,width = 5,height = 5)

COR_EXP<-data.frame(ce_FPKM_result[,4:6],CIBERSORT_result[,c(4,11,12)])
colnames(COR_EXP)
#######PARD6G.AS1---Monocytes
d <- ggplot(COR_EXP, aes(x = PARD6G.AS1, y = Monocytes))
d<-d + geom_point(color="black")+
  geom_smooth(method = "lm", color = "red")+ 
  scale_color_manual(values = "#00AFBB")+
  stat_cor(method = "spearman", label.x = 1.8, label.y = 0.7) + 
  theme_bw()
d
myfile<-"F:/HIRI/immune/PARD6G.AS1_Monocytes.pdf"
ggsave(myfile,d,width = 3.2,height = 3)

#######FGF5---T.cells.CD8
d <- ggplot(COR_EXP, aes(x = FGF5, y = T.cells.CD8))
d<-d + geom_point(color="black")+
  geom_smooth(method = "lm", color = "red")+ 
  scale_color_manual(values = "#00AFBB")+
  stat_cor(method = "spearman", label.x = 0.4, label.y = 0.1) + 
  theme_bw()
d
myfile<-"F:/HIRI/immune/FGF5_T.cells.CD8.pdf"
ggsave(myfile,d,width = 3.2,height = 3)


#######miRNA---NK
d <- ggplot(COR_EXP, aes(x = MIR125B1, y = NK.cells.activated))
d<-d + geom_point(color="black")+
  geom_smooth(method = "lm", color = "red")+ 
  scale_color_manual(values = "#00AFBB")+
  stat_cor(method = "spearman", label.x = 1.5, label.y = 0.1) + 
  theme_bw()
d
myfile<-"F:/HIRI/immune/MIR125B1_NK.pdf"
ggsave(myfile,d,width = 3.5,height = 3)
