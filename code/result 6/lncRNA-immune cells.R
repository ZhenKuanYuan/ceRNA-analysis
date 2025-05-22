
options(stringsAsFactors = FALSE)

rm(list=ls())
#library(data.table)

Lnc_cell<-read.table("F:/HIRI/ceRNA subnetwork/lncRNA/Lnc_Immunecell_Sig.txt",
                    header=T,sep = "\t")

ceRNA<-read.table("F:/HIRI/ceRNA/result/UpmRNA_DownmiRNA_UplncRNA_sanggi.txt",
                  header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
lncRNA<-unique(ceRNA$lncRNA)
lncRNA_cell<-Lnc_cell[Lnc_cell$lncRNA_symbol%in%lncRNA,]
lncRNA_cell1<-paste(lncRNA_cell$immune_cell)

lncRNA_cell1<-aggregate(lncRNA_cell,by=list(lncRNA_cell$lncRNA_symbol,lncRNA_cell$immune_cell), FUN=mean)

lncRNA_cell_long<-lncRNA_cell1[,c(1,2,8)]
colnames(lncRNA_cell_long)[1:2]<-c("lncRNA_symbol","immune_cell")
class(lncRNA_cell_long)
#lncRNA_Path1$log10P=-log10(lncRNA_Path1$p_value)

library(tidyr)
lncRNA_cell_wide<-spread(lncRNA_cell_long,immune_cell,cor_R_value) 

#lncRNA_cell_wide<-spread(lncRNA_cell_long,lncRNA_symbol,cor_R_value) 
#install.packages("corrplot")
rownames(lncRNA_cell_wide)<-lncRNA_cell_wide$lncRNA_symbol

lncRNA_cell_wide<-as.matrix(lncRNA_cell_wide[,-1])

#install.packages("ggcorrplot")
library(ggcorrplot)
ggcorrplot(lncRNA_cell_wide,method = "circle",outline.color = "white",
           ggtheme = theme_bw(),colors = c("#839EDB", "white", "#FF8D8D"),lab = T,lab_size=2,
          insig= "blank", pch.col = "red", pch.cex = 3, tl.cex = 12)
####5 3




# Libraries

library(ggplot2)
# Horizontal version
ggplot(lncRNA_Path1, aes(x=immune_pathway, y=log10P)) +
  geom_segment( aes(x=immune_pathway, xend=immune_pathway, y=0, yend=log10P), color="skyblue",size=1.5) +
  geom_point( color="darkblue", size=3, alpha=0.6) +
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )