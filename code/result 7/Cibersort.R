rm(list = ls())
BiocManager::install("preprocessCore")
library(e1071)
library(parallel)
library(preprocessCore)
setwd("F:/HIRI/immune/Cibersort/")
source('Cibersort.R')
result1 <- CIBERSORT('LM22.txt','FPKM_exp.txt', perm = 1000, QN = T)

