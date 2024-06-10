setwd("/cluster/huanglab/xye/project/eNetX/Online/Code/Dynamics")
rm(list = ls())


cre.mat <- readRDS("/cluster/huanglab/hlin/project/eNetX/R2/Brain/Signac/atac_sparse_matrix1.rds")
dim(cre.mat) # 295326  31304
rna.mat <- readRDS("/cluster/huanglab/hlin/project/eNetX/R2/Brain/Cicero/All_1.ga.mat.RDS")
dim(rna.mat) # 16206 31304

metadata <- readRDS("/cluster/huanglab/hlin/project/eNetX/R2/Brain/Signac/metadata.rds")
metadata$Age <- gsub("pcw16", "stage1",metadata$Age)
metadata$Age <- gsub("pcw20", "stage2",metadata$Age)
metadata$Age <- gsub("pcw21", "stage3",metadata$Age)
metadata$Age <- gsub("pcw24", "stage4",metadata$Age)
table(metadata$Age)
metadata$cell <- rownames(metadata)
metadata$stage <- metadata$Age
metadata <- metadata[,c(29,30)]
metadata <- metadata[c(1:300,5001:5300,20001:20300,31001:31300),]
table(metadata$stage)

cre.mat <- cre.mat[1:100000,metadata$cell]
rna.mat <- rna.mat[1:10000,metadata$cell]
cre.mat <- as.matrix(cre.mat)
rna.mat <- as.matrix(rna.mat)
write.table(cre.mat, file = '0.cre.csv', row.names = T, col.names = T, sep = ',')
write.table(rna.mat, file = '0.rna.csv', row.names = T, col.names = T, sep = ',')
metadata$cell <- gsub('-','.',metadata$cell)
write.table(metadata, file = '0.metadata.csv', row.names = T, col.names = T, sep = ',')
