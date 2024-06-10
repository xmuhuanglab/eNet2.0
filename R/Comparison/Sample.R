setwd("/cluster/huanglab/xye/project/eNetX/Online/Code/Comparison")
rm(list = ls())

cre.mat <- readRDS("/cluster/huanglab/hlin/project/eNetX/Blood/R3/GPPair/0.cre.mat.rds")
dim(cre.mat) # 177449  70072
rna.mat <- readRDS("/cluster/huanglab/hlin/project/eNetX/Blood/R3/GPPair/0.exp.mat.rds")
dim(rna.mat) # 20287 70072

load("/cluster/huanglab/hlin/project/Database/Work/MPAL/ATAC_RNA_matrix/atac_metadata.Rdata")
metadata$ProjectClassification <- ifelse(metadata$ProjectClassification == "Reference", "healthy", metadata$ProjectClassification)
metadata$ProjectClassification <- ifelse(metadata$ProjectClassification != "healthy", "disease", metadata$ProjectClassification)
metadata$stage <- metadata$ProjectClassification
metadata <- metadata[,c(10,11)]
metadata <- metadata[c(1:500,59501:60000),]
metadata$cell <- gsub('-','.',metadata$cell)
cre.mat <- cre.mat[1:100000,metadata$cell]
rna.mat <- rna.mat[1:10000,metadata$cell]

write.table(cre.mat, file = '0.cre.csv', row.names = T, col.names = T, sep = ',')
write.table(rna.mat, file = '0.rna.csv', row.names = T, col.names = T, sep = ',')
write.table(metadata, file = '0.metadata.csv', row.names = T, col.names = T, sep = ',')
