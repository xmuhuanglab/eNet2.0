# Rscript  Step0_Input.R --task='xx' 


library("argparse")
parser <- ArgumentParser()
parser$add_argument("--task", help="The task number.")
argsAll <- parser$parse_args()
load(paste0('args_',argsAll$task,'.RData'))

setwd(dir)
source('./genePeakCorr.R')

library(dplyr)
library("Matrix")
library("igraph")
library(Signac)
library(Seurat)
library(reshape2)
library("SummarizedExperiment")



# load refernce genome
genome <- readRDS(paste0(gdir,'refseq/',genomeArg,'TSSRanges.rds'))

message('-------------- Loading and checking input scATAC-seq data --------------')
# load scATAC-seq data
cre.mat <- read.csv(paste0(inputdir, atac.file))
cre.mat <- as.matrix(cre.mat)
cre.mat <- as(cre.mat, "sparseMatrix")
# check enhancer name
peak <- rownames(cre.mat)
tmp <- grepl('^chr.{1,}-\\d{1,}-\\d{1,}', peak) # match chr*-integer-integer
if(!all(tmp)==T){
  stop('Please check whether the peaks in the scATAC-seq matrix are in the format of "chr-start-end"! ')
}


if(rna.exists=='NO'){
  saveRDS(cre.mat, file = paste0(dir,task,'/out/0.cre.mat.RDS'))
}


if(rna.exists=='YES'){
  message('-------------- Loading and checking input scRNA-seq data --------------')
  # load scRNA-seq data
  rna.mat <- read.csv(paste0(inputdir, rna.file))
  rna.mat <- as.matrix(rna.mat)
  message('Input loading done!')
  # check gene name
  gene <- length(intersect(rownames(rna.mat), as.character(genome$gene_name)))
  if(length(gene)==0)
    stop('Please check whether the genes in the scRNA-seq matrix are in Gene Symbol format')
  
  message('-------------- Checking input data --------------')
  # check cell number
  cells <- intersect(colnames(cre.mat),colnames(rna.mat))
  if(length(cells)==0)
    stop('The columns of the two input matrices are different, please make them identical!')
  if(length(cells)<100)
    warning('Too few cells may result in inaccurate results!')
  # prepare final data
  cre.mat <- cre.mat[,cells]
  rna.mat <- rna.mat[,cells]
  saveRDS(cre.mat, file = paste0(dir,task,'/out/0.cre.mat.RDS'))
  saveRDS(rna.mat, file = paste0(dir,task,'/out/0.rna.mat.RDS'))
  message('Input checking done!')
}


######### Normalize (Optional: YES represents need normalization)
if(atac.raw=='YES'){
  message('-------------- Start normalizing scATAC-seq data --------------')
  cre.mat <- centerCounts(cre.mat)
  saveRDS(cre.mat, file = paste0(dir,task,"/out/0.cre.mat.Normalized.RDS"))
}
if(rna.raw=='YES'){
  message('-------------- Start normalizing scRNA-seq data --------------')
  rna.mat <- NormalizeData(rna.mat, normalization.method = "LogNormalize", scale.factor = 10000) # log1p((x/sum(x))*10000)
  saveRDS(rna.mat, file = paste0(dir,task,"/out/0.rna.mat.Normalized.RDS"))
}



