# Rscript  Step2_Node.R --task='000' 

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



# load atac and rna matrices
if(atac.raw=='YES')
  cre.mat <- readRDS(paste0(dir,task,'/out/0.cre.mat.Normalized.RDS'))
if(atac.raw=='NO')
  cre.mat <- readRDS(paste0(dir,task,'/out/0.cre.mat.RDS'))
if(rna.exists=='YES'){
  if(rna.raw=='YES'){
    rna.mat <- readRDS(paste0(dir,task,'/out/0.rna.mat.Normalized.RDS'))
  }
  if(rna.raw=='NO'){
    rna.mat <- readRDS(paste0(dir,task,'/out/0.rna.mat.RDS'))
  }
}
if(rna.exists=='NO'){
  rna.mat <- readRDS(paste0(dir,task,'/out/1.ga.mat.RDS'))
}
# load genome
genome <- readRDS(paste0(gdir,'refseq/',genomeArg,'TSSRanges.rds'))
assign(paste0(genomeArg,'TSSRanges'), genome)



message('-------------- Start calculating gene-peak correlation --------------')
# prepare SummarizedExperiment object
cells <- colnames(cre.mat)
colData <- data.frame(row.names = cells, cell=cells)
peaks <- StringToGRanges(rownames(cre.mat))
peaks$Peak <- rownames(cre.mat)
names(peaks) <- rownames(cre.mat)
se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = cre.mat),
  rowRanges = peaks
)

# genePeak overlap
Ov <- genePeakOv(ATAC.se = se,
                 RNAmat = rna.mat,
                 genome = genomeArg,
                 windowPadSize = 100000,
                 proPadSize = 2000)
# save(Ov,file = paste0(dir,task,"/out/1.genePeakOv.Rdata"))

# calculate gene peak correlation
peaks <- StringToGRanges(rownames(cre.mat))
peaks$Peak <- rownames(cre.mat)
names(peaks) <- rownames(cre.mat)
genes <- as.character(unique(Ov$gene_name))
library(doParallel)
getDoParRegistered()
registerDoParallel(nCores) # registerCores
GPTab <- foreach(g=genes,.combine = 'rbind',.inorder=TRUE,
                 .errorhandling = 'remove') %dopar% {
                   cat("Running: ",g,which(genes == g),"\n")
                   Ovd <- Ov %>% filter(gene_name == g)
                   ObsCor   <- PeakGeneCor(ATAC = cre.mat,
                                           RNA = rna.mat,
                                           peakRanges = peaks,
                                           OV = Ovd)

                 }
save(GPTab, file = paste0(dir,task,"/out/2.GenePeakcorr.Rdata"))
closeAllConnections() # closed cores



message('-------------- Keep significant correlations --------------')
GPTabFilt <- GPTab %>% filter(estimate != "NA" & class == "corr" & estimate > 0 & FDR < 0.05)
cat("Keeping max correlation for multi-mapping peaks ..\n")
GPTabFilt <- GPTabFilt %>% group_by(Peak) %>% dplyr::filter(estimate==max(estimate))
save(GPTabFilt, file = paste0(dir,task,"/out/2.genePeakTabFilt.Rdata"))

GPPair <- list()
for(i in unique(GPTabFilt$Gene)){
  peak <- GPTabFilt[which(GPTabFilt$Gene==i),]$Peak
  peak <- peak[order(peak)]
  GPPair[[i]] <- peak
}
save(GPPair, file = paste0(dir,task,'/out/2.GPPair.Rdata'))

