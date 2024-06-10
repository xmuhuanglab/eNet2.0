# Rscript  Step1_Edge.R --task='xx' 

library("argparse")
parser <- ArgumentParser()
parser$add_argument("--task", help="The task number.")
argsAll <- parser$parse_args()
load(paste0('args_',argsAll$task,'.RData'))

setwd(dir)
source('./genePeakCorr.R')

library(dplyr)
library("igraph")
library("monocle3")
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(cicero)
library(reshape2)



# read in matrix data 
cre.mat <- readRDS(paste0(dir,task,'/out/0.cre.mat.RDS'))
if(ncol(cre.mat)>10000){
  cols <- ncol(cre.mat)
  k1 <- round(cols/3)
  k2 <- k1*2
  indata <- cbind(as.matrix(cre.mat[,1:k1]), as.matrix(cre.mat[,(k1+1):k2]), as.matrix(cre.mat[,(k2+1):cols]))
}
if(ncol(cre.mat)<=10000){
  indata <- as.matrix(cre.mat)
}
# binarize the matrix
indata[indata>0] <- 1

# format cell info
cellinfo <- as.data.frame(colnames(indata), row.names = colnames(indata))
names(cellinfo) <- "cells"

# format peak info
peakinfo <- colsplit(rownames(indata), pattern = '-', names = c("chr", "bp1", "bp2"))
peakinfo$site_name <- paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep="_")
row.names(peakinfo) <- peakinfo$site_name

rownames(indata) <- row.names(peakinfo)
colnames(indata) <- row.names(cellinfo)

# make CDS
input_cds <-  suppressWarnings(new_cell_data_set(indata,
                                                 cell_metadata = cellinfo,
                                                 gene_metadata = peakinfo))

input_cds <- monocle3::detect_genes(input_cds)

# Ensure there are no peaks included with zero reads
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,] 
# find UMAP coordinates for input_cds
set.seed(2023)
input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method = "LSI")
input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP', 
                              preprocess_method = "LSI")
# saveRDS(input_cds, file = paste0(dir,task,'/out/1.input_cds.RDS'))
# Next, we access the UMAP coordinates from the input CDS object where they are stored by Monocle 3 and run make_cicero_cds:
umap_coords <- reducedDims(input_cds)$UMAP
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)
# To run run_cicero, you need a cicero CDS object (created above) and a genome coordinates file, which contains the lengths of each of the chromosomes in your organism. 
genome <- read.delim(paste0(gdir,'chromsize/',genomeArg,'.chrom.sizes'), header = F)
## run cicero
conns <- run_cicero(cicero_cds, genome, sample_num = 100) 
conns$Peak1 <- gsub('_','-',conns$Peak1)
conns$Peak2 <- gsub('_','-',conns$Peak2)
save(conns, file = paste0(dir,task,'/out/1.conns.Rdata'))



if(rna.exists=='NO'){
  message('-------------- Calculating Cicero gene activity scores --------------')
  # load gtf file
  gene_anno <- readRDS(paste0(gdir,'gtf/gtf_',genomeArg,'.RDS'))
  # rename some columns to match requirements
  gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
  gene_anno$gene <- gene_anno$gene_id
  gene_anno$transcript <- gene_anno$transcript_id
  gene_anno$symbol <- gene_anno$gene_name
  #### Add a column for the pData table indicating the gene if a peak is a promoter ####
  # Create a gene annotation set that only marks the transcription start sites of the genes. We use this as a proxy for promoters.
  # To do this we need the first exon of each transcript
  pos <- subset(gene_anno, strand == "+")
  pos <- pos[order(pos$start),] 
  # remove all but the first exons per transcript
  pos <- pos[!duplicated(pos$transcript),] 
  # make a 1 base pair marker of the TSS
  pos$end <- pos$start + 1 
  
  neg <- subset(gene_anno, strand == "-")
  neg <- neg[order(neg$start, decreasing = TRUE),] 
  # remove all but the first exons per transcript
  neg <- neg[!duplicated(neg$transcript),] 
  neg$start <- neg$end - 1
  
  gene_annotation_sub <- rbind(pos, neg)
  # Make a subset of the TSS annotation columns containing just the coordinates and the gene name
  gene_annotation_sub <- gene_annotation_sub[,c("chromosome", "start", "end", "symbol")]
  # Rename the gene symbol column to "gene"
  names(gene_annotation_sub)[4] <- "gene"
  input_cds <- annotate_cds_by_site(input_cds, gene_annotation_sub)
  
  #### Generate gene activity scores ####
  # generate unnormalized gene activity matrix
  unnorm_ga <- build_gene_activity_matrix(input_cds, conns)
  
  # remove any rows/columns with all zeroes
  unnorm_ga <- unnorm_ga[!Matrix::rowSums(unnorm_ga) == 0, 
                         !Matrix::colSums(unnorm_ga) == 0]
  
  # make a list of num_genes_expressed
  num_genes <- pData(input_cds)$num_genes_expressed
  names(num_genes) <- row.names(pData(input_cds))
  
  # normalize
  ga.mat <- normalize_gene_activities(unnorm_ga, num_genes)
  saveRDS(ga.mat, file = paste0(dir,task,'/out/1.ga.mat.RDS'))
}









