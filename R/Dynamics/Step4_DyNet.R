# Rscript  Step4_DyNet.R --task='222' 

library("argparse")
parser <- ArgumentParser()
parser$add_argument("--task", help="The task number.")
argsAll <- parser$parse_args()
load(paste0('args_',argsAll$task,'.RData'))

setwd(dir)
source('./diffMod_func.R')

library("WGCNA")
library("igraph")
library("doParallel")
library(reshape2)
library("flashClust")
library("matrixStats")


message('-------------- matching peak-peak coaccessibility among all stages  --------------')
load(paste0(dir,task,'/out/2.GPPair.Rdata'))
stages <- strsplit(stage, split = ' ')[[1]]
for(n in 1:length(stages)){
  i <- stages[n]
  load(paste0(dir,task,'/out/1.',i,'.conns.Rdata'))
  conns$Peak1 <- as.character(conns$Peak1)
  conns$Peak2 <- as.character(conns$Peak2)
  conns <- conns[which(conns$Peak1 %in% unlist(GPPair) & conns$Peak2 %in% unlist(GPPair) & conns$coaccess>0),]
  quantile <- quantile(conns$coaccess,seq(0,1,0.05))
  cutoff1 <- round(quantile[[paste0(cutoff,'%')]], 1)
  if(cutoff1==0)
    cutoff1 <- 0.05
  conns[conns$coaccess<cutoff1,'coaccess'] <- 0
  tmp <- do.call(rbind, strsplit(conns$Peak1, split = '-'))
  conns$mid1 <- as.integer(tmp[,2])
  tmp <- do.call(rbind, strsplit(conns$Peak2, split = '-'))
  conns$mid2 <- as.integer(tmp[,2])
  conns$pair <- paste(conns$Peak1, conns$Peak2, sep = ':')
  conns[which(conns$mid1>conns$mid2),]$pair <- paste(conns[which(conns$mid1>conns$mid2),]$Peak2, 
                                                        conns[which(conns$mid1>conns$mid2),]$Peak1, sep = ':')
  tmp <- do.call(rbind, strsplit(conns$pair, split = ':'))
  conns$Peak1 <- tmp[,1]
  conns$Peak2 <- tmp[,2]
  conns <- unique(conns[,-c(4:5)])
  assign(paste0('conns',n),conns)
}

####### cal
getDoParRegistered()
registerDoParallel(nCores)
genes <- names(GPPair)
data <- foreach(g=genes,.inorder=TRUE,
                .errorhandling = 'remove') %dopar% {
                  cat("Running gene: ",g,which(genes == g),"\n")
                  tmp <- matc(gene = g)
                }
closeAllConnections()
data <- do.call(rbind, data)
data$mean <- rowMeans(data[,3:(2+length(stages))])
# save(data, file = paste0(dir,task,'/out/4.coaccess.matched.Rdata'))





message('-------------- Identifing differentially co-expressed modules --------------')
# load(paste0(dir,task,'/out/4.coaccess.matched.Rdata'))
data$peak1 <- colsplit(data$pair, pattern = ':', names = 1:2)[,1]
data$peak2 <- colsplit(data$pair, pattern = ':', names = 1:2)[,2]
stages <- strsplit(stage, split = ' ')[[1]]
data$sum <- 0
for(n in 1:length(stages)){
  data$sum <- data$sum + abs((data[,paste0('coaccess',n)])^2 - (data$mean)^2)
}
data$delta <- sqrt(data$sum/length(stages)) # Compute the adjacency difference matrix
getDoParRegistered()
registerDoParallel(nCores)
genes <- as.character(unique(data$gene))
diffmod <- foreach(g=genes,.inorder=TRUE,
                  .errorhandling = 'remove') %dopar% {
                    cat("Running gene: ",g,which(genes == g),"\n")
                    tmp <- diffCoEx(i = g)
                  }
names(diffmod) <- genes
closeAllConnections()
save(diffmod, file = paste0(dir,task,'/out/4.DyNet.Rdata'))





message('-------------- Calculating network change scores --------------')
stages <- strsplit(stage, split = ' ')[[1]]
load(paste0(dir,task,'/out/2.GPPair.Rdata'))
load(paste0(dir,task,'/out/4.DyNet.Rdata'))
diffNetwork <- data.frame(gene = names(diffmod), diffenhancer = NA)
rownames(diffNetwork) <- diffNetwork$gene
for(n in 1:length(stages)){
  diffNetwork[, paste0('NetworkScore',n)] <- 0
}
for(i in names(diffmod)){
  p <- GPPair[i]
  mod <- diffmod[[i]]
  mod <- mod[mod!=0]
  enh <- names(mod)
  p <- p[[1]][as.integer(enh)]
  if(length(p)!=0){
    diffNetwork[i, 'diffenhancer'] <- paste(p, collapse = ';')
    for(n in 1:length(stages)){
      s <- stages[n]
      load(paste0(dir,task,'/out/3.',s,'.NetworkList.Rdata'))
      inet <- NetworkList[[i]]
      if(length(inet)!=0){
        subnet <- subgraph(inet, enh)
        diffNetwork[i, paste0('NetworkScore',n)] <- length(E(subnet))*2/length(V(subnet)) # network connectivity
      }
    }
  }
} 
diffNetwork$NetworkDynamicsScore <- rowSds(as.matrix(diffNetwork[,3:(2+length(stages))]))

save(diffNetwork, file = paste0(dir,task,'/out/4.NetworkDynamicsInfo.Rdata'))
