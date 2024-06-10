# Rscript  Step4_DiffMod.R --task='111' 

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


message('-------------- matching peak-peak coaccessibility between the two stages  --------------')
load(paste0(dir,task,'/out/2.GPPair.Rdata'))
stages <- strsplit(stage, split = ' ')[[1]]
i1 <- stages[1]
i2 <- stages[2]
#### conns1
load(paste0(dir,task,'/out/1.',i1,'.conns.Rdata'))
conns$Peak1 <- as.character(conns$Peak1)
conns$Peak2 <- as.character(conns$Peak2)
conns1 <- conns[which(conns$Peak1 %in% unlist(GPPair) & conns$Peak2 %in% unlist(GPPair) & conns$coaccess>0),]
quantile <- quantile(conns1$coaccess,seq(0,1,0.05))
cutoff1 <- round(quantile[[paste0(cutoff,'%')]], 1)
if(cutoff1==0)
  cutoff1 <- 0.05
conns1[conns1$coaccess<cutoff1,'coaccess'] <- 0
tmp <- do.call(rbind, strsplit(conns1$Peak1, split = '-'))
conns1$mid1 <- as.integer(tmp[,2])
tmp <- do.call(rbind, strsplit(conns1$Peak2, split = '-'))
conns1$mid2 <- as.integer(tmp[,2])
conns1$pair <- paste(conns1$Peak1, conns1$Peak2, sep = ':')
conns1[which(conns1$mid1>conns1$mid2),]$pair <- paste(conns1[which(conns1$mid1>conns1$mid2),]$Peak2, 
                                                      conns1[which(conns1$mid1>conns1$mid2),]$Peak1, sep = ':')
tmp <- do.call(rbind, strsplit(conns1$pair, split = ':'))
conns1$Peak1 <- tmp[,1]
conns1$Peak2 <- tmp[,2]
conns1 <- unique(conns1[,-c(4:5)])

#### conns2
load(paste0(dir,task,'/out/1.',i2,'.conns.Rdata'))
conns$Peak1 <- as.character(conns$Peak1)
conns$Peak2 <- as.character(conns$Peak2)
conns2 <- conns[which(conns$Peak1 %in% unlist(GPPair) & conns$Peak2 %in% unlist(GPPair) & conns$coaccess>0),]
quantile <- quantile(conns2$coaccess,seq(0,1,0.05))
cutoff2 <- round(quantile[[paste0(cutoff,'%')]], 1)
if(cutoff2==0)
  cutoff2 <- 0.05
conns2[conns2$coaccess<cutoff2,'coaccess'] <- 0
tmp <- do.call(rbind, strsplit(conns2$Peak1, split = '-'))
conns2$mid1 <- as.integer(tmp[,2])
tmp <- do.call(rbind, strsplit(conns2$Peak2, split = '-'))
conns2$mid2 <- as.integer(tmp[,2])
conns2$pair <- paste(conns2$Peak1, conns2$Peak2, sep = ':')
conns2[which(conns2$mid1>conns2$mid2),]$pair <- paste(conns2[which(conns2$mid1>conns2$mid2),]$Peak2, 
                                                      conns2[which(conns2$mid1>conns2$mid2),]$Peak1, sep = ':')
tmp <- do.call(rbind, strsplit(conns2$pair, split = ':'))
conns2$Peak1 <- tmp[,1]
conns2$Peak2 <- tmp[,2]
conns2 <- unique(conns2[,-c(4:5)])

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
data$mean <- rowMeans(data[,3:4])
# save(data, file = paste0(dir,task,'/out/4.coaccess.matched.Rdata'))





message('-------------- Identifing differentially co-expressed modules --------------')
# load(paste0(dir,task,'/out/4.coaccess.matched.Rdata'))
data$peak1 <- colsplit(data$pair, pattern = ':', names = 1:2)[,1]
data$peak2 <- colsplit(data$pair, pattern = ':', names = 1:2)[,2]
data$delta <- abs((data$coaccess2)^2-(data$coaccess1)^2) # Compute the adjacency difference matrix

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
save(diffmod, file = paste0(dir,task,'/out/4.DiffNetwork.Rdata'))





message('-------------- Calculating network change scores --------------')
stages <- strsplit(stage, split = ' ')[[1]]
i1 <- stages[1]
i2 <- stages[2]
load(paste0(dir,task,'/out/2.GPPair.Rdata'))
load(paste0(dir,task,'/out/3.',i1,'.NetworkList.Rdata'))
net1 <- NetworkList
load(paste0(dir,task,'/out/3.',i2,'.NetworkList.Rdata'))
net2 <- NetworkList
load(paste0(dir,task,'/out/4.DiffNetwork.Rdata'))
diffNetwork <- data.frame(gene = names(diffmod), diffenhancer = NA, NetworkScore1 = 0, NetworkScore2 = 0, NetworkChangeScore = 0)
rownames(diffNetwork) <- diffNetwork$gene

for(i in names(diffmod)){
  p <- GPPair[i]
  mod <- diffmod[[i]]
  mod <- mod[mod!=0]
  enh <- names(mod)
  p <- p[[1]][as.integer(enh)]
  if(length(p)!=0){
    diffNetwork[i, 'diffenhancer'] <- paste(p, collapse = ';')
    inet1 <- net1[[i]]
    inet2 <- net2[[i]]
    if(length(inet1)!=0){
      subnet1 <- subgraph(inet1, enh)
      diffNetwork[i, 'NetworkScore1'] <- length(E(subnet1))*2/length(V(subnet1)) # network connectivity1
    }
    if(length(inet2)!=0){
      subnet2 <- subgraph(inet2, enh)
      diffNetwork[i, 'NetworkScore2'] <- length(E(subnet2))*2/length(V(subnet2)) # network connectivity2
    }
  }
} 
diffNetwork$NetworkChangeScore <- diffNetwork$NetworkScore2 - diffNetwork$NetworkScore1

save(diffNetwork, file = paste0(dir,task,'/out/4.NetworkChangeInfo.Rdata'))
