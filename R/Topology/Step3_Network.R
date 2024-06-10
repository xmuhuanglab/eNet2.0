# Rscript  Step3_Network.R --task='000' 


library("argparse")
parser <- ArgumentParser()
parser$add_argument("--task", help="The task number.")
argsAll <- parser$parse_args()
load(paste0('args_',argsAll$task,'.RData'))

setwd(dir)
source('./genePeakCorr.R')

library(dplyr)
library("igraph")
library("doParallel")


# load data
load(paste0(dir,task,'/out/2.GPPair.Rdata'))
load(paste0(dir,task,"/out/1.conns.Rdata")) # conns

message('-------------- Determining co-accessibility cutoff --------------')
tmp <- conns[which(conns$Peak1 %in% unlist(GPPair) & conns$Peak2 %in% unlist(GPPair)),]
quantile <- quantile(tmp$coaccess[which(tmp$coaccess>0)],seq(0,1,0.05),na.rm = T)
cutoff <- round(quantile[[paste0(cutoff,'%')]], 1)
if(cutoff==0)
  cutoff <- 0.05

message(paste0('Using co-accessibility value ', cutoff, ' to build enhancer networks'))
# apply the cutoff and keep significant co-accessible pairs
conns$Peak1 <- as.character(conns$Peak1)
conns$Peak2 <- as.character(conns$Peak2)
conns <- conns %>% filter(Peak1 %in% unlist(GPPair) & Peak2 %in% unlist(GPPair) & coaccess >= cutoff)



library("doParallel", lib.loc="/cluster/apps/anaconda/201903/envs/R-3.6.2/lib/R/library")
getDoParRegistered()
registerDoParallel(nCores) # registerCores
message('-------------- Building enhancer network --------------')
load(paste0(dir,task,"/out/2.genePeakTabFilt.Rdata"))
genes <- unique(GPTabFilt$Gene)
NetworkList <- foreach(g=genes,.inorder=TRUE,
                       .errorhandling = 'remove') %dopar% {
                         cat("Running gene: ",g,which(genes == g),"\n")
                         GPTabFilt_g <- GPTabFilt %>% dplyr::filter(Gene == g)
                         eNet <- BuildNetwork(conns = conns,
                                              GPTab = GPTabFilt_g,
                                              cutoff = cutoff)
                       }
names(NetworkList) <- genes
save(NetworkList, file = paste0(dir,task,"/out/3.NetworkList.Rdata"))



message('-------------- Calculating the complexity of enhancer networks --------------')
NetworkInfo <- foreach(g=genes,.inorder=TRUE,.combine = 'rbind',
                       .errorhandling = 'remove') %dopar% {
                         cat("Running gene: ",g,which(genes == g),"\n")
                         GPTabFilt_g <- GPTabFilt %>% dplyr::filter(Gene == g)
                         Tab <- NetworkComplexity(net = NetworkList[[g]],
                                                  gene = g,
                                                  GPTab = GPTabFilt_g)
                       }
save(NetworkInfo, file = paste0(dir,task,"/out/3.NetworkInfo.Rdata"))
closeAllConnections() # closed cores



