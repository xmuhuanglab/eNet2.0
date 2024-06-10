# Rscript  Step4_Module.R --task='000' 


library("argparse")
parser <- ArgumentParser()
parser$add_argument("--task", help="The task number.")
argsAll <- parser$parse_args()
load(paste0('args_',argsAll$task,'.RData'))

setwd(dir)
source(paste0(dir,'genePeakCorr.R'))


library(reshape2)
library(igraph)

message('-------------- Detecting modules and classifying enhancer pairs in each network --------------')
load(paste0(dir,task,"/out/2.GPPair.Rdata"))
load(paste0(dir,task,'/out/3.NetworkInfo.Rdata'))
load(paste0(dir,task,'/out/3.NetworkList.Rdata'))

ModEnhInfo <- data.frame(enhancer=NA, gene=NA, module=NA, degree=NA, modsize=NA, normalized_degree=NA, type=NA)
EnhPair <- data.frame(enhancer1=NA, enhancer2=NA, module1=NA, gene=NA, module2=NA, type=NA)

i <- 'SKI'
for(i in rownames(NetworkInfo)){
  if(NetworkInfo[i,'EdgeNum']>0){
    net <- NetworkList[[i]]
    fc <- cluster_louvain(net) 
    mod <- names(sizes(fc)[which(sizes(fc)>1)]) # keep modules consisting of >=2 enhancers
    
    if(length(mod)>1){ # keep networks consisting of >=2 modules
      peaks <- unlist(GPPair[[i]])
      df <- data.frame(enhancer = peaks, gene =i, module=NA, degree=NA, normalized_degree=NA, type='NonHub')
      df$module <- membership(fc)
      for(m in mod){
        v <- rownames(df[which(df$module==m),])
        subnet <- induced_subgraph(net, v)
        df[v,'degree'] <- degree(subnet)
        df[v,'modsize'] <- vcount(subnet)
      }
      rownames(df) <- df$enhancer
      
      mat <- matrix(i, nrow = length(peaks), ncol = length(peaks))
      mat <- reshape2::melt(mat, varnames = c('enhancer1', 'enhancer2'), value.name = 'gene')
      mat <- mat[mat$enhancer1 > mat$enhancer2,]
      mat[,'enhancer1'] <- peaks[as.integer(mat$enhancer1)]
      mat[,'enhancer2'] <- peaks[as.integer(mat$enhancer2)]
      mat[, 'module1'] <- df[mat$enhancer1, 'module']
      mat[, 'module2'] <- df[mat$enhancer2, 'module']
      mat[which(mat$module1 %in% mod & mat$module2 %in% mod),'type'] <- 'Inter'
      mat[which(mat$module1==mat$module2),'type'] <- 'Intra'
      mat[which(!(mat$module1 %in% mod) & !(mat$module2 %in% mod)),'type'] <- 'Outer' 
      mat[which(is.na(mat$type)),'type'] <- 'Outer' 

      ModEnhInfo <- rbind(ModEnhInfo, df)
      EnhPair <- rbind(EnhPair, mat)
    }
  }
}
ModEnhInfo <- ModEnhInfo[-1,]
EnhPair <- EnhPair[-1,]
save(EnhPair, file = paste0(dir,task,"/out/4.EnhPair.Rdata"))

message('-------------- Defining hub enhancers --------------') # quantile parameter, something like [[70]], by default 75
# scale degree: (enhancer degree)/(module size)
ModEnhInfo$normalized_degree <- ModEnhInfo$degree/ModEnhInfo$modsize
# use a cutoff to classify enhancers into Hub and NonHub
tmp <- quantile(ModEnhInfo$normalized_degree, seq(0,1,0.05), na.rm = T)
Hscore <- tmp[[paste0(Hscore,'%')]]
ModEnhInfo$type[which(ModEnhInfo$normalized_degree >= Hscore)] <- "Hub"
save(ModEnhInfo, file = paste0(dir,task,"/out/4.ModEnhInfo.Rdata"))
