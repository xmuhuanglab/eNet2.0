
matc <- function(gene){
  peaks <- GPPair[[gene]]
  allpair <- c()
  for(n in 1:length(stages)){
    conns <- get(paste0('conns',n))
    tmp <- unique(conns[which(conns$Peak1 %in% peaks & conns$Peak2 %in% peaks),])
    allpair <- unique(c(tmp$pair, allpair))
  }
  if(length(allpair)>0){
    df <- data.frame(row.names = allpair, pair=allpair, gene=gene)
    for(n in 1:length(stages)){
      df[,paste0('coaccess',n)] <- 0
      conns <- get(paste0('conns',n))
      for(j in rownames(df)){
        tmp <- conns[which(conns$pair==j),]
        if(nrow(tmp)>0)
          df[j,paste0('coaccess',n)] <- mean(tmp$coaccess)
      }
    }
    return(df)
  }
}





diffCoEx <- function(i){
  sub <- data[which(data$gene==i),]
  peaks <- GPPair[[i]] 
  peaks <- peaks[order(peaks)]
  mat <- matrix(0, nrow = length(peaks), ncol = length(peaks))
  rownames(mat) <- colnames(mat) <- peaks
  for(x in rownames(mat)){
    for(y in colnames(mat)){
      tmp1 <- as.integer(strsplit(x, split = '-')[[1]][2])
      tmp2 <- as.integer(strsplit(y, split = '-')[[1]][2])
      if(tmp1<tmp2)
        paste <- paste0(x,':',y)
      if(tmp1>tmp2)
        paste <- paste0(y,':',x)
      if(tmp1==tmp2){
        mat[x,y] <- 0
        next
      }
      if(tmp1!=tmp2 & length(intersect(paste, rownames(sub)))>0){
        mat[x,y] <- sub[paste,'delta']
        mat[y,x] <- sub[paste,'delta']
      }
    }
  }
  dissTOM <- TOMdist(abs(mat))
  geneTree <- flashClust(as.dist(dissTOM), method = "average")
  dynamicModsHybrid <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                     method="hybrid", cutHeight=0.99, deepSplit = T,
                                     pamRespectsDendro = FALSE, minClusterSize = 1
  ) 
  names(dynamicModsHybrid) <- 1:length(peaks)
  return(dynamicModsHybrid)
}
