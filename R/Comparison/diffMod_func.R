matc <- function(gene){
  peaks <- GPPair[[gene]]
  tmp1 <- unique(conns1[which(conns1$Peak1 %in% peaks & conns1$Peak2 %in% peaks),])
  tmp2 <- unique(conns2[which(conns2$Peak1 %in% peaks & conns2$Peak2 %in% peaks),])
  allpair <- unique(c(tmp1$pair, tmp2$pair))
  if(length(allpair)>0){
    df <- data.frame(row.names = allpair, pair=allpair, gene=gene, coaccess1=0, coaccess2=0)
    for(j in rownames(df)){
      tmp1 <- conns1[which(conns1$pair==j),]
      if(nrow(tmp1)>0)
        df[j,'coaccess1'] <- mean(tmp1$coaccess)
      tmp2 <- conns2[which(conns2$pair==j),]
      if(nrow(tmp2)>0)
        df[j,'coaccess2'] <- mean(tmp2$coaccess)
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
