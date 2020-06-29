#' @title Find clusters
#' 
#' @description  desc
#' 
#' @import pheatmap
#' 
#' @export
findClusters <- function(data, fun=mad, 
                         nClust=c(1,20),
                         clustering_method="ward.D2",
                         plot=T) {
  ## init cluster
  pm <- pheatmap(data, clustering_method=clustering_method, silent=T)
  
  i <- min(length(data[1,]), nClust[2])
  print(paste("Trying clusters: ", nClust[1], " - ", i))
  print("###################")
  
  ## Store number of clusters + cluster dispersion vals
  metr <- list()
  dispVal <- NULL
  
  repeat {
    print(paste("#Clusters: ", i))
    
    ## Create different numbers of clusters
    pm.clust <- cutree(pm$tree_col, i)
    
    ## Determine variability within clusters
    dv <- list()
    dvSum <- 0
    for (j in 1:i) {
      tmpDispVal <- apply(data[,which(pm.clust == j),drop=F], 1, mad)
      dv[[j]] <- tmpDispVal
      dvSum <- dvSum+sum(tmpDispVal, na.rm=T)
    }
    dvSum <- dvSum/i
    
    ## Compare  total cluster variability to cutoff
    print(paste("Overall variability: ", dvSum))
    
    ##store values
    metr[[length(metr)+1]] <- list(i=i, dispSum=dvSum, disp=dv, pm.clust=pm.clust)
    
    i <- i-1
    if (!(i >= nClust[1])) { break }
  }
  
  metr <- do.call(rbind, metr)

  return(metr)
}
