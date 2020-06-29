#' @title Plot clusters
#' @description test
#' 
#' @export
plotClusters <- function(data, metr, cut, clustering_method="ward.D2", silent=F) {
  pm <- pheatmap(data, clustering_method=clustering_method, silent=T)
  
  metr <- metr[order(unlist(metr[,"dispSum"])),]

  nCluster <- metr[which(metr[,"dispSum"] >= cut),"i"][1]
  print(paste("nCluster: ", nCluster))
  
  pm.clust <- cutree(pm$tree_col, nCluster)
  df <- data.frame(CLUSTER=factor(pm.clust))
  pheatmap(data, clustering_method=clustering_method, annotation_col=df, silent=silent)
  
  return(nCluster)
}
