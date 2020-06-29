#' @title Calculate Synthetic Variables
#' 
#' @param  P
#' 
#' @export
calculateSyntheticVars <- function(data, mapping, fun="mean") {
  vars <- list()
  anno <- list()
  for (lv in unique(mapping)) {
    vars[[length(vars)+1]] <- apply(data[,which(colnames(data) %in% names(mapping)[which(mapping == lv)]),drop=F], 1, fun)
    anno[[length(anno)+1]] <- list(used=colnames(data)[which(colnames(data) %in% names(mapping)[which(mapping == lv)])],
                                   full=names(mapping)[which(mapping == lv)])
  }
  vars <- do.call(cbind, vars)
  colnames(vars) <- paste("SYN_VAR_", 1:length(unique(mapping)), sep="")
  names(anno) <- colnames(vars)
  
  return(list(data=vars, anno=anno))
}
