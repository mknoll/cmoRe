#' @title Assigns Cell-Cycle Status to each cell
#' 
#' @description A
#' 
#' @param data data
#' @param calc Data calculated with cellCycleFractIntegrDNAInt()
#' @param type Method on how to assign cell cycle status
#' to each cell. CIexcl assigns only cells with values 
#' outside the confidence interval for the separating value 
#' determined in cellCycleFractIntegrDNAInt() by bootstrapping.
#' CImedian uses the median as cutoff.
#' 
#' @import parallel
#' @import foreach
#' @import doParallel
#' 
#' @export
assignNcFilterResSimple <- function(obj, calc) {
    nc <- rep(NA, length(obj@data[,1]))
    dat0 <- obj@data$AreaShape_Area.cell/obj@data$AreaShape_Area.nucl

    calcDF <- do.call(rbind, calc)                      
    v <- unlist(calcDF[,"versuch"])                                 
    p <- unlist(calcDF[,"platte"])                      
    for (i in 1:length(calcDF[,1])) {    
	cat(".")               
	w <- which(obj@data$VERSUCH == v[i] & obj@data$PLATTE == p[i])    
	treat <- obj@data$TREATMENT[w]    
	for (tr in unique(unlist(lapply(calc[[i]]$data, function(x) x$treatment)))) {    
	    ww <-which(treat == tr)    
	    val <- lapply(calc[[i]]$data, function(x)  if (x$treatment == tr) {  x$estim })    
	    val <- Filter(length,val)[[1]]                                                             
	    nc[w][ww][which(dat0[w][ww] < val[["minLeftPos"]])] <- 0
	    nc[w][ww][which(dat0[w][ww] >= val[["minLeftPos"]])] <- 1 
	}    
    }
    nc
}
