#' @title Add cutoff variables
#' 
#' @description Cutoff right of the global maximum
#' @import Hmisc
#'
#' @export
addCutoffVars <- function(obj, vars, fun=identity, type="mG", ...) {
    ## Identify first minimum right of global maximum
    if (type != "mG") {
	stop("Not implemented yet!")
    }

    ## Transform data (e.g. log)
    res <- transformVars(obj@data, vars, fun)
    obj@data <- res$data

    ## Find and assign cutoff for each variable
    cutoffData <- list()
    addData <- list()
    addDataNames <- c()
    for (var in res$vars) {
	#### retain only necessary columns
	obj2 <- obj
	### TODO: define requ vars elsewhere
	obj2@data <- obj2@data[,which(colnames(obj2@data) %in% c(var, "TREATMENT", "VERSUCH", "PLATTE"))]

	### get cutoffs
	cu <- cellCycleFractIntegrDNAInt(obj2@data, var=var,xMinGlobMax=min(obj2@data[,var], na.rm=T), ...)
	cu <- imputeMC(cu)
	cutoffData[[length(cutoffData)+1]] <- cu

	### assign single cells 
	newD <- assignCutoffSimple(obj@data, cu, var=var)
	addData[[length(addData)+1]] <- newD 
	addDataNames <- c(addDataNames, paste0("GRP_",var))
	high <- ifelse(newD == "HIGH", 1, 0)
	addData[[length(addData)+1]] <- high
	addDataNames <- c(addDataNames, paste0("GRP_",var,"_BINARY_HIGH"))
	low <- ifelse(newD == "LOW", 1, 0)
	addData[[length(addData)+1]] <- low
	addDataNames <- c(addDataNames, paste0("GRP_",var,"_BINARY_LOW"))

	#obj@data <- assignCutoff(obj@data, cu, var=var, outvar=paste("GRP_", var,sep=""))
    }
    names(cutoffData) <- res$vars
    obj@dataCutoff <- cutoffData

    ## add newly created vars
    addData <- do.call(cbind, addData)
    colnames(addData) <- addDataNames
    obj@data <- cbind(obj@data, addData)

    return(obj)
}



#' @title Create transformed vars
#' @export
transformVars <- function(data, vars, fun=identity) {
    ## TODO check if varsin colnames(data)
    vn <- c()
    add <- list()
    vn <- c()
    for (var in vars) {
        #vn <- c(vn, paste("CUTOFFVAR", var, as.character(substitute(fun)), sep="_"))
        vn <- c(vn, paste("CUTOFFVAR", var, as.character(substitute(fun)), sep="_"))
        add[[length(add)+1]] <- sapply(data[,var], fun)
    }
    add <- do.call(cbind, add)
    colnames(add) <- vn

    return(list(data=cbind(data, add), vars=vn))
}
