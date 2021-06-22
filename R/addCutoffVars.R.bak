#' @title Add cutoff variables
#' 
#' @description Cutoff right of the global maximum
#' @import Hmisc
#'
#' @export
addCutoffVars <- function(obj, vars, fun=identity, type="mG", plotOnly=F, ...) {
    if (!plotOnly) {
        ## Identify first minimum right of global maximum
        if (type != "mG") {
            stop("Not implemented yet!")
        }
        ## Transform data (e.g. log)
        res <- transformVars(obj@data, vars, fun)
        obj@data <- res$data

        ## Find and assign cutoff for each variable
        cutoffData <- list()
        for (var in res$vars) {
	    #### retain only necessary columns
	    obj2 <- obj
	    ### TODO: define requ vars elsewhere
	    obj2@data <- obj2@data[,which(colnames(obj2@data) %in% c(var, "TREATMENT", "VERSUCH", "PLATTE"))]
	    ####

            cu <- cellCycleFractIntegrDNAInt(obj2@data, var=var,xMinGlobMax=min(obj2@data[,var], na.rm=T), ...)
	    cu <- imputeMC(cu)
            cutoffData[[length(cutoffData)+1]] <- cu
            obj@data <- assignCutoff(obj@data, cu, var=var, outvar=paste("GRP_", var,sep=""))
        }
        names(cutoffData) <- res$vars
        obj@dataCutoff <- cutoffData
        plotV <- paste("GRP_", res$vars,sep="")
    } else {
        plotV <- colnames(obj@data)[which(grepl("GRP_", colnames(obj@data)) & !grepl("BINARY_", colnames(obj@data)) )]
        plotV <- plotV[which(substr(plotV, nchar(plotV)-4, nchar(plotV)) != "Width")]
    }

    ## Create plots for single cell data
    for (var in plotV) {
        ## Non assigned 
        tbl <- table(obj@data[,var], obj@data$TREATMENT, exclude=F)
        barplot(apply(tbl, 2, function(x) x/sum(x)), las=2, main=paste(var, "Non-assigned"))
        ## Assigned
        tbl <- table(obj@data[,var], obj@data$TREATMENT)
        barplot(apply(tbl, 2, function(x) x/sum(x)), las=2, main=var)
        ### TODO: Wilson confidence intervals
    }

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
