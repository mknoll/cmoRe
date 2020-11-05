# An S4 Helper class to allow NULL values
setClassUnion("characterOrNULL", c("character","NULL"))
# An S4 Helper class to allow NULL values
setClassUnion("listOrNULL", c("list", "NULL"))
# An S4 Helper class to allow NULL values
setClassUnion("dataframeOrNULL", c("data.frame", "NULL"))

#' An S4 class representing a single experiment
#' 
#' @slot expID Experiment ID
#' @slot experiment list, each element corresponds to a
#' single run; each run contains 1 to n folders as list
#' @slot data data.frame containing the sinlge measurements;
#' altered during filtering steps
#' @slot dataCC result of the cellCycleFractIntegrDNAInt() function,
#' which determines the cell cycle
#' @slot dataNC result of the ncAreaFilter() function, used 
#' to filter cells by cell/nuleus area ratio
#' @slot dataFB result of the cellCycleFractIntegrDNAInt(), used 
#' to filter fibroblasts 
#' @slot dataCutoff TODO
#' @slot dataAgg per well aggregated data
#' @slot ddataZ z-transformed dataAgg data
#' @slot log event-log of the different processing steps
imgExp <- setClass("imageExp", 
                   slots=c(expID="character",
                           experiment="list",
                           data="dataframeOrNULL",
                           dataCC="listOrNULL",
                           dataNC="listOrNULL",
                           dataFB="listOrNULL",
                           dataCutoff="listOrNULL",
                           dataAgg="dataframeOrNULL",
                           dataZ="dataframeOrNULL",
                           log="listOrNULL"
                   ))

#' @title Constructor imageExp
#'
#' @description Constructor for imageExp
#' 
#' @param .Object imageExp object
#' @param expID experiment ID
#' @param experiment two level list containing experimental
#' data folders
#' 
#' @return imageExp instance
setMethod("initialize", "imageExp",
          function(.Object,
                   expID=character,
                   experiment=list) {
              .Object@expID <- expID
              .Object@experiment <- experiment
              .Object@data <- NULL
              .Object@dataCC <- NULL
              .Object@dataNC <- NULL
              .Object@dataCutoff <- NULL
              .Object@dataAgg <- NULL
              .Object@dataZ <- NULL
              ## Log
              l <- list()
              l[[length(l)+1]] <- data.frame(Sys.time(), "Created Instance")
              l[[length(l)+1]] <- data.frame(Sys.time(), paste("On", paste(Sys.info(), collapse=";")))
              l[[length(l)+1]] <- data.frame(Sys.time(), paste("Basedir: ", getwd()))
              .Object@log <- l 
              .Object
          })

#' @title Load experiment data from
#' @description reads cell profiler output and treatment files
#' @param obj imgExp instance
#' @export
getData <- function(obj, force=F, ...) {
    if (is.null(obj@data) || force) {
        print("Reading data ... ")
        l <- obj@log
        l[[length(l)+1]] <- data.frame(Sys.time(), "Started reading data...")
        obj@data <- readData(obj@experiment, ...)
        l[[length(l)+1]] <- data.frame(Sys.time(), "Finished reading data...")
        obj@log <- l
    } else {
        stop("Data already loaded. Exiting. Set force=T to re-read data")
    }
    return(obj)
}

#' @title Determine Cell Cycle, filter for Cell/Nucleus ratio and fibroblasts
#' @description Calculates derived variables
#' @param obj imgExp instance
#' @param fun function to calculate, can be cc (cell cycle) or nc
#' (nucler/cell ratio) or fb (gibrobaslten). defaults to NULL 
#' (all calculations are performed)
#' @export
calc <- function(obj, fun=NULL, varFB="Intensity_MedianIntensity_DNA.nucl", 
				varCC="Intensity_IntegratedIntensity_DNA.nucl",
				xMinGlobMaxFB=0.2, xMaxGlobMaxFB=0.6,
				...) {
    check(obj)
    if (is.null(fun) || fun == "cc") {
	if (!varCC %in% colnames(obj@data)) {
	    warning(paste("CC Variable (", varCC, ") not found!"))
	} else {
	    l <- obj@log
	    l[[length(l)+1]] <- data.frame(Sys.time(), "Started calculating CellCycle ...")
	    l[[length(l)+1]] <- data.frame(Sys.time(), paste("Parameters: ", ..., collapse="|"))
	    obj@dataCC <- cellCycleFractIntegrDNAInt(obj@data, log=T, xMinGlobMax=log(100), var=varCC, ...)
	    l[[length(l)+1]] <- data.frame(Sys.time(), paste("Impute cutoffs ... "))
	    obj@dataCC <- imputeMC(obj@dataCC, ...)
	    l[[length(l)+1]] <- data.frame(Sys.time(), "Assign CellCycle ...")
	    obj@data <- assignCellCycle(data=obj@data, calc=obj@dataCC, ...)
	    l[[length(l)+1]] <- data.frame(Sys.time(), "Calculate fractions ...")
	    obj@dataCC <- calcFract(data=obj@data, cuts=obj@dataCC, ...)
	    l[[length(l)+1]] <- data.frame(Sys.time(), "Finished CellCycle ...")
	    obj@log <- l
	}
    } 
    if (is.null(fun) || fun == "nc") {
	#TODO: Check if vars are there!
        l <- obj@log
        l[[length(l)+1]] <- data.frame(Sys.time(), "Started calculating Nucleus/Cell ...")
        l[[length(l)+1]] <- data.frame(Sys.time(), paste("Parameters: ", ..., collapse="|"))
        obj@dataNC <- ncAreaFilter(obj@data, ...)
        l[[length(l)+1]] <- data.frame(Sys.time(), paste("Impute cutoffs ... "))
        obj@dataNC <- imputeMC(obj@dataNC, ...)
        l[[length(l)+1]] <- data.frame(Sys.time(), "Assign MC results ...")
        obj@data <- assignNcFilterRes(obj@data, obj@dataNC, ...)
        l[[length(l)+1]] <- data.frame(Sys.time(), "Finished NC ...")
        obj@log <- l
    }
    if (is.null(fun) || fun == "fb") {
	if (!varFB %in% colnames(obj@data)) {
	    warning(paste("FB Variable (", varFB, ") not found!"))
	} else {
	    l <- obj@log
	    l[[length(l)+1]] <- data.frame(Sys.time(), "Started calculating Fibroblasten ...")
	    l[[length(l)+1]] <- data.frame(Sys.time(), paste("Parameters: ", ..., collapse="|"))
	    ## TODO: change manual var parameter?
	    obj@dataFB <- cellCycleFractIntegrDNAInt(obj@data, var=varFB, xMinGlobMax=xMinGlobMaxFB, xMaxGlobMax=xMaxGlobMaxFB, cutUpper=1, ...)
	    l[[length(l)+1]] <- data.frame(Sys.time(), paste("Impute cutoffs ... "))
	    obj@dataFB <- imputeMC(obj@dataFB, ...)
	    l[[length(l)+1]] <- data.frame(Sys.time(), "Assign Fibroblast status ...")
	    obj@data <- assignFibroblast(data=obj@data, calc=obj@dataFB, ...)
	    l[[length(l)+1]] <- data.frame(Sys.time(), "Finished Fibroblast identification ...")
	    obj@log <- l
	}
    }
    return (obj)
}

#' @title Filter cell measurments
#' @description Filters cells 
#' @export
filter <- function(obj, fun=NULL, cellCycle=c("G1", "G2/S")) {
    check(obj)
    if ((is.null(fun) || fun == "cc") && !is.null(obj@dataCC) && "CellCycle" %in% colnames(obj@data) ) {
        print("Filter Cell Cycle ...")
        l <- obj@log
        l[[length(l)+1]] <- data.frame(Sys.time(), paste("N measurments before CC filtering: ", length(obj@data[,1])))
        l[[length(l)+1]] <- data.frame(Sys.time(), paste("Retain: ", paste(cellCycle, collapse=" ")))
        obj@data <- obj@data[which(obj@data$CellCycle %in% cellCycle),,drop=F]
        l[[length(l)+1]] <- data.frame(Sys.time(), paste("N measurments after CC filtering: ", length(obj@data[,1])))
        obj@log <- l
    }
    if ( (is.null(fun) || fun == "nc") && !is.null(obj@dataNC) && "ncArea" %in% colnames(obj@data)) {
        print("Filter N/C Ratio...")
        l <- obj@log
        l[[length(l)+1]] <- data.frame(Sys.time(), paste("N measurments before NC filtering: ", length(obj@data[,1])))
        obj@data <- obj@data[which(obj@data$ncArea == 1),,drop=F]
        l[[length(l)+1]] <- data.frame(Sys.time(), paste("N measurments after NC filtering: ", length(obj@data[,1])))
        obj@log <- l
    }
    if ((is.null(fun) || fun == "fb") && !is.null(obj@dataFB) && "FIBROBLAST" %in% colnames(obj@data)) {
        print("Filter Fibroblasts ...")
        l <- obj@log
        l[[length(l)+1]] <- data.frame(Sys.time(), paste("N measurements before Fibroblast filtering: ", length(obj@data[,1])))
        obj@data <- obj@data[which(obj@data$FIBROBLAST == 0),,drop=F]
        l[[length(l)+1]] <- data.frame(Sys.time(), paste("N measurments after Fibroblast filtering: ", length(obj@data[,1])))
        obj@log <- l
    }
    return(obj)
}

#' @title Analyze Cell-Cycle results
#' @export
analyzeCC <- function(obj, ...) {
    if (is.null(obj@dataCC)) {
        stop("No cell cycle cutoffs calculated yet. Run calc()")
    }
    analysisCCFractions(obj@dataCC, ...)
}

#' @title Aggregate data
#' @export
agg <- function(obj, ...) {
    check(obj)
    ## FIXME
    if (any(grepl(":", obj@data$TREATMENT))) {
        stop(": not allowed in TREATMENT!")
    }
    l <- obj@log
    l[[length(l)+1]] <- data.frame(Sys.time(), paste("Start aggregating data per well - nRow: ", length(obj@data[,1])))
    l[[length(l)+1]] <- data.frame(Sys.time(), paste("Parameters: ", ..., collapse="|"))
    obj@dataAgg <- aggregateByWell(obj@data, ...)
    l[[length(l)+1]] <- data.frame(Sys.time(), paste("Finished aggregating data - nRow: ", length(obj@dataAgg[,1])))
    obj@log <- l
    return(obj)
}

#' @title Get Log
#' @export 
getLog <- function(obj) {
    df <- lapply(obj@log, function(x) data.frame(TIME=x[[1]][1], MSG=x[[2]][1]))
    return(do.call(rbind, df))
}

#' @title Remove raw data
#' @export
removeRawData <- function(obj) {
    if (is.null(obj@dataAgg)) {
        stop("No aggregated data found! Run agg()! Exiting!")
    }
    obj@data <- NULL
    return(obj)
}

#' @title QC Plots
#' @export 
qcPlots <- function(obj, folder=NULL, ...) {
    ## FIXME: no filename given -> ERROR
    if (is.null(obj@data)) {
        print("No data found! Re-loading original data!")
        obj <- getData(obj)
    }
    ret <- qc(obj@data, uid=obj@expID, folder=folder, ...)
    print(ret[[1]])

    l <- obj@log 
    l[[length(l)+1]] <- data.frame(Sys.time(), paste("Created QC Plots: ", ret[[1]]))
    obj@log <- l
    return(obj)
}

check <- function(obj) {
    if (is.null(obj@data)) {
        stop("No data available! Load experiment data with getData()")
    }
}

#' @title Create Ratio/Diff Vars
#' @export 
addVars <- function(obj, addVars) {
    check(obj)
    l <- obj@log
    l[[length(l)+1]] <- data.frame(Sys.time(), paste("Add Ratio/Diff Vars. Cols: ", length(obj@data[1,])))
    l[[length(l)+1]] <- data.frame(Sys.time(), paste(addVars, collapse="|"))
    obj@data <- createAdditionalVars(obj@data, vars=addVars)
    l[[length(l)+1]] <- data.frame(Sys.time(), paste("Finished Ratio/Diff Vars. Cols:", length(obj@data[1,])))
    obj@log <-l 
    return(obj)
}

#' @title Z-Transform data
#' @export 
zTrans <- function(obj, versuche=NULL) {
    if (is.null(obj@dataAgg)) {
        stop("Only aggregated data are transformed. Not found. Run agg())")
    }
    if (is.null(versuche)) {
        versuche <- unique(obj@dataAgg$VERSUCH)
    } 
    l <- obj@log
    l[[length(l)+1]] <- data.frame(Sys.time(), paste("Z transform data. Versuche: ", paste(versuche, collapse="|")))
    obj@dataZ <- do.call(rbind, zTransform(obj@dataAgg, versuch=versuche))
    l[[length(l)+1]] <- data.frame(Sys.time(), paste("Z transform data finished."))
    obj@log <- l
    return(obj)
}


#' @title Check and filter treatments
#' @param filter remove treated categories with small number of cells,
#' defaults to NULL (no filtering)
#' @export 
checkFilter <- function(obj, filter=NULL) {
    tbl <- table(TREATMENT=obj@data$TREATMENT, VERSUCH=obj@data$VERSUCH, PLATTE=obj@data$PLATTE)

    #### remove NA treatments
    wNA <- which(is.na(obj@data$TREATMENT ))
    if (length(wNA) > 0) {
	warning(paste("Found ", length(wNA), " NA treatments! Removing!", sep=""))
	obj@data <- obj@data[-wNA,,drop=F]
    }

    df <- data.frame(tbl)
    print(df)
    if (!is.null(filter)) {
	df <- df[which(df$Freq < filter),,drop=F]
	#TODO: check dim
	w <- c()
	for (i in 1:length(df[,1])) {
	    ww <- which(obj@data$VERSUCH == df$VERSUCH[i] &
			obj@data$TREATMENT == df$TREATMENT[i] &
			obj@data$PLATTE == df$PLATTE[i])
	    w<-c(w, ww)
	}
	if (length(w) > 0) {
	    obj@data <- obj@data[-w,,drop=F]
	}
    }
    return(obj)
}
