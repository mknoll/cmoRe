#' @title Adds separately calculated cutoff data
#'
#' @export
addCutoffs <- function(obj, meta, type=c("nc", "fb", "cc")) {
    ### check for missing data
    if (length(meta$dataCC) > 0 && "cc" %in% type) {
	len <- lapply(meta$dataCC, function(x) length(x))
	if (any(unlist(len) == 1)) { 
	    stop("Missing data! Please rerun calcCutoffs with already retrieved data: 
		 calcCutoff(obj, cutoffs)!")
	}
    }
    if (length(meta$dataNC) > 0 && "nc" %in% type) {
	len <- lapply(meta$dataNC, function(x) length(x))
	if (any(unlist(len) == 1)) { 
	    stop("Missing data! Please rerun calcCutoffs with already retrieved data: 
		 calcCutoff(obj, cutoffs)!")
	}
    }
    if (length(meta$dataFB) > 0 && "fb" %in% type) {
	len <- lapply(meta$dataFB, function(x) length(x))
	if (any(unlist(len) == 1)) { 
	    stop("Missing data! Please rerun calcCutoffs with already retrieved data: 
		 calcCutoff(obj, cutoffs)!")
	}
    }


    if (length(meta$dataCC) > 0 && "cc" %in% type) {
	obj@dataCC <- lapply(meta$dataCC, function(x) x[[1]])
	obj@dataCC <- imputeMC(obj@dataCC)
	#obj@data <- assignCellCycle(data=obj@data, calc=obj@dataCC)
	cc <- assignCellCycleSimple(obj, obj@dataCC)
	#obj@data <- data.frame(obj@data, CellCycle=cc)
	obj@data <- data.frame(obj@data, 
			       CellCycle=cc, 
			       CellCycle_DEAD=ifelse(cc == "DEAD", 1, 0), 
			       CellCycle_G1=ifelse(cc == "G1", 1, 0), 
			       CellCycle_G2S=ifelse(cc == "G2/S", 1,0))
	#obj@dataCC <- calcFract(data=obj@data, cuts=obj@dataCC)
    }
    if (length(meta$dataNC) > 0 && "nc" %in% type) {
	obj@dataNC <- lapply(meta$dataNC, function(x) x[[1]])
	obj@dataNC <- imputeMC(obj@dataNC)
	nc <- assignNcFilterResSimple(obj, obj@dataNC)
	obj@data <- data.frame(obj@data, ncArea=nc)
	#obj@data <- assignNcFilterRes(obj@data, obj@dataNC)
    }
    if (length(meta$dataFB) > 0 && "fb" %in% type) {
	obj@dataFB <- lapply(meta$dataFB, function(x) x[[1]])
	obj@dataFB <- imputeMC(obj@dataFB)
	fb <- assignFibroSimple(obj, obj@dataFB)
	obj@data <- data.frame(obj@data, FIBROBLAST=fb)
	#obj@data <- assignFibroblast(data=obj@data, calc=obj@dataFB)
    }
    return(obj)
}
