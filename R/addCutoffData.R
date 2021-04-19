#' @title Adds separately calculated cutoff data
#'
#' @export
addCutoffs <- function(obj, meta) {
    if (length(meta$dataCC) > 0) {
	obj@dataCC <- meta$dataCC
	obj@dataCC <- imputeMC(obj@dataCC)
	obj@data <- assignCellCycle(data=obj@data, calc=obj@dataCC)
	obj@dataCC <- calcFract(data=obj@data, cuts=obj@dataCC)
    }
    if (length(meta$dataNC) > 0) {
	obj@dataNC <- meta$dataNC
	obj@dataNC <- imputeMC(obj@dataNC)
	obj@data <- assignNcFilterRes(obj@data, obj@dataNC)
    }
    if (length(meta$dataFB) > 0) {
	obj@dataFB <- imputeMC(obj@dataFB)
	obj@data <- assignFibroblast(data=obj@data, calc=obj@dataFB)
    }
    return(obj)
}
