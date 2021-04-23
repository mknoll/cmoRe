#' @title Read CP file single column
#' 
#' @param file CP file
#' @param var CP column name
#' 
#' @import data.table
#'
#' @useDynLib cmoRe
#'
#' @export
readSingleCol <- function(file, var, nrow=NULL, type="numeric") {
    if (!file.exists(file) || dir.exists(file)) {
	warning("Not a regular file!") ### FIXME
	return(NA)
    }

    #get rows in file
    if (is.null(nrow)) {
	nrow <- getNRows(file)
	if (is.na(row)) {
	    warning("Not a regular file!") ### FIXME
	    return(NA)
	}
    }

    #find col pos
    tmp <- data.frame(fread(file, nrow=10))
    ncol <- as.integer(which(colnames(tmp) == var))
    print(paste0("ncol: ", ncol))
    if (length(ncol) == 0) {
	stop(paste0("Could not find ", var))
    }


    if (type == "numeric") {
	data <- numeric(nrow)
	data <- .C("readCol", file, nrow, ncol, data)[[4]][-c(1:2)]
    } else if (type == "character") {
	data <- character(nrow)
	data <- .C("readColChar", file, nrow, ncol, data)[[4]][-c(1:2)]
    }
    return(data)
}

#' @title Get number of rows from file
#' @title Number of rows (including header)
#' @param file Filename
#' @export
getNRows <- function(file) {
    if (!file.exists(file) || dir.exists(file)) {
	warning("Not a regular file!") ## FIXME
	return(NA)
    }

    #get rows in file
    nrow <- as.integer(1)
    nrow <- .C("getNRows", file, nrow)[[2]]
    return(as.integer(nrow))
}
