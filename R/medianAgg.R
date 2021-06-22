#' @title median aggregation
#'
#' @useDynLib cmoRe
#'
#' @export
medianAgg <-function(obj) {
    ### cols to transform
    ww <- which(!is.na(as.numeric(obj@data[1,])) & !colnames(obj@data) %in% c("PLATTE", "VERSUCH") & !grepl("CellCycle", colnames(obj@data)))

    id <-paste(droplevels(factor(obj@data$VERSUCH)),droplevels(factor(obj@data$PLATTE)),droplevels(factor(obj@data$WELL)))
    pos <- list()    
    for (i in unique(id)) {    
	w <- which(id == i)    
	pos[[length(pos)+1]] <- w 
    }    

    coll <- list()
    tr <- list()
    i <- 1
    for (w in ww) {
	if (i %% 10 == 0) {
	    cat(paste("\r     "), round(i/length(ww)*100, 2), "%")
	}
	i<-i+1
	#print(colnames(obj@data)[w])
	tmp <- list()
	for (q in 1:length(pos)) {
	    res <-integer(1)    
	    val <- obj@data[pos[[q]],w]
	    fact <- 1
	    if (length(val) > 2 && max(val[which(!is.na(val) & !is.infinite(val))]) < 1000) {
		fact <- 100000
	    }
	    val <- val*fact ###ACHTUNG!!!! Kann zu infinite fuehren!
	    val <- val[which(!is.na(val) & !is.infinite(val))]
	    n<-as.integer(length(val))    
	    if (length(val) > 2) {
		y<-as.integer(val)    
		## remove NA in y
		if (any(is.na(y))) {
		    rmY <- which(is.na(y))
		    y <- y[-rmY]
		    n <- as.integer(length(y))
		}
		if (length(y) > 2) {
		    ret <- .C("median", n, y, res)    
		    tmp[[length(tmp)+1]] <- ret[[3]]/fact
		} else {
		    tmp[[length(tmp)+1]] <- NA
		}
	    } else {
		tmp[[length(tmp)+1]] <- NA
	    }
	}
	coll[[length(coll)+1]] <- unlist(tmp)
    }
    df <- do.call(cbind, coll)
    colnames(df) <- colnames(obj@data)[ww]
    add <- data.frame(do.call(rbind, strsplit(unique(id), " ")))
    colnames(add) <- c("VERSUCH", "PLATTE", "WELL")
    add$TREATMENT <- obj@data$TREATMENT[match(paste(add$VERSUCH,add$PLATTE, add$WELL), id)]
    data.frame(df, add)
}
