#' @title median aggregation
#'
#' @useDynLib cmoRe
#'
#' @export
medianAgg <-function(obj) {
    ### cols to transform
    ww <- which(!is.na(as.numeric(obj@data[1,])) & !colnames(obj@data) %in% c("PLATTE", "VERSUCH"))

    id <-paste(droplevels(factor(obj@data$VERSUCH)),droplevels(factor(obj@data$PLATTE)),droplevels(factor(obj@data$WELL)))
    pos <- list()    
    for (i in unique(id)) {    
	w <- which(id == i)    
	pos[[length(pos)+1]] <- w 
    }    

    coll <- list()
    tr <- list()
    for (w in ww) {
	#print(colnames(obj@data)[w])
	tmp <- list()
	for (q in 1:length(pos)) {
	    res <-integer(1)    
	    val <- obj@data[pos[[q]],w]
	    val <- val[which(!is.na(val))]
	    n<-as.integer(length(val))    
	    if (length(val) > 2) {
		y<-as.integer(val*10000)    
		ret <- .C("median", n, y, res)    
		tmp[[length(tmp)+1]] <- ret[[3]]/10000
	    } else {
		tmp[[length(tmp)+1]] <- NA
	    }
	}
	coll[[length(coll)+1]] <- unlist(tmp)
    }
    df <- do.call(cbind, coll)
    colnames(df) <- colnames(obj@data)[ww]
    print(tr)
    add <- data.frame(do.call(rbind, strsplit(unique(id), " ")))
    colnames(add) <- c("VERSUCH", "PLATTE", "WELL")
    add$TREATMENT <- obj@data$TREATMENT[match(paste(add$VERSUCH,add$PLATTE, add$WELL), id)]
    data.frame(df, add)
}
