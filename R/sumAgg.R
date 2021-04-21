#' @title sum aggregation to calculate fraction
#'
#' @useDynLib cmoRe
#'
#' @export
sumAgg <-function(obj,var="CellCycle_") {
    ### cols to transform
    ww <- which(grepl(var, colnames(obj@data)))

    id <-paste(droplevels(factor(obj@data$VERSUCH)),droplevels(factor(obj@data$PLATTE)),droplevels(factor(obj@data$WELL)))
    pos <- list()    
    nc <- list() #number of cells
    for (i in unique(id)) {    
	w <- which(id == i)    
	pos[[length(pos)+1]] <- w 
	nc[[length(nc)+1]] <- length(w)
    }    

    coll <- list()
    i <- 1
    for (w in ww) {
	if (i %% 10 == 0) {
	    cat(paste("\r     "), round(i/length(ww)*100, 2), "%")
	}
	i<-i+1
	tmp <- list()
	for (q in 1:length(pos)) {
	    res <-integer(1)    
	    val <- obj@data[pos[[q]],w]
	    val <- val[which(!is.na(val) & !is.infinite(val))]
	    if (length(val) > 2) {
		tmp[[length(tmp)+1]] <- sum(val)/length(val)
	    } else {
		tmp[[length(tmp)+1]] <- NA
	    }
	}
	coll[[length(coll)+1]] <- unlist(tmp)
    }
    df <- data.frame(do.call(cbind, coll))
    colnames(df) <- colnames(obj@data)[ww]
    colnames(df) <-paste0(colnames(df), "_FRAC")
    df$NCELLS <- unlist(nc)
    add <- data.frame(do.call(rbind, strsplit(unique(id), " ")))
    colnames(add) <- c("VERSUCH", "PLATTE", "WELL")
    add$TREATMENT <- obj@data$TREATMENT[match(paste(add$VERSUCH,add$PLATTE, add$WELL), id)]
    data.frame(df, add)
}
