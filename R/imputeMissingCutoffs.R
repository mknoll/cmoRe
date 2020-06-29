#' @title Impute missing cutoffs
#'
#' @description Fills missing data
#'
#' @param cuts output of the cellCycleFractIntegrDNAInt function
#' @param plot plot estimates
#' 
#' @export
imputeMC <- function(cuts, plot=F, ...) {
    ## TODO: Update n and fract 
    ## Imput lower /upper 

    for (i in 1:length(cuts)) {
	print(i)
	tmp <- data.frame(do.call(rbind, lapply(cuts, function(x) lapply(x$data, function(y) y$estim))[[i]]), do.call(rbind, lapply(cuts, function(x) lapply(x$data, function(y) y$treatment))[[i]]))
	if (length(tmp[1,]) == 3) {
	    colnames(tmp)[3] <- "TREATMENT"
	    to <-2
	} else {
	    colnames(tmp)[4] <- "TREATMENT"
	    to <- 3
	}
	if (length(tmp[,1]) == 1) {
	    tmp[,1:to] <- data.frame(t(apply(tmp[,1:to], 2, function(x) as.numeric(as.character(unlist(x, recursive=T))))))
	} else {
	    tmp[,1:to] <- data.frame(apply(tmp[,1:to], 2, function(x) as.numeric(as.character(unlist(x, recursive=T)))))
	}

	## fill NA
	wNA <- which(apply(tmp, 1, function(x) any(is.na(x))))
	if (length(wNA) > 0) {
	    print(paste("Found NAs: ", length(wNA)))
	    vals <- apply(tmp[-wNA,1:to,drop=F], 2, function(x) median(x))
	    for (z in 1:to) {
		tmp[wNA,z] <- vals[z]
	    }
	}

	if (plot) {
	    matplot(tmp[,1:to], type="b")
	}

	##insert
	for (j in 1:length(cuts[[i]]$data)) {
	    #print(cuts[[i]]$data[[j]]$estim)
	    for (k in 1:to) {
		#print(tmp[j,k])
		cuts[[i]]$data[[j]]$estim[[k]] <- tmp[j,k]
	    }
	}
    }

    return(cuts)
}
