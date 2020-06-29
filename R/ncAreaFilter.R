#' @title Nucl/Cell Ratio Filter
#' 
#' @description desca
#' 
#' @param data data
#' 
#' @import parallel
#' @import foreach
#' @import doParallel
#' 
#' @export
ncAreaFilter <- function(data, varC="AreaShape_Area.cell",
                         varN="AreaShape_Area.nucl",
                         lowerCL=0.1, upperCL=0.9,
                         nBoot=100, xCut=15, bw=0.1,
                         xMinGlobMax=1.5, ...) {
  #retain only necessary columns
    data <- data[,which(colnames(data) %in% c(varC, varN, "TREATMENT", "VERSUCH", "PLATTE"))]

  #parallelization
  no_cores <- parallel::detectCores() - 1
  no_cores <- ifelse(no_cores == 0, 1, no_cores)
  no_cores <- ifelse(is.na(no_cores), 1, no_cores)
  doParallel::registerDoParallel(no_cores)
  
  ##Calculations for each plate separately
  plateLv <- unique(paste(data$VERSUCH, data$PLATTE))
  ret <- list()
  
  for (plate in plateLv) {
    print(paste("Processing", plate))
    ## selektiere plate for each experiment
    sub <- data[which(paste(data$VERSUCH, data$PLATTE) == plate),]
    
    ## calculate for each treatment - level
    out <- list()
    lvs <- unique(sub$TREATMENT)
    for(trLv in lvs) {      
        cat(paste("\r                ", trLv, "  ", round(which(trLv == lvs)/length(lvs)*100, 2), "%", sep=""))
      subTr0 <- sub[which(sub$TREATMENT == trLv),]
      
      ## permutation
      coll <- foreach (i=1:nBoot) %dopar% {
        if (i %% 100 == 0) { cat(paste("\r   ", round(i/nBoot*100,2), "%        ",sep="")) }
        sel <- sample(1:length(subTr0[,1]), length(subTr0[,1]), replace=T)
        subTr <- subTr0[sel,]
        
        rat <- subTr[,varC]/subTr[,varN]
        d <- density(rat[which(rat < xCut)], bw=bw)
        #plot(d)
        ## get extremal points (max, min)
        d1 <- diff(d$y)
        d2 <- diff(d1)
        v1 <- d1[-1]
        v2 <- d1[-length(d1)]
        maxCandPos <- which(v1 > 0 & v2 < 0 | v1 < 0 & v2 > 0) 
        
        ##alle maxima / minima bei auto-bandwidth
        maximaPos <- maxCandPos[which(d2[maxCandPos] < 0)]+1
        minimaPos <- maxCandPos[which(d2[maxCandPos] > 0)]+1
        
        ##identifizierte globales maximum mit xVal > xMinGlobMax
        maxPos <- which(d$y == max(d$y[which(d$x > xMinGlobMax)], na.rm=T))
        
        ## identifzierte minimum links davon, gt 1
        minLeftPos <- which(d$y == min(d$y[which(d$x > 1 & d$x < d$x[maxPos])]))
        
        ## collect
	v1 <- ifelse(length(maxPos) == 0, NA, d$x[maxPos])
	v2 <- ifelse(length(minLeftPos) == 0, NA, d$x[minLeftPos])
        #data.frame(d$x[maxPos], d$x[minLeftPos])
        data.frame(v1, v2)
      }
      coll <- do.call(rbind, coll)
      colnames(coll) <- c("maxPos", "minLeftPos")
      
      ##cehck for NAs -> extremal point could not be found!
      nas <- apply(coll, 2, function(x) any(is.na(x)))
      if (any(is.na(nas))) {
        warning("NAs detected!")
      }
      
      ##confidence intervals
      lower <- apply(coll, 2, function(x) quantile(x, probs=lowerCL, na.rm=T)[[1]])
      upper <- apply(coll, 2, function(x) quantile(x, probs=upperCL, na.rm=T)[[1]])
      estim <- apply(coll, 2, function(x) quantile(x, probs=0.5, na.rm=T)[[1]])
      
      ##Total cells
      nTotal <- length(subTr0[,1])
      
      ## FIXME
      if (F) {
	  ##Filtered Cells
	  nFiltered <- c(length(which(subTr0[,varC]/subTr0[,varN] < c(lower[2]))),
			 length(which(subTr0[,varC]/subTr0[,varN] < c(upper[2]))),
			 length(which(subTr0[,varC]/subTr0[,varN] < c(estim[2]))))
	  names(nFiltered) <- paste("nFiltered", c("lower", "upper", "median"), sep="_")
      }
      
      
      out[[length(out)+1]] <- list(treatment=trLv,
                                   nTotal=nTotal,
                                   #nFiltered=nFiltered, 
                                   #fracFiltered=nFiltered/nTotal, 
                                   lower=lower,
                                   upper=upper,
                                   estim=estim,
                                   widthCI=upper-lower)
      
    }
    pv <- do.call(rbind, strsplit(as.character(plate), " "))
    ret[[length(ret)+1]] <- list(data=out, versuch=pv[1,1], platte=pv[1,2])
  }
  
  return(ret)
  
  doParallel::stopImplicitCluster()
}
