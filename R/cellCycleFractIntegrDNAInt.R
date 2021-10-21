#' @title Derive cell cycle fractions from nuclear 
#' DNA intensity measurements
#' 
#' @description Calculates cell-cycle fractions (G1, G2) 
#' of cells based on a nuclear DNA intensity variable; 
#' defaults to "Intensity_IntegratedIntensity_DNA".
#' 
#' @param data data.frame containing at least the columns 
#'  TREATMENT, VERSUCH, PLATTE and the variable used for 
#'  DNA cell cycle inference (see parameter var)
#' @param var Variable in which the nuclear integrated dna 
#'  intensity is stored. Default is 
#'  "Intensity_IntegratedIntensity_DNA.nucl"
#' @param nBoot number of runs to determine minima/maxima
#'  to derive cutoffs
#' @param lowerCL lower confidence level, default: 0.1
#' @param upperCL upper confidence level, default: 0.9
#' @param cutUpper use only data below a certain cutoff, 
#' default: NULL -> not used
#' @param nMin minimum required number of cells
#' 
#' @import parallel
#' @import foreach
#' @import doParallel
#' 
#' @export
#' 
#' @return  data.frame containing the computed values
cellCycleFractIntegrDNAInt <- function(data, var="Intensity_IntegratedIntensity_DNA.nucl",
                                       nBoot=100, lowerCL=0.1, upperCL=0.9, log=F, 
				       cutUpper=NULL,
				       xMinGlobMax=0, xMaxGlobMax=NULL, 
				       no_cores=NULL, nMin=100, ..) {
  #retain only necessary columns
    data <- data[,which(colnames(data) %in% c(var, "TREATMENT", "VERSUCH", "PLATTE"))]
    print("-----")
  print(paste("Processing:", var))

  #parallelization
  if (is.null(no_cores)) {
    no_cores <- parallel::detectCores() - 1
  }
  no_cores <- ifelse(no_cores == 0, 1, no_cores)
  no_cores <- ifelse(is.na(no_cores), 1, no_cores)
  doParallel::registerDoParallel(no_cores)
  
  ##Calculations for each plate separately
  plateLv <- unique(paste(data$VERSUCH, data$PLATTE))
  ret <- list()
  
  for (plate in plateLv) {
    print(paste("Processing", plate))
    ## selektiere plate for each experiment
    sub <- data[which(paste(data$VERSUCH, data$PLATTE) == plate),,drop=F]
    ## apply cutoff for upper data
    if (!is.null(cutUpper)) {
	sub <- sub[which(sub[,var] <= cutUpper),,drop=F]
    }
    
    ## calculate for each treatment - level
    out <- list()
    lvs <- unique(sub$TREATMENT)
    lvs <- lvs[which(!is.na(lvs))]
    for(trLv in lvs) {      
	cat(paste("\r             ", trLv, "  ", round(which(trLv == lvs)/length(lvs)*100, 2), "%", sep=""))
	subTr0 <- sub[which(sub$TREATMENT == trLv),var,drop=F]
	if (length(subTr0[,1]) < nMin) {
	    v1 <- NA
	    v2 <- NA
	    v3 <- NA
	} else {

	    ##### TODO: check if sufficient data is available!
	    if (length(subTr0[,1]) <= 1) {
		### FIXME -> FAIL
		warning("Not enough data!")
	    }
	    if (log) {
		subTr0 <- data.frame(apply(subTr0, 2, log))
	    }

	    ## permutation
	    coll <- foreach (i=1:nBoot) %dopar% {
		if (i %% 100 == 0) { cat(paste("\r   ", round(i/nBoot*100,2), "%        ",sep="")) }
		sel <- sample(1:length(subTr0[,1]), length(subTr0[,1]), replace=T)
		subTr <- subTr0[sel,,drop=F]

		d <- density(subTr[,var], na.rm=T)
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

		##identifizierte globales maximum -> i.d.R. G1 peak
		if (is.null(xMaxGlobMax)) {
		    g1MaxPos <- which(d$y == max(d$y[which(d$x > xMinGlobMax)]))
		} else {
		    g1MaxPos <- which(d$y == max(d$y[which(d$x > xMinGlobMax & d$x < xMaxGlobMax)]))
		}

		## identifiziere 1. min links/rechts davon
		g1MinLeftPos <- rev(minimaPos[which(g1MaxPos-minimaPos > 0)])[1]
		g1MinRightPos <- minimaPos[which(g1MaxPos-minimaPos < 0)][1]

		##identifiziere lokales maximum links von glob -> tote zellen
		deadMaxPos <- which(d$y == max(d$y[d$x < d$x[g1MinLeftPos]] ))

		##identifizierte lokales maximum rechts von glob -> G2 peak
		g2MaxPos <- which(d$y == max(d$y[d$x > d$x[g1MinRightPos]]))

		## collect
		#data.frame(d$x[g1MaxPos], d$x[g1MinLeftPos], d$x[g1MinRightPos])
		v1 <- ifelse(length(g1MaxPos) == 0, NA, d$x[g1MaxPos])
		v2 <- ifelse(length(g1MinLeftPos) == 0, NA, d$x[g1MinLeftPos])
		v3 <- ifelse(length(g1MinRightPos) == 0, NA, d$x[g1MinRightPos])
	    }
	  data.frame(v1, v2, v3)
      }


      coll <- do.call(rbind, coll)
      if (log) {
	  coll <- data.frame(apply(coll, 2, exp))
      }

      colnames(coll) <- c("g1MaxPos", "g1MinLeftPos", "g1MinRightPos")

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

      #### TODO
      if (F) {
	  ### FIXME -> update nach imputation!
	  ##Dead cells
	  nDead <- c(length(which(subTr0[,var] < c(lower[2]))),
		     length(which(subTr0[,var] < c(upper[2]))),
		     length(which(subTr0[,var] < c(estim[2]))))
	  names(nDead) <- paste("nDead", c("lower", "upper", "median"), sep="_")

	  ##Cells in G2
	  nG2 <- c(length(which(subTr0[,var] > c(lower[3]))),
		   length(which(subTr0[,var] > c(upper[3]))),
		   length(which(subTr0[,var] > c(estim[3]))))
	  names(nG2) <- paste("nG2", c("lower", "upper", "median"), sep="_")


	  ##Cells in G1
	  nG1 <- nTotal-nDead-nG2
	  names(nG1) <- paste("nG1", c("lower", "upper", "median"), sep="_")
      }
      
        
      out[[length(out)+1]] <- list(treatment=trLv,
            nTotal=nTotal,
            #nDead=nDead, 
            #nG2=nG2, 
            #nG1=nG1, 
            #fracDead=nDead/nTotal, 
            #fracG2=nG2/nTotal, 
            #fracG1=nG1/nTotal,
            lower=lower,
            upper=upper,
            estim=estim)
            #widthCI=upper-lower)
          
      }
    pv <- do.call(rbind, strsplit(as.character(plate), " "))
    ret[[length(ret)+1]] <- list(data=out, versuch=pv[1,1], platte=pv[1,2])
  }
  
  doParallel::stopImplicitCluster()
      
  return(ret)
  
}
