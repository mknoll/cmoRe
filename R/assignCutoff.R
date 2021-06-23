#' @title Assigns Cutoff result to each cell
#' 
#' @description A
#' 
#' @param data data
#' @param calc Data calculated with cellCycleFractIntegrDNAInt()
#' @param type Method on how to assign cell cycle status
#' to each cell. CIexcl assigns only cells with values 
#' outside the confidence interval for the separating value 
#' determined in cellCycleFractIntegrDNAInt() by bootstrapping.
#' CImedian uses the median as cutoff.
#' 
#' @import parallel
#' @import foreach
#' @import doParallel
#' 
#' @export
assignCutoff <- function(data, calc, type="CIexcl",
                            var="Intensity_IntegratedIntensity_DNA.nucl", 
                            nCores=NULL, outvar=NULL, ...) {
  #parallelization
    if (is.null(nCores)) {
        nCores <- parallel::detectCores() - 1
        nCores <- ifelse(nCores == 0, 1, nCores)
        nCores <- ifelse(is.na(nCores), 1, nCores)
    } 
  doParallel::registerDoParallel(nCores)

  calcDF <- do.call(rbind, calc)
  lvPV <- paste(calcDF[,"versuch"], calcDF[,"platte"])
  
  if (is.null(outvar)) {
      outvar <- "Cutoff"
  }
  outvarW <- paste(outvar, "Width", sep="")
  add <- data.frame(rep(NA, length(data[,1])),
                    rep(NA, length(data[,1])))
  colnames(add) <- c(outvar, outvarW)
  data <- cbind(data, add)
  #data$Cutoff <- NA
  #data$CutoffWidth <- NA

  out <- foreach(lv=lvPV) %dopar% {
      dataSub <- data[which(paste(data$VERSUCH, data$PLATTE) == lv),]
      ##by experiment / plate
      #for (lv in lvPV) {
      print(paste("Processing", lv))
      cc <- do.call(rbind, calcDF[which(paste(calcDF[,"versuch"], calcDF[,"platte"]) == lv),1][[1]])
      ## by treatment
      for (tr in unique(unlist(cc[,"treatment"]))) {
          #print(tr)
          val <- cc[which(cc[,"treatment"] == tr),]
          sel <- which(dataSub$TREATMENT == tr)

          ##Type of assignment
          if (type=="CIexcl") {
              ##left
              #dataSub$Cutoff[sel] <- ifelse(dataSub[sel,var] < val$lower[["g1MinRightPos"]], "LOW", dataSub$Cutoff[sel]) 
              dataSub[sel,outvar] <- ifelse(dataSub[sel,var] < val$lower[["g1MinRightPos"]], "LOW", dataSub[sel,outvar]) 
              ##rigth
              #dataSub$Cutoff[sel] <- ifelse(dataSub[sel,var] > val$upper[["g1MinRightPos"]], "HIGH", dataSub$Cutoff[sel]) 
              dataSub[sel,outvar] <- ifelse(dataSub[sel,var] > val$upper[["g1MinRightPos"]], "HIGH", dataSub[sel, outvar]) 
              ##width  - check for max width
              #dataSub$CellCycle[sel] <- ifelse(dataSub[sel,var] > val$upper[["g1MinLeftPos"]] & 
              #                                 dataSub[sl,var] < val$lower[["g1MinRightPos"]], "G1", dataSub$CellCycle[sel]) 
              #dataSub$CutoffWidth[sel] <- val$upper[["g1MinLeftPos"]] - val$lower[["g1MinRightPos"]]
              dataSub[sel,outvarW] <- val$upper[["g1MinLeftPos"]] - val$lower[["g1MinRightPos"]]
          } else if (type == "CImedian") {
              ##left
              #dataSub$Cutoff[sel] <- ifelse(dataSub[sel,var] < val$estim[["g1MinRightPos"]], "LOW", dataSub$Cutoff[sel]) 
              dataSub[sel,outvar] <- ifelse(dataSub[sel,var] < val$estim[["g1MinRightPos"]], "LOW", dataSub[sel,outvar]) 
              ##right
              #dataSub$Cutoff[sel] <- ifelse(dataSub[sel,var] > val$estim[["g1MinRightPos"]], "G2/S", dataSub$Cutoff[sel]) 
              dataSub[sel,outvar] <- ifelse(dataSub[sel,var] > val$estim[["g1MinRightPos"]], "HIGH", dataSub[sel,outvar]) 
              ##g1 cells
              #dataSub$CellCycle[sel] <- ifelse(is.na(dataSub$CellCycle[sel]), "G1", dataSub$CellCycle[sel])
          } else {
              stop("Unknown type! Can be CIexcl or CImedian!")
          }
      }
      dataSub
  }
  data <- do.call(rbind, out)

  ## Additional variables
  df <- data.frame(high=ifelse(data[,outvar] == "HIGH", 1, 0),
                   low=ifelse(data[,outvar] == "LOW", 1, 0))
  cn <- paste(outvar, "_BINARY_", c("HIGH", "LOW"), sep="")
  colnames(df) <- cn
  data <- cbind(data, df)

  doParallel::stopImplicitCluster()

  return (data)
}


#' @title assign cutoffs imple
#' @export
assignCutoffSimple <- function(data, calc, var) {
    ## out vector
    group <- rep(NA, length(data[,1]))
    calcDF <- do.call(rbind, calc)
    dat0 <- data[,var]
    v <- unlist(calcDF[,"versuch"])      
    p <- unlist(calcDF[,"platte"])                          
    for (i in 1:length(calcDF[,1])) {           
	cat(".")    
	w <- which(data$VERSUCH == v[i] & data$PLATTE == p[i])    
	treat <- data$TREATMENT[w]    
	for (tr in unique(unlist(lapply(calc[[i]]$data, function(x) x$treatment)))) {    
	    ww <-which(treat == tr)        
	    val <- lapply(calc[[i]]$data, function(x)  if (x$treatment == tr) {  x$estim })        
	    val <- Filter(length,val)[[1]]                                                                 
	    #group[w][ww][which(dat0[w][ww] < val[["g1MinLeftPos"]])] <- "DEAD"         
	    group[w][ww][which(dat0[w][ww] < val[["g1MinRightPos"]])] <- "LOW"        
	    group[w][ww][which(dat0[w][ww] > val[["g1MinRightPos"]])] <- "HIGH"        
	}        
    }    
    group
}
