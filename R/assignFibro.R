#' @title Assigns Fibroblast status to each cell
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
assignFibroblast <- function(data, calc, type="CImedian",
                            var="Intensity_MedianIntensity_DNA.nucl", 
                            nCores=NULL, ...) {
  #parallelization
    if (is.null(nCores)) {
        nCores <- parallel::detectCores() - 1
        nCores <- ifelse(nCores == 0, 1, nCores)
        nCores <- ifelse(is.na(nCores), 1, nCores)
    } 
  doParallel::registerDoParallel(nCores)

  calcDF <- do.call(rbind, calc)
  lvPV <- paste(calcDF[,"versuch"], calcDF[,"platte"])
  data$FIBROBLAST <- 0

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

          ##Type of cellcycle assignment
          if (type=="CIexcl") {
	      ## fibroblast
	      dataSub$FIBROBLAST[sel] <- ifelse(dataSub[sel,var] < val$lower[["g1MinLeftPos"]], 1, dataSub$FIBROBLAST[sel]) 
          } else if (type == "CImedian") {
	      ## fibroblast 
              dataSub$FIBROBLAST[sel] <- ifelse(dataSub[sel,var] < val$estim[["g1MinLeftPos"]], 1, dataSub$FIBROBLAST[sel]) 
          } else {
              stop("Unknown type! Can be CIexcl or CImedian!")
          }
      }
      dataSub
  }
  data <- do.call(rbind, out)
  
  doParallel::stopImplicitCluster()

  return (data)
}
