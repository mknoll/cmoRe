#' @title Assigns Cell-Cycle Status to each cell
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
assignCellCycle <- function(data, calc, type="CImedian",
                            var="Intensity_IntegratedIntensity_DNA.nucl", 
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
  data$CellCycle <- NA

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
              ##dead cells
              dataSub$CellCycle[sel] <- ifelse(dataSub[sel,var] < val$lower[["g1MinLeftPos"]], "DEAD", dataSub$CellCycle[sel]) 
              ##g2 cells
              dataSub$CellCycle[sel] <- ifelse(dataSub[sel,var] > val$upper[["g1MinRightPos"]], "G2/S", dataSub$CellCycle[sel]) 
              ##g1 cells
              dataSub$CellCycle[sel] <- ifelse(dataSub[sel,var] > val$upper[["g1MinLeftPos"]] & 
                                               dataSub[sel,var] < val$lower[["g1MinRightPos"]], "G1", dataSub$CellCycle[sel]) 
          } else if (type == "CImedian") {
              ##dead cells
              dataSub$CellCycle[sel] <- ifelse(dataSub[sel,var] < val$estim[["g1MinLeftPos"]], "DEAD", dataSub$CellCycle[sel]) 
              ##g2 cells
              dataSub$CellCycle[sel] <- ifelse(dataSub[sel,var] > val$estim[["g1MinRightPos"]], "G2/S", dataSub$CellCycle[sel]) 
              ##g1 cells
              dataSub$CellCycle[sel] <- ifelse(is.na(dataSub$CellCycle[sel]), "G1", dataSub$CellCycle[sel])
          } else {
              stop("Unknown type! Can be CIexcl or CImedian!")
          }
      }
      dataSub
  }
  data <- do.call(rbind, out)

  ##code as separate vars
  data$CellCycle_DEAD <- NA
  data$CellCycle_DEAD <- ifelse(data$CellCycle == "DEAD", 1, 0)
  data$CellCycle_G1 <- NA
  data$CellCycle_G1 <- ifelse(data$CellCycle == "G1", 1, 0)
  data$CellCycle_G2S <- NA
  data$CellCycle_G2S <- ifelse(data$CellCycle == "G2/S", 1, 0)
  
  doParallel::stopImplicitCluster()

  return (data)
}
