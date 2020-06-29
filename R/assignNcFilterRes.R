#' @title Assign ncAreaFilter Results
#'
#' @description des
#'
#' @import parallel
#' @import foreach
#' @import doParallel
#'
#' @export
assignNcFilterRes <- function(data, calc, type="CImedian",
                              varC="AreaShape_Area.cell", 
                              varN="AreaShape_Area.nucl", 
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
  
  ##by experiment / plate
  data$ncArea <- NA

  out <- foreach(lv=lvPV) %dopar% { 
    dataSub <- data[which(paste(data$VERSUCH, data$PLATTE) == lv),]

    print(paste("Processing", lv))
    cc <- do.call(rbind, calcDF[which(paste(calcDF[,"versuch"], calcDF[,"platte"]) == lv),1][[1]])
    ## by treatment
    for (tr in unique(unlist(cc[,"treatment"]))) {
      #print(tr)
      val <- cc[which(cc[,"treatment"] == tr),]
      sel <- which(dataSub$TREATMENT == tr)
      rat <- dataSub[sel,varC]/dataSub[sel,varN]
      
      ##Type of cellcycle assignment
      if (type=="CIexcl") {
        ##cells with ratio below cutoff
        dataSub$ncArea[sel] <- ifelse(rat < val$estim[["minLeftPos"]], 0, dataSub$ncArea[sel]) 
        ##cells - ok
        dataSub$ncArea[sel] <- ifelse(rat > val$estim[["minLeftPos"]], 1, dataSub$ncArea[sel]) 
      } else if (type == "CImedian") {
        ##cells with ratio below cutoff
        dataSub$ncArea[sel] <- ifelse(rat < val$estim[["minLeftPos"]], 0, dataSub$ncArea[sel]) 
        dataSub$ncArea[sel] <- ifelse(rat >= val$estim[["minLeftPos"]], 1, dataSub$ncArea[sel]) 
      } else {
        stop("Unknown type! Can be CIexcl or CImedian!")
      }
    }
    dataSub
    #print(length(which(is.na(data$ncArea))))
  }
  data <- do.call(rbind, out)
  doParallel::stopImplicitCluster()

  return (data)
}
