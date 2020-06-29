#' @title Class assignment
#' 
#' @description  desc
#' 
#' @param fun Aggregation per row
#' 
#' @import randomForest
#' 
#' @export
classAssignment <- function(testData, trainData, mapVars) {
  coll <- list()
  for (kat in unique(names(mapVars))) {
    if (is.null(mapVars[[kat]])) { 
      print(paste("Skipping", kat))
      next 
    }
    
    #FIXME
    ##Select relevant treatment from train data
    trainDat <- trainData[which(grepl(paste(kat, "_", sep=""),trainData$TREATMENT)),] 
    train <- calculateSyntheticVars(trainDat, mapVars[[kat]])
    tra <- apply(train$data,2, function(x) median(x, na.rm=T))

    test <- calculateSyntheticVars(testData, mapVars[[kat]])
    tes <- apply(test$data,2, function(x) median(x, na.rm=T))

    coll[[length(coll)+1]] <- list(kat=kat, SS=sqrt(sum((tra-tes)^2, na.rm=T))/length(tra))
  }
  coll <- do.call(rbind, coll)
  
  coll <- data.frame(apply(coll, 2, "unlist"))
  colnames(coll) <- c( "REF", "DIST")
  coll$DIST <- as.numeric(as.character(coll$DIST))
  
  return(coll)
}


