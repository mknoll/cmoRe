#' @title Z-Transform data
#' @description desc
#' 
#' @export
zTransform <- function(data, versuche) {
  if (length(versuche) < 2) {
    stop("Mindestens 2 Versuche notwendig!")
  }
  
  ##Intersect all treatments
  intTr <- levels(droplevels(factor(data$TREATMENT[which(data$VERSUCH == versuche[1])])))
  for (i in 2:length(versuche)) {
    intTr <- intersect(intTr, levels(factor(data$TREATMENT[which(data$VERSUCH == versuche[i])])))
  }
  print("Used treatments:")
  print(intTr)
  
  #TODO: Check if number of measurements are equal in all datasets
  ret <- list()
  
  ##Find numeric columns
  num <- c()
  nn <- c()
  for (i in 1:length(data[1,])) {
    if (is.numeric(data[,i])) {
      num <- c(num, i)
    } else {
      nn <- c(nn, i)
    }
  }

  for (i in 1:length(versuche)) {
    print(versuche[i])
    sel <- which(data$VERSUCH == versuche[i] & data$TREATMENT %in% intTr)
    
    df <- data.frame(scale(data[sel,num]), data[sel,nn])
    df$VERSUCH <- versuche[i]

    ret[[length(ret)+1]] <- df
  }
  
  return(ret)
}
