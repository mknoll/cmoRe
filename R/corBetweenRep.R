#' @title correlation between replicates
#' 
#' @description Only numeric columns will be used;
#' Irrespective of substance concentration
#' 
#' @param data specifies the input data
#' @param filter drops columns
#' 
#' @export
corBetweenRep <- function(data, filter=c("Meta", "Number", "VERSUCH", "TREAT")) {
  ## Check for 2 levels
  if (!length(levels(droplevels((factor(data$VERSUCH))))) == 2) {
    stop("Exakt zwei Versuche notwendig!")
  }
  
  ## Filter for numeric and vars
  ind <- c()
  for (i in 1:length(data[1,])) {
    add <- T
    if (!is.numeric(data[,i])) { next }
    for (j in 1:length(filter)) {
      if (grepl(filter[j], colnames(data)[i])) {
        add <- F
      }
    }
    if (add) {
      ind <- c(ind, i)
    }
  }
  
  #extract substanz w/o conc
  ## TODO: was TREATMENT.1
  data$SUBSTANCE <- do.call(rbind, strsplit(as.character(data$TREATMENT), "_"))[,1] 
  data$CONC <- as.numeric(do.call(rbind, strsplit(as.character(data$TREATMENT), "_"))[,2]) 
  
  ##Korrelationen
  cors <- list()
  for (substance in levels(factor(data$SUBSTANCE))) {
    subData <- data[which(data$SUBSTANCE == substance),]
    subData <- subData[order(subData$VERSUCH, subData$SUBSTANCE, subData$CONC),]
    
    lvV <- levels(factor(subData$VERSUCH))
    for (v in 1:length(lvV)) {
      subV <- subData[which(subData$VERSUCH == lvV[v]),]
      for (v2 in 2:length(lvV)) {
        if (v2 == v) { next }
        subV2 <- subData[which(subData$VERSUCH == lvV[v2]),]
        ##match
        if (dim(subV)[1] != dim(subV2)[1]) {
          stop("Inkompatible Dimensionen!")
        }
        
        dt <- diag(cor(subV[,ind], subV2[,ind]))
        cors[[length(cors)+1]] <- list(V1=lvV[v], V2=lvV[v2], Substanz=substance, cor=dt)
      }
    }
  }
  
  return (cors)
}






