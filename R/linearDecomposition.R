#' @title Linear Decomposition
#' 
#' @description This function correlates the residuals of 
#' 
#' @param data data
#' @param versuche Zu vergleichende Versuche
#' @param vars output of intersectLmeCor
#' @param startvar Startvariable, standard: AreaShape_Area.x
#' @param maxIterations maximum of allowed Iterations
#' @param continueIfNA If the residual cannot be calculated, NA values occur.
#' The calculation continues if set to TRUE. If set to FALSE it stops.
#' 
#' @return data.frame containing bla ...
#' 
#' @export
linearDecomp <- function(data, versuche, vars=NULL, 
                         startvar="AreaShape_Area.x", maxIterations=40,
                         continueIfNA=TRUE) {
  if (length(unique(versuche)) != 2) {
    stop("Genau zwei Versuche muessen angegeben werden!")
  }
  if (is.null(vars)) {
    stop("Bitte vars angeben!")
  }
  collRet <- list()
  
  ##intersectLmeCor output
  mtch <- vars
  
  #TODO -> was TREATMENT.1
  for (tr in levels(factor(unlist(vars[,7])))) {
    print(paste("Treatment:", tr))
    
    substance <- strsplit(tr, "_")[[1]][1]    
    names <- mtch[which(mtch[,2] == substance),5]
    ind <- which(colnames(data) %in% names$vars)
    
    if (length(names$vars) == 0) { next }
    
    ##Startvariable
    selVar <- startvar
    ##Anzahl positiv Korrelierter Variabnle
    posCor <- c()
    
    maxIt <- 0
    while (maxIt < maxIterations) {
      print(paste("Iteration:", maxIt))
      maxIt <- maxIt + 1
      
      ##############################################
      ### Residuen berechnen
      ##############################################
      ##Versuch 1
      rv1 <- getResiduals(data=data[which(data$VERSUCH == versuche[1] & metr$TREATMENT == tr),], vars=names$vars, selVar=selVar)
      if (is.null(rv1)) {
        print("Could not calculate residuals!")
        break
      }
      rv1 <- do.call(rbind, rv1)
      rownames(rv1) <- colnames(data)[ind]
      
      ## Versuch 2
      rv2 <- getResiduals(data=data[which(data$VERSUCH == versuche[2] & metr$TREATMENT == tr),], vars=names$vars, selVar=selVar)
      if (is.null(rv1)) {
        print("Could not calculate residuals!")
        break
      }
      rv2 <- do.call(rbind, rv2)
      rownames(rv2) <- colnames(data)[ind]
      
      
      if (length(rv1) != length(rv2)) {
        warning("Non-matching length of residual-vectors!")
        print("Could not calculate residuals!")
        break
      }
      
      
      ###############################################
      ### Residuen korrelieren
      ###############################################
      corRes <- NULL
      for (i in 1:length(rv1[,1])) {
        sel <- which(!is.na(rv1[i,]) & !is.na(rv2[i,]))
        corRes <- c(corRes, cor(rv1[i,][sel], rv2[i,][sel]))
      }
      names(corRes) <- rownames(rv1)
      
      ## Wie viele negative / positive Korrelationen gibt es?
      grp <- ifelse(corRes > 0, "pos", "neg" )
      tbl <- table(grp)   
      
      ##Abbruch, wenn irgendein cor coef = NA
      if (any(is.na(corRes)) && !continueIfNA) {
        print(selVar)
        print(paste("Iterationen: ", maxIt))
        print("Min. ein NA-Korrelationskoeffizient! STOP!")
        break
      }
      ##Abbruch wenn nur negative Korrelationskoeffizienten 
      if (all(grp[!is.na(grp)] == "neg")) {
        print(selVar)
        print(paste("Iterationen: ", maxIt))
        print("Nur negative Korrelationskoeffizienten! STOP!")
        break
      }
      if (!all(grp[!is.na(grp)] == "pos")) {
        ##Abbruch wenn mehr neg > pos Korrelationsloffizienten
        if (all(is.na(grp))) {
          print(selVar)
          print(paste("Iterationen: ", maxIt))
          print("Nur NA Korrelationskoeffizienten! STOP!")
          break
        }
        
        ##Abbruch wenn mehr neg > pos Korrelationsloffizienten
        if (tbl[[1]] >= tbl[[2]]) {
          print(selVar)
          print(paste("Iterationen: ", maxIt))
          print("Zu viele neg. Korrelationskoeffizienten! STOP!")
          break
        }
      }
      posCor <- c(posCor, length(grp[which(grp == "pos")])/length(grp))
      #posCor <- c(posCor, length(grp))
      #posCor <- c(posCor, length(grp[!is.na(grp)]))
      #posCor <- c(posCor, length(grp[which(grp == "pos")])/length(grp[which( grp == "neg")]))
      
      
      ##Selektion neuer Variable
      it <- 0
      corRes[which(is.na(corRes))] <- corRes[which.min(corRes)]
      corRes <- corRes[rev(order(corRes))]
      while (names(corRes)[1+it] %in% selVar) {
        it <- it+1
      }
      ## weitere Variable hinzufuegen
      if (!is.na(names(corRes)[1+it])) {
        selVar <- c(selVar, names(corRes)[1+it])
      }
    }
    
    collRet[[length(collRet)+1]] <- list(treatment=tr, vars=selVar, posCor=posCor)
  
  }
  
  return (do.call(rbind, collRet))
}


#' @title Get Residuals
#' 
#' @description D
#' 
#' @param data d
#' @param var names of the vars to use
#' @param selVar already selected variables
#' 
#' @return  a
#' 
#' @export
getResiduals <- function(data, vars, selVar) {
  if (length(data[,1]) <= 3) {
    warning("Cannot calculate residuals: less than three observations!")
    return (NULL)
  }
  ##Selecte the variables of interest
  ind <- which(colnames(data) %in% vars)
  
  ## Build the right side of the fitting formula
  form <- ""
  for (j in 1:length(selVar)) {
    if (j == 1) {
      form <- selVar[j]
    } else {
      form <- paste(form, "+", selVar[j])
    }
  }  
  
  ## Get residuals from a linear model fit 
  resTotal <- list()
  for (var in colnames(data)[ind]) {
    formTmp <- paste(var, " ~ ", form)
    formTmp <- as.formula(formTmp)
    

    fit <- lm(formTmp, data=data)
    res <- fit$residuals
    if (is.na(fit$coefficients[length(fit$coefficients)][[1]])) {
      res <- rep(NA, length(res))
    }
    resTotal[[length(resTotal)+ 1 ]] <-  res
  }
  
  return(resTotal)
}


