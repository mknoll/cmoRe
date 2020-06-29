#' @title A Intersect lme/cor results
#' 
#' @description A
#' 
#' @param corR Ergebnisse der corBetweenRep() Funktion
#' @param lmeR Ergebnisse der lmeKoVsTreatment() Funktion
#' @param corCut Cutoff for a minimum required correlation (pearson)
#' @param p-Value Cutoff for lme derived pval (maximum)
#' 
#' @export
intersectLmeCor <- function(corR, lmeR, corCut=NULL, lmePCut=NULL, nRes=NULL) {
  ## Ergebnisse des linear mixed effects models
  trKDf <- do.call(rbind, strsplit(levels(factor(lmeR$VGL)), "-"))
  trKDf <- cbind(trKDf,levels(factor(lmeR$VGL)))
  trKDf <- cbind(trKDf, do.call(rbind, strsplit(trKDf[,1], "_")))
  trKDf <- cbind(trKDf, do.call(rbind, strsplit(trKDf[,2], "_")))
  trKDf <- cbind(trKDf, SUBST=ifelse(trKDf[,4] == "Ko", trKDf[,6], trKDf[,4]))
  trKDf <- apply(trKDf, 2, "trimws")
  
  ## Ergbenisse korrelation
  corDf <- do.call(rbind, corR)
  
  intCorMod <- list()
  for (i in 1:length(trKDf[,1])) {
    trK <- trKDf[i,3]
    
    ##LME
    sub <- lmeR[which(lmeR$VGL == trK),]
    sub <- sub[order(sub$p.value),]
    sub <- sub[which(sub$p.value <= lmePCut),]
    
    ##COR
    corI <- which(corDf[,3] == trKDf[i,8])
    corDat <- cors[[corI]]$cor
    corDatTmp <- corDat[which(corDat >= corCut)]
    int <- intersect(sub$PARAM, names(corDatTmp))
    
    ##best gerankte 
    if (!is.null(nRes)) {
      if (length(int) > 0) {
        to <- ifelse(length(int) > nRes, nRes, length(int))
        int <- int[1:to]
      }
    }
    
    
    intCorMod[[length(intCorMod)+1]] <- 
      list(S1=trK, S2=cors[[corI]]$Substanz, 
           pval=lmePCut, cor=corCut, vars=int, 
           nVars=length(int), nRes=nRes)
  }
  
  return(do.call(rbind, intCorMod))
}



#' @title Int Cor dyn
#' @description dfs
#' 
intersectLmeCorDynamic <- function(corR, lmeR, minCorCut=0.2, lmePCut=0.05, n=10) {
  
  ## Ergebnisse des linear mixed effects models
  trKDf <- do.call(rbind, strsplit(levels(factor(lmeR$VGL)), "-"))
  trKDf <- cbind(trKDf,levels(factor(lmeR$VGL)))
  trKDf <- cbind(trKDf, do.call(rbind, strsplit(trKDf[,1], "_")))
  trKDf <- cbind(trKDf, do.call(rbind, strsplit(trKDf[,2], "_")))
  trKDf <- cbind(trKDf, SUBST=ifelse(trKDf[,4] == "Ko", trKDf[,6], trKDf[,4]))
  trKDf <- apply(trKDf, 2, "trimws")
  
  ## Ergbenisse korrelation
  corDf <- do.call(rbind, corR)
  
  intCorMod <- list()
  for (i in 1:length(trKDf[,1])) {
    trK <- trKDf[i,3]
    
    ##LME; fixed
    sub <- lmeR[which(lmeR$VGL == trK),]
    sub <- sub[order(sub$p.value),]
    sub <- sub[which(sub$p.value <= lmePCut),]
    
    ##COR; dynamic
    corI <- which(corDf[,3] == trKDf[i,8])
    corDat <- cors[[corI]]$cor
    
    corCur <- 1
    delta <- 0.001
    int <- c()
    while (corCur >= minCorCut && length(int) <= n) {
      corCur <- corCur-delta
      corDatTmp <- corDat[which(corDat > corCur)]
      int <- intersect(sub$PARAM, names(corDatTmp))
    }

    
    intCorMod[[length(intCorMod)+1]] <- 
      list(S1=trK, S2=cors[[corI]]$Substanz, 
           pval=lmePCut, cor=corCur, vars=int, 
           nVars=length(int), nRes=n)
  }
  
  return(do.call(rbind, intCorMod))
}
