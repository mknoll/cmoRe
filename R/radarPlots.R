#' @title Draw Radarplots
#' 
#' @param data specifies the data for plotting
#' @param vars specifies which variables are used for plotting and the order of plotting
#' @param agg statistical method of data aggregation (median or mean)
#' @param main Title
#' @import grDevices
#' 
#' @return no return value
#' 
#' @export
drawRadarplots <- function(data, vars, agg="median", 
                           res=250, file="", width=1000, height=1000,
                           labels=F, col=NULL, ctrlLevel="Ko_1",
                           main=NULL, pPos=1.05, pSym="*",
                           pTest="t.test", pAdj="none",
                           pCol=list("0.05"="red", "0.1"="blue"),
                           pCex=3, na.rm=T, scale=T) {
  ## Test auf Referenzlevel
  if (!ctrlLevel %in% levels(factor(data$TREATMENT))) {
    stop("Reference level does not exist! Please provide via ctrlLevel.")
  }
  if (!all(c("VERSUCH", "TREATMENT") %in% colnames(data))) {
      stop("VERSUCH and/or TREATMENT columns are missing!")
  }
  
  ## Testen?
  if (is.null(pTest)) {
    test <- F
  } else {
    test <- T
  }
  
  ## Aggregiere daten pro Versuch und Treatment
  #agg <- function(x) { agg(x[which(!is.na(x))]) }
  print("1. lksdhfdsf")
  aggData <- aggregate(data[,which(colnames(data) %in% vars)], by=list(data$VERSUCH, data$TREATMENT), function(x) agg(x, na.rm=T))
  print("lksdhfdsf")
  grps <- aggData[,c(1:2),drop=F]
  aggData <- aggData[,-c(1:2),drop=F]
  print(aggData)####
  if (scale) {
     ## Skaliere Daten fuer Starplot
     metrScaled <- apply(aggData, 2, function(x) (x-min(x, na.rm=na.rm))/(max(x, na.rm=na.rm)-min(x, na.rm=na.rm)))
  } else {
     metrScaled <- aggData
  }
  
  ###  Data for the control Level
  ctrlData <- NULL
  subMetr <- metrScaled[which(grps$Group.2 == ctrlLevel),,drop=F]
  subMetr <- subMetr[,which(colnames(subMetr) %in% vars),drop=F]
  subMetr <- rbind(subMetr, apply(subMetr, 2, agg))
  subMetr <- subMetr[,match(vars, colnames(subMetr)),drop=F]
  subMetr <- subMetr[length(subMetr[,1]),,drop=F]
  ctrlData <- subMetr

  ############################################
  
  ##FIXME
  r <- seq(0, 1, 1/length(levels(factor(data$TREATMENT))))
  i <- 1
 
  ret <- list()
  ## Plot fuer jedes Treatment gegen Kontrolle
  for (tr in levels(factor(data$TREATMENT))) {      
      ## Selektiere Treatment
      subMetr <- metrScaled[which(grps$Group.2 == tr),,drop=F]
      ctrlMetr <- metrScaled[which(grps$Group.2 == ctrlLevel),,drop=F]
      
      ## Nur gegebene Variablen
      subMetr <- subMetr[,which(colnames(subMetr) %in% vars),drop=F]
      ctrlMetr <- ctrlMetr[,which(colnames(ctrlMetr) %in% vars),drop=F]
      
      #Teste auf Unterschiede gegen Kontrolle
      if (test) {
        grp <- c(rep("TR", length(subMetr[,1])), rep("CTRL", length(ctrlMetr[,1])))
        pVal <- apply(rbind(subMetr, ctrlMetr), 2, function(x) {
          ifelse(pTest == "t.test",  t.test(x~grp)$p.value, wilcox.test(x~grp)$p.value)
        })
        names(pVal) <- colnames(subMetr)
      }
      
      ##Aggregiere daten
      subMetr <- rbind(subMetr, apply(subMetr, 2, agg))
      
      ##Sortierung anpassen an vorgegebene Reihenfolge
      if (test) {
        pVal <- pVal[match(vars, names(pVal))]
      }
      subMetr <- subMetr[,match(vars, colnames(subMetr)),drop=F]
      
      ## P value adjustment
      if (test) {
        pValAdj <- p.adjust(pVal, pAdj)
      }
      
      ## Nur die aggregierten Werte plotten
      subMetr <- subMetr[length(subMetr[,1]),,drop=F]
      
      ## Add Control Level Values 
      if (tr != ctrlLevel) {
        subMetr <- rbind(ctrlData, subMetr)
      }
    
      #######################
      cl <- rep(rgb(0.3,0.3,0.3,1), length(subMetr[,1]))
      if (is.null(col)) {
        cl[length(subMetr[,1])] <- rgb(r[i],0.6,0.1,0.5)
      } else {
        cl <- col[[i]]
      }
      i<-i+1
      
      ##Titel umbauen # FIXME - REMOVE
      
      if (file != "") {
        png(paste(file,"__",tr,".png",sep=""), res=res, width=width, height=height)
      }
      #umgebuaten titel als ueberschrift verwenden
      if (!is.null(main)) {
        tr <- main
      }
      
      
      if (labels) {
        stars(subMetr, locations = c(0, 0), radius = FALSE, key.loc = c(0, 0), lty = 1, scale=F, lwd=3,
              col.stars=cl, col.lines=cl, xlim=c(-1,1), ylim=c(-1,1), main = tr, cex=0.7) 
      } else {
        stars(subMetr, locations = c(0, 0), radius = FALSE, key.loc = c(0, 0), lty = 1, scale=F, lwd=3,
              col.stars=cl, col.lines=cl, xlim=c(-1,1), ylim=c(-1,1), key.labels=NULL, main = tr) 
      }
      
      ###################################
      ## polarkoordinaten / pvalue color
      if (test) {
        par(xpd=T)
        deg <- 2*pi/length(subMetr[1,])
        rad <- pPos
        cut <- as.numeric(as.character(names(pCol)))
        pCol <- pCol[order(cut)]
        cut <- as.numeric(as.character(names(pCol)))
        cl <- rep("white", length(pValAdj))
        for (u in 1:length(pValAdj)) {
          cl[u] <- pCol[which(pValAdj[u] < cut)[1]]
        }
        cl <- unlist(lapply(cl, function(x) ifelse(is.null(x), "white", x[[1]])))
        for (z in 1:length(subMetr[1,])) {
          phi <- (z-1)*deg
          x <- rad*cos(phi)
          y <- rad*sin(phi)
          text(x,y,pSym, font=2, col=cl[z], cex=pCex)
        }
        par(xpd=F)
      }
      
      
      if (file != "") {
        dev.off()
      }

      ret[[length(ret)+1]] <- list(subMetr, tr)
  }
  
  return(ret)
}
