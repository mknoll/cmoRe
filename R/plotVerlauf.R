#' @title Plotte Verlaufskurve
#' 
#' @description Vlerauf
#' 
#' @param column zu plottende Spalte der Daten. 
#' Wenn null, Ausgabe der Zellzahl (=Anzahl Zeilen / Well)
#' @param fun Funktion zur Aggregation der Daten, 
#' Standard: Median
#' 
#' @import sciplot
#' 
#' @export
#' 
#' @return A
plotVerlauf <- function(dat, ylim=NULL, column=NULL, main="", fun=median, lwd=1) {
  if (!is.null(column) && !(column %in% colnames(dat))) {
    stop("Angegebene Spalte [column] nicht in [dat] enthalten!")
  }
  
  ##verschiedene behandlungen
  genTr <- levels(factor(dat$TREATMENT))
  genTr <- cbind(do.call(rbind, strsplit(genTr, "_")), genTr)
  
  ##bestimme ylims
  if (is.null(ylim)) {
    if (is.null(column)) {
      ylim = c(0, 1.05*max(table(paste(dat$TREATMENT, dat$Metadata_Well))))
    } else {
      ylim = c(0, 1.05*unlist(quantile(dat[,column], na.rm=T, 0.8))[[1]][1])
    }
  } 
  
  ## kurven
  add <- F
  i <- 1
  leg <- c()
  cl <- c()
  ty  <- c()
  par(mar=c(3,3,3,8))
  for (l in levels(factor(genTr[,1]))) {
    selTr <- genTr[which(genTr[,1] == l & genTr[,3] != "PE"),4]
    selTr2 <- genTr[which(genTr[,1] == l & genTr[,3] == "PE"),4]
    
    if (is.null(column)) {
      ##Zellzahl
      all <- dat[which(dat$TREATMENT %in% c(selTr)),]
      all$CONC <- as.numeric(do.call(rbind, strsplit(as.character(all$TREATMENT), "_"))[,2])
      lineplot.CI(all$CONC, rep(1, length(all[,1])), all$Metadata_Well,
                  ylim=ylim, fun=sum, add=add, 
                  col=i, legend=F, lwd=lwd, 
                  pch=rep(19, length(levels(factor(paste(all$CONC, all$Metadata_Well))))),
                  main=paste("Zellzahl", main, data[[sel]][3]))
      add <- T 
      
      all <- dat[which(dat$TREATMENT %in% c(selTr2)),]
      all$CONC <- as.numeric(do.call(rbind, strsplit(as.character(all$TREATMENT), "_"))[,2])
      lineplot.CI(all$CONC, rep(1, length(all[,1])), all$Metadata_Well,
                  ylim=ylim, fun=sum, add=add, 
                  col=i, legend=F, lwd=lwd,
                  pch=rep(21, length(levels(factor(paste(all$CONC, all$Metadata_Well))))),
                  main=paste("Zellzahl", main,data[[sel]][3]))
      
    } else {
      #Beliebige Spalten
      all <- dat[which(dat$TREATMENT %in% c(selTr, selTr2)),]
      all$CONC <- as.numeric(do.call(rbind, strsplit(as.character(all$TREATMENT), "_"))[,2])
      all$PE <- do.call(rbind, strsplit(as.character(all$TREATMENT), "_"))[,3]
      all$PE <- factor(all$PE)
      all$PE <- relevel(all$PE, "PE")
      
      
      #lineplot.CI(all$CONC, all$AreaShape_Area.x, all$PE, ylim=yl, fun=median, add=add, col=i, legend=F, lwd=1)
      lineplot.CI(all$CONC, all[,column], all$PE, 
                  ylim=ylim, fun=fun, add=add, 
                  col=i, legend=F, lwd=lwd,
                  main=paste(column, main,data[[sel]][3]))
    }
    
    leg <- c(leg,l , paste(l, "+ PE"))
    cl <- c(cl, rep(i, 2))
    ty  <- c(ty, c(1,2))
    add <- T
    i <- i+1
  }
  
  ##legende
  par(xpd=T)
  legend(5.5, 0.7*ylim[2], leg, fill=cl, lty=ty, cex=0.8, lwd=lwd)
}
