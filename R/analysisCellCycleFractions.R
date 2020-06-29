#' @title Analysis of cell cycle fractions
#' 
#' @description desc
#' 
#' @param data data
#' 
#' 
#' @export
analysisCCFractions <- function(data, aggregate=F) {
  
  if (!aggregate) {
    #separate for each plate
    for (z in 1:length(data)) {
      tmp <- do.call(rbind, do.call(rbind, data)[,1][[z]])
      tmp <- cbind(tmp, 
                   SUBST=do.call(rbind, strsplit(as.character(tmp[,"treatment"]), "_"))[,1],
                   CONC=as.numeric(do.call(rbind, strsplit(as.character(tmp[,"treatment"]), "_"))[,2])
                   )
      plotFract(tmp, main=paste("V",data[[z]][["versuch"]], " P",data[[z]][["platte"]]))
    }
  } else {
    ##aggregate over plates
    tmpAll <- list()
    for (z in 1:length(data)) {
      tmp <- do.call(rbind, do.call(rbind, data)[,1][[z]])
      tmp <- cbind(tmp, 
                   SUBST=do.call(rbind, strsplit(as.character(tmp[,"treatment"]), "_"))[,1],
                   CONC=as.numeric(do.call(rbind, strsplit(as.character(tmp[,"treatment"]), "_"))[,2],
                   VERSUCH=data[[z]][["versuch"]],
                   PLATTE=data[[z]][["platte"]])
      )
      tmpAll[[length(tmpAll)+1]] <- tmp
    }
    tmpAll <- do.call(rbind, tmpAll)
    plotFract(tmpAll, agg=T)
    return(tmpAll)
  }
}


#' @title Calculate fractions (cellcycle)
#' 
#' @export
calcFract <- function(data, cuts, var="Intensity_IntegratedIntensity_DNA.nucl", ...) {
    data <- data[,which(colnames(data) %in% c(var, "TREATMENT", "VERSUCH", "PLATTE"))]

    plateLv <- unique(paste(data$VERSUCH, data$PLATTE))
    ret <- list()
    pI <- 0
    for (plate in plateLv) {
	pI <- pI+1
	sub <- data[which(paste(data$VERSUCH, data$PLATTE) == plate),,drop=F]
	out <- list()

	lvs <- unique(sub$TREATMENT)
        tI <- 0
	for(trLv in lvs) {
	    tI <- tI+1

	    subTr0 <- sub[which(sub$TREATMENT == trLv),var,drop=F]

	    ##Total cells    
	    nTotal <- length(subTr0[,1]) 
	    if (nTotal != cuts[[pI]]$data[[tI]]$nTotal) {
		stop("Something went terribly wrong! Order of plates / treatments not equal!")
	    }
            ##previously calculated
 	    lower <- cuts[[pI]]$data[[tI]]$lower ###FIXME: imput
 	    upper <- cuts[[pI]]$data[[tI]]$upper #### FIXME: imput
 	    estim <- cuts[[pI]]$data[[tI]]$estim

	    ##Dead cells    
	    nDead <- c(length(which(subTr0[,var] < c(lower[2]))),    
		       length(which(subTr0[,var] < c(upper[2]))),    
		       length(which(subTr0[,var] < c(estim[2]))))    
	    names(nDead) <- paste("nDead", c("lower", "upper", "median"), sep="_")    

	    ##Cells in G2    
	    nG2 <- c(length(which(subTr0[,var] > c(lower[3]))),    
		     length(which(subTr0[,var] > c(upper[3]))),    
		     length(which(subTr0[,var] > c(estim[3]))))    
	    names(nG2) <- paste("nG2", c("lower", "upper", "median"), sep="_")    


	    ##Cells in G1    
	    nG1 <- nTotal-nDead-nG2    
	    names(nG1) <- paste("nG1", c("lower", "upper", "median"), sep="_")    


	    ####
	    out[[length(out)+1]] <- list(treatment=trLv,    
					 nTotal=nTotal,    
					 nDead=nDead,     
					 nG2=nG2,     
					 nG1=nG1,     
					 fracDead=nDead/nTotal,     
					 fracG2=nG2/nTotal,     
					 fracG1=nG1/nTotal,    
					 lower=lower,    
					 upper=upper,    
					 estim=estim,
					widthCI=upper-lower)    
	}
	pv <- do.call(rbind, strsplit(as.character(plate), " "))    
	ret[[length(ret)+1]] <- list(data=out, versuch=pv[1,1], platte=pv[1,2])  
    }
 	return(ret)
}


plotFract <- function(tmp, agg=F, main="") {
  ## stacked graph 
  stck <- data.frame(TREATMENT=unlist(tmp[,"treatment"]),
                     DEAD=unlist(tmp[,"fracDead"])[which(names(unlist(tmp[,"fracDead"])) == "nDead_median")],
                     DEAD_LOW=unlist(tmp[,"fracDead"])[which(names(unlist(tmp[,"fracDead"])) == "nDead_lower")],
                     DEAD_UP=unlist(tmp[,"fracDead"])[which(names(unlist(tmp[,"fracDead"])) == "nDead_upper")],
                     
                     G1=unlist(tmp[,"fracG1"])[which(names(unlist(tmp[,"fracG1"])) == "nG1_median")],
                     G1_LOW=unlist(tmp[,"fracG1"])[which(names(unlist(tmp[,"fracG1"])) == "nG1_lower")],
                     G1_UP=unlist(tmp[,"fracG1"])[which(names(unlist(tmp[,"fracG1"])) == "nG1_upper")],
                     
                     G2=unlist(tmp[,"fracG2"])[which(names(unlist(tmp[,"fracG2"])) == "nG2_median")],
                     G2_LOW=unlist(tmp[,"fracG2"])[which(names(unlist(tmp[,"fracG2"])) == "nG2_lower")],
                     G2_UP=unlist(tmp[,"fracG2"])[which(names(unlist(tmp[,"fracG2"])) == "nG2_upper")]
                     )
  
    ### TODO: do not jhardcode!
  addInfo <- data.frame(do.call(rbind, strsplit(as.character(stck$TREATMENT), "_")))
  colnames(addInfo)[1:2] <- c("SUBST", "CONC")
  addInfo <- data.frame(cbind(TREATMENT=as.character(stck$TREATMENT), addInfo))
  addInfo$CONC <- as.numeric(as.character(addInfo$CONC))
  stck <- stck[order(addInfo$SUBST, addInfo$CONC),]
  
  #par(mar=c(4,4,4,6))
  if (!agg) {
    barplot(data.matrix(t(stck[,c(1,4,7)+1])), ylab="%Cells", names.arg=stck[,1], las=2, main=main)
  } else {
    barplotCI(stck)
    print(stck)
  }
  par(xpd=T)
  legend(length(unique(stck[,1]))+6, 1, c("Dead", "G1", "G2/S"), fill=c("darkgray", "lightgray", "white"),
         cex=0.8, bty='n')
  par(xpd=F)
}


se <- function(x) { return (sd(x)/sqrt(length(x))) }


barplotCI <- function(data, funVar=se, funEW=mean) {
  ##Average over measurements
  data.agg <- aggregate(data[,-1], by=list(data[,"TREATMENT"]), FUN=funEW)
  colnames(data.agg)[1] <- "TREATMENT"
  par(xpd=T)
  barplot(data.matrix(t(data.agg[,c(1,4,7)+1])), ylab="%Cells", names.arg=data.agg[,1], las=2,
          col=c("darkgray", "lightgray", "white"))
  
  ##Errors
  print(data)
  err <- aggregate(data[,-1], by=list(data[,"TREATMENT"]), FUN=funVar)
  
  ## add error bars
  for (i in 1:length(err[,1])) {
    for (j in c(1,2,3)) {
      varInd <- (j-1)*3+2
      
      x <- (i-1)+0.3*j+0.2*(i-1) #0.2 deafult in barplot
      y <- data.agg[i,varInd]
      z<-j
      while (z>1) {
        y <- y+data.agg[i,(varInd-(z-1)*3)]
        z<-z-1
      }
      
      lw <- 1
      #upper
      segments(x,y,x,y+err[i,varInd], lwd=lw)
      segments(x-0.1, y+err[i, varInd], x+0.1, y+err[i, varInd], lwd=lw)
      ##lower
      segments(x,y,x,y-err[i, varInd], lwd=lw)
      segments(x-0.1, y-err[i, varInd], x+0.1, y-err[i, varInd], lwd=lw)
    }
  }
  par(xpd=F)
  
}




