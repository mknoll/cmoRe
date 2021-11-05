#' @title Calculate cutoffs for filters
#' 
#' @export
calcCutoffs <- function(obj, fun=c("cc", "nc", "fb"), delim="\t", rmTreatCol=c(1), prev=NULL) {
    # FIXME -> alles hardcoded!
    ##TODO: use only specific filters

    dataCC <- list()    
    dataFB <- list()    
    dataNC <- list()    
    varCC <- "Intensity_IntegratedIntensity_DNA.nucl"    
    varFB="Intensity_MedianIntensity_DNA.nucl"    
    xMinGlobMaxFB=0.2    
    xMaxGlobMaxFB=0.6    

    cnt <- 0
    for (v in 1:length(obj@experiment)) {    
	for (p in 1:length(obj@experiment[[v]])) {    
	    cnt <- cnt+1
	    print(paste("VERSUCH: ", v, " PLATTE: ", p))    
	    base <- obj@experiment[[v]][p]    

	    ## fix previously calulcated cutoffs object
	    fixFB<- F
	    fixNC<-F
	    fixCC<-F
	    if (!is.null(prev)) {
		if ("fb" %in% fun) {
		    if (length(prev$dataFB[[cnt]][[1]]$data) > 0) { next }
		    warning("Found missing data: FB!")
		    fixFB <- T
		}
		if ("nc" %in% fun) {
		    if (length(prev$dataFB[[cnt]][[1]]$data) > 0) { next }
		    warning("Found missing data: NC!")
		    fixNC <- T
		}
		if ("cc" %in% fun) {
		    if (length(prev$dataFB[[cnt]][[1]]$data) > 0) { next }
		    warning("Found missing data: CC!")
		    fixCC <- T
		}
	    }

	    ##### TREATMENTS
	    ### FIXME: well name        
	    treatT <- read.csv(paste0(obj@experiment[[v]][p], "Treatment.csv"), nrow=1)
	    treat <- read.csv(paste0(obj@experiment[[v]][p], "Treatment.csv"), colClasses=rep("character", length(treatT[1,])))    
	    treat <- data.frame(apply(treat,2, function(x) trimws(x)))
	    treat$TREATMENT <- apply(treat, 1, function(x) paste(x[-rmTreatCol], collapse="_"))   
	    treat_bak <- treat

	    #########################    
	    ###CC    
	    if ("cc" %in% fun || "fb" %in% fun || "nc" %in% fun) {
		file <- paste(base, "Primarieswithoutborder.txt", sep="")    
		nrow <- getNRows(file)        
		meta <- readSingleCol(file, "Metadata_Well",nrow=nrow, type="character",delim=delim)             
		treat <- treat[match(meta, treat$well),]    
		
		if ("cc" %in% fun) {
		    dat <- readSingleCol(file, "Intensity_IntegratedIntensity_DNA",nrow=nrow,delim=delim)       
		    df <- data.frame(val=dat, WELL=meta, VERSUCH=v, PLATTE=p, TREATMENT=treat$TREATMENT)
		    colnames(df)[1] <- varCC
		    if (!fixCC) {
			dataCC[[length(dataCC)+1]] <- cellCycleFractIntegrDNAInt(df, log = T, xMinGlobMax = log(100), var = varCC)
		    } else {
			prev$dataCC[[cnt]] <- cellCycleFractIntegrDNAInt(df, log = T, xMinGlobMax = log(100), var = varCC)
		    }
		}
	    }

	    ########################
	    ### FB
	    if ("fb" %in% fun) {
		dat <- readSingleCol(file, "Intensity_MedianIntensity_DNA",nrow=nrow,delim=delim)   
		### FIXME: well name    
		df <- data.frame(val=dat, WELL=meta, VERSUCH=v, PLATTE=p, TREATMENT=treat$TREATMENT)
		colnames(df)[1] <- varFB
		if (!fixFB) {
		    dataFB[[length(dataFB)+1]] <- cellCycleFractIntegrDNAInt(df, var=varFB, xMinGlobMax=xMinGlobMaxFB, xMaxGlobMax=xMaxGlobMaxFB, cutUpper=1)
		} else {
		    prev$dataFB[[cnt]] <- cellCycleFractIntegrDNAInt(df, var=varFB, xMinGlobMax=xMinGlobMaxFB, xMaxGlobMax=xMaxGlobMaxFB, cutUpper=1)
		}
	    }

	    ######################
	    ### NC
	    if ("nc" %in% fun) {
		datNucl <- readSingleCol(file, "AreaShape_Area",nrow=nrow,delim=delim)   
		metaN <- meta
		metaN_ObjId <- readSingleCol(file, "ObjectNumber",nrow=nrow,delim=delim)
		metaN_ImageNumber <- readSingleCol(file, "ImageNumber",nrow=nrow,delim=delim)
		treat <- treat[match(meta, treat$well),]    
		dfN <- data.frame(val=datNucl, WELL=metaN, ON=metaN_ObjId, IN=metaN_ImageNumber,
				  VERSUCH=v, PLATTE=p, TREATMENT=treat$TREATMENT)
		colnames(dfN)[1] <- "AreaShape_Area.nucl"

		file <- paste(base, "Cells.txt", sep="")
		nrow <- getNRows(file)        
		datCell <- readSingleCol(file, "AreaShape_Area",nrow=nrow,delim=delim)   
		metaC <- readSingleCol(file, "Metadata_Well",nrow=nrow, type="character",delim=delim)         
		metaC_ObjId <- readSingleCol(file, "ObjectNumber",nrow=nrow, type="character",delim=delim)         
		metaC_ImageNumber <- readSingleCol(file, "ImageNumber",nrow=nrow, type="character",delim=delim)         
		treat <- treat_bak
		treat <- treat[match(metaC, treat$well),]    
		dfC <- data.frame(val=datCell, WELL=metaC, ON=metaC_ObjId, IN=metaC_ImageNumber,
				  VERSUCH=v, PLATTE=p, TREATMENT=treat$TREATMENT)
		colnames(dfC)[1] <- "AreaShape_Area.cell"
		df <- merge(dfC, dfN, by=c("VERSUCH", "PLATTE", "WELL", "TREATMENT", "ON", "IN"))
		if (!fixNC)  {
		    dataNC[[length(dataNC)+1]] <- ncAreaFilter(df)
		} else {
		    prev$dataNC[[cnt]] <- ncAreaFilter(df)
		}
	    }

	}
    }

    if (is.null(prev)) {
	cutoffs <- list(dataCC=dataCC, dataNC=dataNC, dataFB=dataFB)
    } else {
	cutoffs <- prev
    }
    return(cutoffs)
}
