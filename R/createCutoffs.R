#' @title Calculate cutoffs for filters
#' 
#' @param cutoffPrev returen data from calcCutoffs from previous run
#' 
#' @export
calcCutoffs <- function(obj, cutoffPrev=NULL, delim="\t", rmTreatCol=c(1)) {
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
	    print(paste("VERSUCH: ", v, " PLATTE: ", p))    
	    base <- obj@experiment[[v]][p]    
	    ### counter for comparion with previous data
	    cnt <- cnt+1
	    if (!is.null(cutoffPrev) && 
		length(cutoffPrev$dataCC[[cnt]]) > 0 && !is.na(cutoffPrev$dataCC[[cnt]]) &&
		length(cutoffPrev$dataNC[[cnt]]) >0 && !is.na(cutoffPrev$dataNC[[cnt]]) &&
		length(cutoffPrev$dataFB[[cnt]]) > 0 && !is.na(cutoffPrev$dataFB[[cnt]])) {
		## data intact FIXME!
		dataCC[[cnt]] <- cutoffPrev$dataCC[[cnt]]
		dataNC[[cnt]] <- cutoffPrev$dataNC[[cnt]]
		dataFB[[cnt]] <- cutoffPrev$dataFB[[cnt]]
		next
	    }


	    #########################    
	    ###CC    
	    file <- paste(base, "Primarieswithoutborder.txt", sep="")    
	    nrow <- getNRows(file)        
	    meta <- readSingleCol(file, "Metadata_Well",nrow=nrow, type="character",delim=delim)             
	    dat <- readSingleCol(file, "Intensity_IntegratedIntensity_DNA",nrow=nrow,delim=delim)       
	    if (is.na(nrow) || length(meta) == 1 || length(dat) == 1) {  ## FIXME
		dataCC[[length(dataCC)+1]] <- NA
		dataFB[[length(dataFB)+1]] <- NA
		dataNC[[length(dataNC)+1]] <- NA
		next #### FIXME
	    }
	    ### FIXME: well name        
	    treat <- read.csv(paste0(obj@experiment[[v]][p], "Treatment.csv"))    
	    treat_bak <- treat
	    treat <- treat[match(meta, treat$well),]    
	    treat$TREATMENT <- apply(treat, 1, function(x) paste(x[-rmTreatCol], collapse="_"))   
	    df <- data.frame(val=dat, WELL=meta, VERSUCH=v, PLATTE=p, TREATMENT=treat$TREATMENT)
	    colnames(df)[1] <- varCC
	    dataCC[[length(dataCC)+1]] <- cellCycleFractIntegrDNAInt(df, log = T, xMinGlobMax = log(100), var = varCC)

	    ########################
	    ### FB
	    dat <- readSingleCol(file, "Intensity_MedianIntensity_DNA",nrow=nrow,delim=delim)   
	    if (length(dat) == 1 && is.na(dat)) {
		dataFB[[length(dataFB)+1]] <- NA
		dataNC[[length(dataNC)+1]] <- NA
		next #### FIXME
	    }
	    ### FIXME: well name    
	    df <- data.frame(val=dat, WELL=meta, VERSUCH=v, PLATTE=p, TREATMENT=treat$TREATMENT)
	    colnames(df)[1] <- varFB
	    dataFB[[length(dataFB)+1]] <- cellCycleFractIntegrDNAInt(df, var=varFB, xMinGlobMax=xMinGlobMaxFB, xMaxGlobMax=xMaxGlobMaxFB, cutUpper=1)

	    ######################
	    ### NC
	    datNucl <- readSingleCol(file, "AreaShape_Area",nrow=nrow,delim=delim)   
	    metaN <- meta
	    metaN_ObjId <- readSingleCol(file, "ObjectNumber",nrow=nrow,delim=delim)
	    metaN_ImageNumber <- readSingleCol(file, "ImageNumber",nrow=nrow,delim=delim)
	    if (length(datNucl) == 1 || length(metaN_ObjId) ==1 || length(metaN_ImageNumber) == 1) { 
		dataNC[[length(dataNC)+1]] <- NA
		next #### FIXME
	    }
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
	    treat <- data.frame(apply(treat,2, function(x) trimws(x)))
	    treat$TREATMENT <- apply(treat, 1, function(x) paste(x[-rmTreatCol], collapse="_"))   
	    dfC <- data.frame(val=datCell, WELL=metaC, ON=metaC_ObjId, IN=metaC_ImageNumber,
			      VERSUCH=v, PLATTE=p, TREATMENT=treat$TREATMENT)
	    colnames(dfC)[1] <- "AreaShape_Area.cell"
	    df <- merge(dfC, dfN, by=c("VERSUCH", "PLATTE", "WELL", "TREATMENT", "ON", "IN"))
	    dataNC[[length(dataNC)+1]] <- ncAreaFilter(df)

	}
    }

    cutoffs <- list(dataCC=dataCC, dataNC=dataNC, dataFB=dataFB)
    return(cutoffs)
}
