#' @title Calculate cutoffs for filters
#' 
#' @export
calcCutoffs <- function(obj) {
    # FIXME -> alles hardcoded!
    ##TODO: use only specific filters

    dataCC <- list()    
    dataFB <- list()    
    dataNC <- list()    
    varCC <- "Intensity_IntegratedIntensity_DNA.nucl"    
    varFB="Intensity_MedianIntensity_DNA.nucl"    
    xMinGlobMaxFB=0.2    
    xMaxGlobMaxFB=0.6    

    for (v in 1:length(obj@experiment)) {    
	for (p in 1:length(obj@experiment[[v]])) {    
	    print(paste("VERSUCH: ", v, " PLATTE: ", p))    
	    base <- obj@experiment[[v]][p]    

	    #########################    
	    ###CC    
	    file <- paste(base, "Primarieswithoutborder.txt", sep="")    
	    nrow <- getNRows(file)        
	    meta <- readSingleCol(file, "Metadata_Well",nrow=nrow, type="character")             
	    dat <- readSingleCol(file, "Intensity_IntegratedIntensity_DNA",nrow=nrow)       
	    ### FIXME: well name        
	    treat <- read.csv(paste0(obj@experiment[[v]][p], "Treatment.csv"))    
	    treat <- treat[match(meta, treat$well),]    
	    treat$TREATMENT <- apply(treat, 1, function(x) paste(x[-c(1,2)], collapse="_"))   
	    df <- data.frame(val=dat, WELL=meta, VERSUCH=v, PLATTE=p, TREATMENT=treat$TREATMENT)
	    colnames(df)[1] <- varCC
	    dataCC[[length(dataCC)+1]] <- cellCycleFractIntegrDNAInt(df, log = T, xMinGlobMax = log(100), var = varCC)

	    ########################
	    ### FB
	    dat <- readSingleCol(file, "Intensity_MedianIntensity_DNA",nrow=nrow)   
	    ### FIXME: well name    
	    df <- data.frame(val=dat, WELL=meta, VERSUCH=v, PLATTE=p, TREATMENT=treat$TREATMENT)
	    colnames(df)[1] <- varFB
	    dataFB[[length(dataFB)+1]] <- cellCycleFractIntegrDNAInt(df, var=varFB, xMinGlobMax=xMinGlobMaxFB, xMaxGlobMax=xMaxGlobMaxFB, cutUpper=1)

	    ######################
	    ### NC
	    datNucl <- readSingleCol(file, "AreaShape_Area",nrow=nrow)   
	    metaN <- meta
	    metaN_ObjId <- readSingleCol(file, "ObjectNumber",nrow=nrow)
	    metaN_ImageNumber <- readSingleCol(file, "ImageNumber",nrow=nrow)
	    dfN <- data.frame(val=datNucl, WELL=metaN, ON=metaN_ObjId, IN=metaN_ImageNumber,
			      VERSUCH=v, PLATTE=p, TREATMENT=treat$TREATMENT)
	    colnames(dfN)[1] <- "AreaShape_Area.nucl"

	    file <- paste(base, "Cells.txt", sep="")
	    datCell <- readSingleCol(file, "AreaShape_Area",nrow=nrow)   
	    metaC <- readSingleCol(file, "Metadata_Well",nrow=nrow, type="character")         
	    metaC_ObjId <- readSingleCol(file, "ObjectNumber",nrow=nrow, type="character")         
	    metaC_ImageNumber <- readSingleCol(file, "ImageNumber",nrow=nrow, type="character")         
	    dfC <- data.frame(val=datCell, WELL=metaC, ON=metaC_ObjId, IN=metaC_ImageNumber,
			      VERSUCH=v, PLATTE=p, TREATMENT=treat$TREATMENT)
	    colnames(dfC)[1] <- "AreaShape_Area.cell"
	    df <- merge(dfC, dfN, by=c("VERSUCH", "PLATTE", "WELL", "TREATMENT", "ON", "IN"))
	    dataNC[[length(dataNC)+1]] <- ncAreaFilter(df)

	}
    }

    cutoffs <- list(dataCC=dataCC, dataNC=dataNC, dataFB=dataFB)
    return(cutoffs)

    #save(cutoffs, file=paste0("../analysis/cutoffs_",uid))
}
