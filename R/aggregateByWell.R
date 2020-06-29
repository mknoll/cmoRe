#' @title Aggregate data per well
#' 
#' @param data
#' 
#' @export 
#' @return dataframe containing aggregated data per well
#' 
aggregateByWell <- function (data, fun=median, onlyCellNumber=F) {
    ######################################
    data$WELL <- data$Metadata_Well
    ##cells/well
    cnt<- rep(1, length(data[,1]))
    print("Calculate numbers of cells/well ...")
    cnAgg <- aggregate(cnt, 
		       by=list(paste(data$VERSUCH, data$PLATTE, data$WELL, data$TREATMENT, sep=":")), 
		       function(x) sum(x, na.rm=T))
    cn <- data.frame(CELLNUMBER=cnAgg[,2], do.call(rbind, strsplit(as.character(cnAgg[,1]), ":")))
    colnames(cn)[-1] <- c("VERSUCH", "PLATTE", "WELL", "TREATMENT")
    ## LEdiglich Anzahl Zellen pro Well bestimmen?
    if (onlyCellNumber) {
	return(cn)
    }

    ######################################
    ##find all numeric & exclude metadata & exclude cutoff
    num <- which(is.numeric(as.numeric(as.character(data[1,]))) 
		 & !grepl("Meta", colnames(data)) 
		 & !colnames(data) %in% c("ImageNumber", "ObjectNumber",
					  "TREATMENT", "VERSUCH", "WELL",
					  "PLATTE", "ncArea") 
		 & !grepl("CellCycle", colnames(data))
		 & !grepl("GRP_CUTOFFVAR", colnames(data)))
    print("Aggregate numeric data ...")
    agg <- aggregate(data[,num], 
		     by=list(paste(data$VERSUCH, data$PLATTE, data$WELL, data$TREATMENT, sep=":")), 
		     fun)
    ##########################################
    ## add VERSUCH, PLATTE, WELL, TREATMENT
    tmp <- do.call(rbind, strsplit(as.character(agg[,1]), ":"))
    colnames(tmp) <- c("VERSUCH", "PLATTE", "WELL", "TREATMENT")
    agg <- cbind(agg[,-1], tmp)

    ###########################################
    ## add cell number
    agg <- cbind(agg, CELLNUMBER=cnAgg[,2])

    ## add cell cycle fractions
    if (any(grepl("CellCycle", colnames(data)))) {
	print("Count cell cycle fractions ...")
	ccAgg <- aggregate(data[,which(grepl("CellCycle_", colnames(data)))], 
			   by=list(paste(data$VERSUCH, data$PLATTE, data$WELL, data$TREATMENT, sep=":")),
			   function(x) sum(x, na.rm=T))
	colnames(ccAgg)[-1] <- colnames(data[,which(grepl("CellCycle_", colnames(data)))])
	agg <- cbind(agg, ccAgg[,-1])
	## calculate Fractions
	ccFrac <- ccAgg[,-1]/agg$CELLNUMBER
	colnames(ccFrac) <- paste(colnames(data[,which(grepl("CellCycle_", colnames(data)))]), "_FRAC", sep="")
	agg <- cbind(agg, ccFrac)
    } else {
	warning("Could not detect CellCycle Data!")
    }

    ############################################
    ## add ncAreaFiltered
    if (any(grepl("ncArea", colnames(data)))) {
	print("Count N/C Area ratio filtered cells ...")
	ncAgg <- aggregate(data[,which(grepl("ncArea", colnames(data)))], 
			   by=list(paste(data$VERSUCH, data$PLATTE, data$WELL, data$TREATMENT, sep=":")),
			   function(x) sum(x, na.rm=T))
	agg <- cbind(agg, ncArea=ncAgg[,-1])
	agg$ncArea_FRAC <- agg$ncArea/agg$CELLNUMBER
    } else {
	warning("Could not detect NC/Area Filtered Data!")
    }

    #############################################
    ## add Fibroblast filtered
    if (any(grepl("FIBROBLAST", colnames(data)))) {
	print("Count fibroblasts ...")
	fbAgg <- aggregate(data[,which(grepl("FIBROBLAST", colnames(data)))], 
			   by=list(paste(data$VERSUCH, data$PLATTE, data$WELL, data$TREATMENT, sep=":")),
			   function(x) sum(x, na.rm=T))
	agg <- cbind(agg, fibro=fbAgg[,-1])
	agg$FIBROBLAST_FRAC <- agg$fibro/agg$CELLNUMBER
    } else {
	warning("Could not detect Fibroblast filtered Data!")
    }

    ############################################
    ## add cutoff variables
    if (any(grepl("GRP_CUTOFFVAR", colnames(data)))) {
	print("Count coutoff [hihg/low] cells ...")
	w <- which(grepl("GRP_CUTOFFVAR_", colnames(data)) & grepl("BINARY", colnames(data)))
	cutoffAgg <- aggregate(data[,w], 
			       by=list(paste(data$VERSUCH, data$PLATTE, data$WELL, data$TREATMENT, sep=":")),
			       function(x) sum(x, na.rm=T))
	colnames(cutoffAgg)[-1] <- colnames(data[,w])
	agg <- cbind(agg, cutoffAgg[,-1])
	## calculate Fractions
	cutoffFrac <- cutoffAgg[,-1]/agg$CELLNUMBER
	colnames(cutoffFrac) <- paste(colnames(data[,w]), "_FRAC", sep="")
	agg <- cbind(agg, cutoffFrac)
    } else {
	warning("Could not detect Cutoff Data!")
    }

    return(agg)
}
