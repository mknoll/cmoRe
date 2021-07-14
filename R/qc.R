#' @title Create QC plots
#' 
#' @description Creates QC plots
#' 
#' @param data data.frame containing the cell profiler data to check. 
#' Expects the columns TREATMENT, PLATTE, VERSUCH, Metadata_Well
#' @param file pdf filenamename to store the analysis at
#' @param vars Variables to analyze, defaults to nCells (number of cells)
#' @param fun Function to obtain aggregated values per well, e.g. median
#' @param na.rm delete na in fun
#'
#' @import pheatmap
#' @import PerformanceAnalytics
#' 
#' @export
qc <- function(data, folder=NULL, file=NULL, uid=NULL, vars="nCells", fun=NULL, na.rm=T) {
    if (!is.null(folder) && folder == "") { folder <- NULL }
    if (is.null(uid)) { uid <- substr(Sys.time(),1,10) }
    if (!is.null(folder)) { 
        if (is.null(file)) {
	    #FIXME: test for other platforms
	    sep0 <- ifelse(.Platform$OS.type == "unix", "/", "\\")
            file <- unlist(strsplit(tempfile(), sep0))
            file <- file[length(file)]
        }
        #file <- paste(folder, sep0, uid, "__", file, ".pdf", sep="") 
        file <- paste(folder, uid, "__", file, ".pdf", sep="") 
    } else {
        if (is.null(file)) { file <- paste(tempfile(),".pdf", sep="") }
    }
    if (is.null(fun)) { fun <- median }

    pdf(file)
    res <- list()

    ## Treatments
    ov <- aggregate(data$TREATMENT, by=list(paste(data$VERSUCH, data$PLATTE, data$Metadata_Well, sep="|")), function(x) x[1])
    ov <- cbind(do.call(rbind, strsplit(ov[,1], "\\|")), ov[,-1,drop=F])
    colnames(ov)[1:3] <- c("VERSUCH", "PLATTE", "WELL")
    ov$WELL <- as.character(ov$WELL)
    ov$ID1 <- substr(ov$WELL, 1, 1)
    ov$ID2 <- substr(ov$WELL, 2, nchar(ov$WELL))
    treat <- ov

    ### Plate-View
    for (var in vars) {
        if (var == "nCells") {
            ## Number of measurements per well
            ov <- aggregate(data[,1], by=list(paste(data$VERSUCH, data$PLATTE, data$Metadata_Well, data$TREATMENT, sep="|")), function(x) length(x))
        } else {
            ov <- aggregate(data[,var], by=list(paste(data$VERSUCH, data$PLATTE, data$Metadata_Well, data$TREATMENT, sep="|")), fun)
        }
        ov <- cbind(do.call(rbind, strsplit(ov[,1], "\\|")), ov[,-1,drop=F])
        colnames(ov)[1:4] <- c("VERSUCH", "PLATTE", "WELL", "TREATMENT")
        ov$WELL <- as.character(ov$WELL)
        ov$ID1 <- substr(ov$WELL, 1, 1)
        ov$ID2 <- substr(ov$WELL, 2, nchar(ov$WELL))
        res[[length(res)+1]] <- list(data=ov, var=var)
        for (v in unique(ov$VERSUCH)) {
            sub <- ov[which(ov$VERSUCH == v),,drop=F]
            for (p in unique(sub$PLATTE)) {
                subP <- sub[which(sub$PLATTE == p),,drop=F]
                m <- matrix(ncol=12, nrow=8, NA)
                rownames(m) <- c("A", "B", "C", "D", "E", "F", "G", "H")
                colnames(m) <- 1:12
                m2 <- m
		subP$ID2 <- as.numeric(as.character(subP$ID2))
                for (i in 1:length(subP[,1])) {
                    m[subP$ID1[i], subP$ID2[i]] <- subP$x[i]
                    m2[subP$ID1[i], subP$ID2[i]] <- as.character(treat$x[which(treat$VERSUCH == v & treat$PLATTE == p & treat$WELL ==subP$WELL[i])])
                }
                title <- paste0(var, "  V:", v, " P:", p,sep="")
                pheatmap(m, cluster_cols=F, cluster_rows=F, main=title, 
                         cellwidth=10,cellheight=10)
                #print(m)
                textplot(m2, main=title)
                #print(m2)
            }
        }
    }

    dev.off()
    return(list(file, res, treat))
}

