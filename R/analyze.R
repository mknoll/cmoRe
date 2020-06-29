#' @title Analyze data
#' 
#' @description TAVI data analysis, consinstig of a linear mixed
#' effect modelling approach with cross validation. 
#' 
#' @param data data.frame as obtained by the aggregateByWell() 
#' or zTransform() functions
#' @param frm Formula to use with lme4::lmer
#' @param frm0 Null-model formula
#' @param ref Reference level of the GRP variable
#' 
#' @import foreach
#' @import parallel
#' @import doParallel
#' @import lme4
#'
#' @export
analyze <- function(data, 
                        frm=as.formula(VAL~GRP+DOSE+(1|VERSUCH/PLATTE/TREATMENT)),
                        frm0=as.formula(VAL~1+(1|VERSUCH/PLATTE/TREATMENT)),
                        ref="Ko_1", nCores=NULL, trSep="_"
                        ) {
    ## Check input
    if (!"VERSUCH" %in% colnames(data)) { stop("VERSUCH Colum required in data!") }
    if (!"TREATMENT" %in% colnames(data)) { stop("TREATMENT Colum required in data!") }
    if (!"PLATTE" %in% colnames(data)) { stop("PLATTE Colum required in data!") }
    if (!"DOSE" %in% colnames(data)) { stop("DOSE Colum required in data!") }
    if (!"GRP" %in% colnames(data)) { stop("GRP Column required in data!") }
    if (!ref %in% levels(factor(data$TREATMENT))) { stop(paste(ref, "not a level in TREATMENT")) }

    ## Musltiple cores
    no_cores <- nCores
    if (is.null(nCores)) {
        no_cores <- parallel::detectCores() - 1
        no_cores <- ifelse(no_cores == 0, 1, no_cores)
        no_cores <- ifelse(is.na(no_cores), 1, no_cores)
    } 
    doParallel::registerDoParallel(no_cores)

    ##nur relevante variablen untersuchen
    ## TODO: check!
    columns  <- which(!colnames(data) %in% c("TREATMENT", "PLATTE", "WELL", "VERSUCH", "GRP") & !grepl("Metadata", colnames(data)) & !colnames(data) %in% c("Group.1", "ImageNumber", "ObjectNumber"))

    ## Get treatments and dose
    lvs <- levels(factor(unlist(do.call(rbind, strsplit(as.character(data$TREATMENT), trSep))[,1])))
    lvs <- lvs[which(lvs != unlist(strsplit(ref, trSep)[[1]][1]))]
    #FIXME!!!
    dose <- levels(factor(unlist(do.call(rbind, strsplit(as.character(data$TREATMENT), trSep))[,2])))
    coll <- list()
    for (lv in lvs) {
        cat(paste("\r                   ", round(which(lv == lvs)/length(lvs)*100, 2), "%", sep=""))
        sel <- which(data$TREATMENT %in% ref | data$TREATMENT %in% paste(lv, dose, sep=trSep))

        res <- foreach(varIndex=columns) %dopar% {
            cat(paste("\r  ", round(which(varIndex == columns)/length(columns)*100, 2), "%", sep=""))
            #print(colnames(data)[varIndex])
            df <- data.frame(GRP=data$GRP[sel], 
                             PLATTE=data$PLATTE[sel],
                             VAL=data[sel,varIndex], 
                             VAR=colnames(data)[varIndex],
                             DOSE=data$DOSE[sel],
                             TREATMENT=data$TREATMENT[sel],
                             VERSUCH=factor(data$VERSUCH[sel]))

            out <- NULL
            tryCatch({
                ## use ML for LRT 
                fit0 <- lmer(frm0, data=df, REML=F)
                fit <- lmer(frm, data=df, REML=F)
                a <- anova(fit, fit0)

                ## use REML as reference
                fit <- lmer(frm, data=df, REML=T)
                R2CV <- NULL

                ##cross validate / per Well
                ##  cross validation: leaving out one treatment level measure from one experiment (e.g. Pe_1)
                cvLv <- levels(factor(paste(data$TREATMENT[sel], data$VERSUCH[sel])))
                pTot <- foreach (lv2=cvLv) %dopar% { 
                    predData <- data[sel,,drop=F][which(paste(data$TREATMENT[sel], data$VERSUCH[sel]) == lv2), ,drop=F]
                    dfPred <- data.frame(GRP=predData$GRP,
                                         PLATTE=predData$PLATTE,
                                         VAL=predData[,varIndex], 
                                         VAR=colnames(predData)[varIndex],
                                         DOSE=predData$DOSE,
                                         TREATMENT=predData$TREATMENT,
                                         VERSUCH=predData$VERSUCH)

                    trainData <- data[sel,,drop=F][which(paste(data$TREATMENT[sel], data$VERSUCH[sel]) != lv2), ,drop=F]
                    dfTrain <- data.frame(GRP=trainData$GRP,
                                          PLATTE=trainData$PLATTE,
                                          VAL=trainData[,varIndex], 
                                          VAR=colnames(trainData)[varIndex],
                                          DOSE=trainData$DOSE,
                                          TREATMENT=trainData$TREATMENT,
                                          VERSUCH=trainData$VERSUCH)

                    #######
                    tryCatch({
                        fit1 <- lmer(frm, data=dfTrain)
                        prdTmp <- predict(fit1, newdata=dfPred,  allow.new.levels=T) 
                    }, error=function(e) {
                        prdTmp <- rep(NA, length(dfPred[,1]))
                    })
                    prdTmp
                }
                pTot <- unlist(pTot)
                R2CV <- 1- (sum(predict(fit)-pTot)^2 / length(pTot) / var(predict(fit)))
                R2CV <- ifelse(R2CV < 0, 0, R2CV) ##check
                out <- data.frame(summary(fit)$coef[-1,,drop=F], 
                                  VAR=colnames(data)[varIndex], 
                                  aP=a[2,8], R2cv=R2CV,
                                  LV=lv)
                out <- cbind(out, CMP=rownames(out))
            },
            error=function(e) { })
            out
        }
        res <- do.call(rbind, res)
        coll[[length(coll)+1]] <- res
    }
    doParallel::stopImplicitCluster()

    return(coll)
}
