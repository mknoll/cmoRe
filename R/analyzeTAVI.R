#' @title Analyze TAVI data
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
analyzeTAVI <- function(data, 
                        frm=as.formula(VAL~GRP+DOSE+(1|VERSUCH/SUBJ/TREATMENT)),
                        frm0=as.formula(VAL~1+(1|VERSUCH/SUBJ/TREATMENT)),
                        ref="preTAVI"
                        ) {
    ## Check input
    if (!"VERSUCH" %in% colnames(data)) { stop("VERSUCH Colum required in data!") }
    if (!"TREATMENT" %in% colnames(data)) { stop("TREATMENT Colum required in data!") }
    if (!"GRP" %in% colnames(data)) { stop("GRP Colum required in data!") }
    if (!"SUBJ" %in% colnames(data)) { stop("SUBJ Colum required in data!") }
    if (!ref %in% levels(factor(data$GRP))) { stop(paste(ref, "not a level in GRP")) }

    ## Multiple cores
    no_cores <- parallel::detectCores() - 1
    no_cores <- ifelse(no_cores == 0, 1, no_cores)
    no_cores <- ifelse(is.na(no_cores), 1, no_cores)
    doParallel::stopImplicitCluster()
    doParallel::registerDoParallel(no_cores)

    ##nur relevante variablen untersuchen
    ## TODO: check!
    columns  <- which(!colnames(data) %in% c("TREATMENT", "PLATTE", "WELL", "VERSUCH", "GRP", "SUBJ") & !grepl("Metadata", colnames(data)) & !colnames(data) %in% c("Group.1", "ImageNumber", "ObjectNumber"))

    ## Unteschiede zw. pre/post Tavi
    taviDiff <- list()
    #taviDiff <- foreach(varIndex=columns) %dopar% {
    for (varIndex in columns) {
        print(colnames(data)[varIndex])
        df <- data.frame(GRP=data$GRP, 
                         SUBJ=data$SUBJ, 
                         VAL=data[,varIndex], VAR=colnames(data)[varIndex],
                         DOSE=log(as.numeric(as.character(do.call(rbind, strsplit(as.character(data$TREATMENT), "_"))[,2]))),
                         TREATMENT=data$TREATMENT,
                         VERSUCH=factor(data$VERSUCH))

        # Center Dose Variable
        df$DOSE <- df$DOSE-mean(df$DOSE)

        out <- NULL
        tryCatch({
            fit0 <- lmer(frm0, data=df, REML=F)
            fit <- lmer(frm, data=df, REML=F)
            a <- anova(fit, fit0)

            fit <- lmer(frm, data=df, REML=T)

            #### cross validation
            R2CV <- NULL
            ##cross validate / per Well
            #https://www.researchgate.net/post/Does_R_code_for_k-fold_cross_validation_of_a_nested_glmer_model_exist
            ##  cross validation: leaving out one treatment level measure from one experiment (e.g. Pe_1)
            cvLv <- levels(factor(paste(data$TREATMENT, data$VERSUCH)))

            pTot <- foreach (lv2=cvLv) %dopar% { 
                predData <- data[which(paste(data$TREATMENT, data$VERSUCH) == lv2), ,drop=F]
                dfPred <- data.frame(GRP=predData$GRP,
                                     SUBJ=predData$SUBJ,
                                     VAL=predData[,varIndex], VAR=colnames(predData)[varIndex],
                                     DOSE=log(as.numeric(as.character(do.call(rbind, strsplit(as.character(predData$TREATMENT), "_"))[,2]))),
                                     TREATMENT=predData$TREATMENT,
                                     VERSUCH=predData$VERSUCH)
                dfPred$DOSE <- dfPred$DOSE-mean(df$DOSE)

                trainData <- data[which(paste(data$TREATMENT, data$VERSUCH) != lv2), ,drop=F]
                dfTrain <- data.frame(GRP=trainData$GRP,
                                      SUBJ=trainData$SUBJ,
                                      VAL=trainData[,varIndex], VAR=colnames(trainData)[varIndex],
                                      DOSE=log(as.numeric(as.character(do.call(rbind, strsplit(as.character(trainData$TREATMENT), "_"))[,2]))),
                                      TREATMENT=trainData$TREATMENT,
                                      VERSUCH=trainData$VERSUCH)
                dfTrain$DOSE <- dfTrain$DOSE-mean(df$DOSE)
                dfTrain$GRP <- relevel(dfTrain$GRP, ref)

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
                              aP=a[2,8], R2cv=R2CV)
            out <- cbind(out, CMP=rownames(out))
        },
        error=function(e) { })
	taviDiff[[length(taviDiff)+1]] <- out
    }
    ret <- do.call(rbind, taviDiff)

    doParallel::stopImplicitCluster()

    return(ret)
}
