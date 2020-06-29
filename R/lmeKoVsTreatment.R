#' @title LME Kontrolle vs Treatment
#' 
#' @description Use VERSUCH as random factor
#' 
#' @param vars can be all or a list of selected colnames
#' @param adjust which posthoc adjustment should be used?
#' Can be "tukey", standard is "none"
#' @param koVar a
#' @param koGrep b
#' 
#' @import nlme
#' @import lsmeans
#' @import plyr
#' @import parallel
#' @import foreach
#' @import doParallel
#' 
#' @export
#' 
#' @param data
#' 
lmeKoVsTreatment <- function(data, vars="all", adjust="none", koVar="Ko_1", koGrep="Ko") {
  #parallelization
  no_cores <- parallel::detectCores() - 1
  no_cores <- ifelse(no_cores == 0, 1, no_cores)
  no_cores <- ifelse(is.na(no_cores), 1, no_cores)
  doParallel::registerDoParallel(no_cores)
  
  
  ## Type III 
  options(contrasts=c("contr.sum", "contr.poly"))
  
  ##Welche Variablen sollen untersucht werden?
  if (vars == "all") {
    vars <- colnames(data)
  } 
  
  ret <- foreach (var=vars) %dopar% {
    df <- data.frame(VAR=data[,var], 
                     TREATMENT=data$TREATMENT, #WAS TREATMENT.1 -> Check! #TODO
                     VERSUCH=data$VERSUCH, 
                     WELL=data$WELL)
    
    ##Numerische werte?
    df2 <- NULL
    if (!all(is.na(df$VAR)) && class(df$VAR) == "numeric" && sd(df$VAR, na.rm=T) > 0)  {
      df$TREATMENT <- relevel(df$TREATMENT, koVar) #TODO -> as parameter
      
      mod <- lme(VAR~TREATMENT, data=df, random= ~1 | factor(VERSUCH))
      posthoc <- lsmeans(mod, pairwise~factor(TREATMENT), data=df, adjust=adjust)
      
      ph <- summary(posthoc$contrasts)
      ph <- ph[which(grepl(koGrep, ph[,1])),]
      
      df2 <- data.frame(ph, PARAM=var, VGL=ph[,1])
    }
    #print(df2)
  }
  doParallel::stopImplicitCluster()
  
  
  return(do.call(rbind, ret))
}





