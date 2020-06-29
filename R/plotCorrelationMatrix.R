#' @title Plot correlation matrix
#' 
#' @param data d
#' @param vars v
#' @param filter f
#' 
#' @import pheatmap
#' 
#' @export 
plotCorrelationMatrix <- function(data, vars="all", selVars=NULL,
                                  filter=c("Meta", "Object", "Parent", 
                                           "PLATTE", "VERSUCH"),
                                  perTreatment=F) {

  for (versuch in levels(factor(data$VERSUCH))) {
    vData <- data[which(data$VERSUCH == versuch),]
    
    ind <- c()
    if (vars=="all") {
      ##get valid vars 
      for (i in 1:length(vData[1,])) {
        if (!is.numeric(vData[,i])) { next }  
        add <- T
        for (j in 1:length(filter)) {
          if (grepl(filter[j], colnames(vData)[i])) {
            add <- F
          }
        }
        if (add) { 
          ind <- c(ind, colnames(vData)[i])
        }
      }
    } else {
      ind <- which(colnames(vData) %in% vars)
    }
    
    cData <- cor(data.matrix(vData[,ind]))
    if (any(is.na(cData))) {
      warning("NAs occured!")
    }
    cData[is.na(cData)] <- 0
    
    anno <- NULL
    if (!is.null(selVars)) {
      anno <- data.frame(SELVAR=factor(ifelse(colnames(cData) %in% selVars, "SEL", "NONSEL")))
      rownames(anno) <- colnames(cData)
    }
    
    anno_colors=list(SELVAR=c("SEL"="tomato", "NONSEL"="lightgray"))

    ##Plot All
    pheatmap(cData, show_rownames=F, show_colnames=F, annotation_col=anno, main=paste("ALL", versuch),
             annotation_colors=anno_colors)
    
    if (perTreatment) {
      ##per treatment
      for (tr in levels(factor(vData$TREATMENT))) {
        subD <- vData[which(vData$TREATMENT == tr),ind]
        cData <- cor(data.matrix(subD))
        if (any(is.na(cData))) {
          warning("NAs occured!")
        }
        cData[is.na(cData)] <- 0
        pheatmap(cData, show_rownames=F, show_colnames=F, annotation_col=anno, main=paste(tr, versuch),
                 annotation_colors=anno_colors)
        
      }
    }
  }
}
