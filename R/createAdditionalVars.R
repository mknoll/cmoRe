#' @title creation of new variables (ratios or differences)
#'  
#' @description creates variables calculated as ratios or differences of two other variables
#' 
#' @param data   a
#' @param type one of "both", "diff", "ratio"
#' @param vars  a
#' @param order which value to use first (nucl/cell; nucl-cell) -> order=c("nucl", "cell")
#' 
#' @export
createAdditionalVars <- function(data, type="both", vars, order=c("nucl","cell")) {  
  dt <- data
  
  cls <- which(grepl(paste("\\.",order[2],sep=""), colnames(dt)))
  ncl <- which(grepl(paste("\\.",order[1],sep=""), colnames(dt)))
  if (length(cls) == 0 || length(ncl) == 0) {
      stop("Could not find nucl / cell variables!")
  }
  
  df <- data.frame(cls=cls, cnCls=colnames(dt)[cls], cn=do.call(rbind, strsplit(colnames(dt)[cls], "\\."))[,1])
  df2 <- data.frame(ncl=ncl, cnNcl=colnames(dt)[ncl], cn=do.call(rbind, strsplit(colnames(dt)[ncl], "\\."))[,1])
  df <- cbind(df, df2[match(df$cn, df2$cn),])
  
  df <- df[which(!is.na(df$cls) & !is.na(df$ncl)),]
  df <- df[which(df$cn %in% vars),]
  
  for (i in 1:length(df[,1])) {
    print(i)
    if (class(dt[,df$cls[i]]) %in% c("numeric","integer") && class(dt[,df$ncl[i]]) %in% c("numeric","integer")) {
      if (type %in% c("both", "ratio")) {
        dt <- cbind(dt, RATIO=dt[,df$ncl[i]]/dt[,df$cls[i]])
        colnames(dt)[length(dt[1,])] <- paste("RATIO_NC_",df$cn[i], sep="")
      }
      
      if (type %in% c("both", "diff")) {
        dt <- cbind(dt, DIFF=dt[,df$ncl[i]]-dt[,df$cls[i]])
        colnames(dt)[length(dt[1,])] <- paste("DIFF_NC_",df$cn[i], sep="")
      }
    }
  }
  
  data <- dt

  return(data)
}
