#' @title Flatten cutoffs
#' @export
flatten <- function(obj) {
    data <- obj@dataCutoff
    out <- list()

    for (i in 1:length(data)) {
        for (j in 1:length(data[[i]])) {
            tmp <-  data[[i]][[j]]$data
            versuch <- data[[i]][[j]]$versuch
            platte <- data[[i]][[j]]$platte
            for (k in 1:length(tmp)) {
                df <- data.frame(versuch=versuch, platte=platte,
                                 treatment=tmp[[k]]$treatment,
                                 glbMaxEstim=tmp[[k]]$estim['g1MaxPos'],
                                 glbMaxEstimLow=tmp[[k]]$lower['g1MaxPos'],
                                 glbMaxEstimUp=tmp[[k]]$upper['g1MaxPos'],
                                 rightMinEstim=tmp[[k]]$estim['g1MinRightPos'],
                                 rightMinEstimLow=tmp[[k]]$lower['g1MinRightPos'],
                                 rightMinEstimUp=tmp[[k]]$upper['g1MinRightPos'],
                                 nTotal=tmp[[k]]$nTotal,
                                 nHigh=tmp[[k]]$nG2['nG2_median'],
                                 VAR=names(data)[i])
                out[[length(out)+1]] <- df
            }
        }
    }
    return(do.call(rbind, out))
}
