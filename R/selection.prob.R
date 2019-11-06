#' @import CMplot
#' @import dplyr
#' @import ggplot2
#' @import glmnet
#' @import utils
#' @importFrom bestNormalize boxcox
#' @importFrom compiler cmpfun
#' @importFrom data.table fread
#' @importFrom gridExtra grid.arrange
#' @importFrom gtools mixedorder
#' @importFrom readxl read_xlsx
#' @importFrom tidyr gather
#' @importFrom tidyr spread
#' @importFrom tidyr unite
#' @importFrom writexl write_xlsx
selection.prob <- function(x=myGD, y=myY, save.path, snp.info=myGM,
                           method=c("lasso", "enet"),
                           family="gaussian",
                           Falsediscovery=c(1,5,10),
                           lambda.min.quantile=0.5,
                           n.lambda=10,
                           K=100,
                           psub=0.5,
                           setseed=1234,
                            ...){

    seq.alpha <- switch( method, "lasso" = 1, "enet" = seq(0.1, 0.9, 0.1) )

    sp <- matrix(NA, nrow=ncol(x), ncol=ncol(y)-1 )
    threshold <- matrix(NA, nrow=length(Falsediscovery), ncol=ncol(y)-1 )
    colnames(sp) <- colnames(threshold) <- colnames(y[,-1,drop=FALSE])
    rownames(threshold) <- Falsediscovery

    for( j in 1:(ncol(y)-1)){
        print(paste0("For the ",j,"-th phenotype,"))
        yj <- y[,j+1,drop=FALSE]
        res <- sp.glmnet( x = x, y = yj,
                          family=family, K=K,
                          seq.alpha=seq.alpha,
                          lambda.min.quantile=lambda.min.quantile,
                          n.lambda=n.lambda,
                          psub=psub,
                          setseed = setseed,
                          ...)

        sp[, j] <- apply( res$sp, 1, max )

        threshold[, j] <- sapply( sort(Falsediscovery, decreasing=FALSE),
                                  function(false_discovery) sp.threshold(res$sp,
                                                                         FD=false_discovery))
    }

    sp.df <- data.frame(snp.info[-1,c(1,3,4)], sp, stringsAsFactors=FALSE) %>%
        .[,!apply(., 2, function(X) all(is.na(X)) )]
    colnames(sp.df) <- c("rs", "chr", "pos", colnames(y[,-1,drop=FALSE]) )

    sp.res <- list(sp.df = sp.df, threshold = threshold)
    save( sp.res, file=paste0(save.path,"/[2]sp.res", ".RData") )

    write_xlsx( sp.df, path=paste0(save.path,"/[2]sp.results.xlsx") )
    write_xlsx( data.frame(FD=Falsediscovery, threshold, stringsAsFactors=FALSE), path=paste0(save.path,"/[2]sp.thresholds.xlsx") )

    return( sp.res )
}
