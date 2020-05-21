selection.prob <- function(x=myGD, y=myY, save.path, snp.info=myGM,
                           method=c("lasso", "enet"),
                           family="gaussian",
                           false.discovery=c(1,5,10),
                           lambda.min.quantile=0.5,
                           permutation=FALSE,
                           nperm=100,
                           n.lambda=10,
                           K=100,
                           psub=0.5,
                           # setseed=1234,
                            ...){

    seq.alpha <- switch( method, "lasso" = 1, "enet" = seq(0.1, 0.9, 0.1) )

    sp <- matrix(NA, nrow=ncol(x), ncol=ncol(y)-1 )
    if(permutation){
        sp.perm <- array(NA, c(ncol(x), ncol(y)-1, nperm))
        dimnames(sp.perm) <- list( NULL, colnames(y[,-1,drop=FALSE]), NULL)
    }
    threshold <- matrix(NA, nrow=length(false.discovery), ncol=ncol(y)-1 )
    colnames(sp) <- colnames(threshold) <- colnames(y[,-1,drop=FALSE])
    rownames(threshold) <- false.discovery
    colnames(threshold) <- paste0("Y", 1:(ncol(y)-1))
    
    if( permutation ){
        threshold.perm <- matrix(0, nrow=length(false.discovery), ncol=ncol(y)-1 )
        rownames(threshold.perm) <- false.discovery
        colnames(threshold.perm) <- paste0("Y", 1:(ncol(y)-1))
    }
    
    for( j in 1:(ncol(y)-1)){
        print(paste0("For the ",j,"-th phenotype,"))
        yj <- y[,j+1,drop=FALSE]
        res <- sp.glmnet( x = x, y = yj,
                          family=family, K=K,
                          seq.alpha=seq.alpha,
                          lambda.min.quantile=lambda.min.quantile,
                          n.lambda=n.lambda,
                          psub=psub,
                          # setseed = setseed,
                          verbose=TRUE,
                          ...)

        sp[, j] <- apply( res$sp, 1, max )

        threshold[, j] <- sapply( sort(false.discovery, decreasing=FALSE),
                                  function(false_discovery) sp.threshold(res$sp,
                                                                         FD=false_discovery))
        
        if( permutation ){
            
            print("Permutation started!")
            pb <- txtProgressBar(min=0, max=nperm, style=3)
            for( perm.i in 1:nperm ){
                # set.seed((nperm+1)*j+perm.i)
                setTxtProgressBar(pb, perm.i)
                res.perm <- sp.glmnet( x = x, y = sample(yj),
                                       family=family, K=K,
                                       seq.alpha=seq.alpha,
                                       lambda.min.quantile=lambda.min.quantile,
                                       n.lambda=n.lambda,
                                       psub=psub,
                                       # setseed=setseed,
                                       verbose=FALSE,
                                       ...)
                
                sp.perm[,j,perm.i] <- apply( res.perm$sp, 1, max )
                threshold.perm[, j] <- threshold.perm[, j] + sort( sp.perm[,j,perm.i], decreasing=TRUE )[false.discovery]/nperm
            }
            print(paste0("For the ", j, "th outcome, permutation was done! \n"))

        }
        
    }
    save( sp.perm, file=paste0(save.path,"/[2]sp.perm.results.RData") )
    

    sp.df <- data.frame(snp.info[-1,c(1,3,4)], sp, stringsAsFactors=FALSE) %>%
        .[,!apply(., 2, function(X) all(is.na(X)) )]
    colnames(sp.df) <- c("rs", "chr", "pos", paste0("Y", seq_len(ncol(y)-1) ) )

    sp.res <- list(sp.df = sp.df, threshold = threshold)
    if( permutation ) sp.res <- append(sp.res, list(threshold.perm = threshold.perm))
    
    threshold.df <- data.frame(Method="Theoretical", FD=false.discovery, threshold, stringsAsFactors=FALSE)
    if( permutation ) threshold.df <- rbind.data.frame(threshold.df, 
                                                       data.frame(Method="Permuted", FD=false.discovery, threshold.perm, stringsAsFactors=FALSE) )
    
    save( sp.res, file=paste0(save.path,"/[2]sp.res", ".RData") )

    write.csv( sp.df, file=paste0(save.path,"/[2]sp.results.csv"), row.names=FALSE, fileEncoding = "UTF-8" )
    write.csv( threshold.df, file=paste0(save.path,"/[2]sp.thresholds.csv"), row.names=FALSE, fileEncoding = "UTF-8" )

    return( sp.res )
}
