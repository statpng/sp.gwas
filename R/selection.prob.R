#' @export selection.prob
selection.prob <- function(x=myGD, y=myY, pop=NULL,
                           save.path, snp.info=myGM,
                           method=c("lasso", "enet"),
                           family="gaussian",
                           false.discovery=c(1,5,10),
                           lambda.min.quantile=0.5,
                           permutation=FALSE,
                           alpha.seq=NULL,
                           nperm=100,
                           n.lambda=10,
                           K=100,
                           psub=0.5,
                           seed,
                           ...){
    
    n <- nrow(x)
    p <- ncol(x)
    q <- ncol(y)-1
    
    
    if( method == "lasso" ){
        alpha.seq <- 1
    }
    
    
    
    if( !is.null(pop) ){
        pop.level <- length( table(pop) )
        penalty.factor <- c( rep(0, pop.level), rep(1, p) )
    } else {
        pop.level <- 1
    }
    
    if( pop.level > 1 ){
        pop.mat <- model.matrix( ~.-1, data=as.data.frame( as.factor(pop) ) )
        
        x_with_pop <- cbind.data.frame( pop.mat, x )
    } else {
        x_with_pop <- x
        penalty.factor <- rep(1, p)
    }
    
    
    
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
        res <- sp.glmnet( x = x_with_pop, y = yj,
                          family=family, K=K,
                          alpha.seq=alpha.seq,
                          lambda.min.quantile=lambda.min.quantile,
                          n.lambda=n.lambda,
                          psub=psub,
                          # setseed = setseed,
                          verbose=TRUE,
                          penalty.factor=penalty.factor,
                          seed = seed,
                          ...)
        
        sp[, j] <- apply( res$sp, 1, max )
        
        threshold[, j] <- sapply( sort(false.discovery, decreasing=FALSE),
                                  function(false_discovery) sp.threshold(res$sp,
                                                                         FD=false_discovery))
        
    }
    
    sp.df <- data.frame(snp.info[-1,c(1,3,4)], sp, stringsAsFactors=FALSE) %>%
        .[,!apply(., 2, function(X) all(is.na(X)) )]
    colnames(sp.df) <- c("rs", "chr", "pos", paste0("Y", seq_len(ncol(y)-1) ) )
    
    if( permutation ){
        threshold.df <- data.frame(Method="Theoretical", FD=false.discovery, threshold, stringsAsFactors=FALSE)
        write.csv( threshold.df, file=paste0(save.path,"/[2]sp.thresholds_before_perm.csv"), row.names=FALSE, fileEncoding = "UTF-8" )
    }
    
    write.csv( sp.df, file=paste0(save.path,"/[2]sp.results.csv"), row.names=FALSE, fileEncoding = "UTF-8" )
    
    
    for( j in 1:(ncol(y)-1)){
        yj <- y[,j+1,drop=FALSE]
        
        if( permutation ){
            
            print("Permutation started!")
            pb <- txtProgressBar(min=0, max=nperm, style=3)
            for( perm.i in 1:nperm ){
                # set.seed((nperm+1)*j+perm.i)
                setTxtProgressBar(pb, perm.i)
                res.perm <- sp.glmnet( x = x_with_pop, y = sample(yj),
                                       family=family, K=K,
                                       alpha.seq=alpha.seq,
                                       lambda.min.quantile=lambda.min.quantile,
                                       n.lambda=n.lambda,
                                       psub=psub,
                                       # setseed=setseed,
                                       verbose=FALSE,
                                       penalty.factor=penalty.factor,
                                       seed=seed,
                                       ...)
                
                sp.perm[,j,perm.i] <- apply( res.perm$sp, 1, max )
                threshold.perm[, j] <- threshold.perm[, j] + sort( sp.perm[,j,perm.i], decreasing=TRUE )[false.discovery]/nperm
            }
            print(paste0("For the ", j, "th outcome, permutation was done! \n"))
            
        }
        
    }
    
    
    
    
    if(permutation){
        save( sp.perm, file=paste0(save.path,"/[2]sp.perm.results.RData") )
    }
    
    
    sp.res <- list(sp.df = sp.df, threshold = threshold)
    if( permutation ) sp.res <- append(sp.res, list(threshold.perm = threshold.perm))
    
    threshold.df <- data.frame(Method="Theoretical", FD=false.discovery, threshold, stringsAsFactors=FALSE)
    if( permutation ) threshold.df <- rbind.data.frame(threshold.df, 
                                                       data.frame(Method="Permuted", FD=false.discovery, threshold.perm, stringsAsFactors=FALSE) )
    
    save( sp.res, file=paste0(save.path,"/[2]sp.res", ".RData") )
    
    write.csv( threshold.df, file=paste0(save.path,"/[2]sp.thresholds.csv"), row.names=FALSE, fileEncoding = "UTF-8" )
    
    return( sp.res )
}
