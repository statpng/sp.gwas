#' @importFrom stats quantile
sp.glmnet <- function(x, y,
                      family="gaussian",
                      seq.alpha=NULL,
                      n.lambda=NULL,
                      lambda.min.quantile=0.5,
                      K=100,
                      psub=0.5,
                      setseed, ...){

    if( NROW(y) != nrow(x) ) stop("x and y should be equal length of row")
    if( NCOL(y)>1 & (family!="mgaussian") ) stop("The family should be 'mgaussian'")
    if( NCOL(y)==1 & (family=="mgaussian") ) stop("The family should not be 'mgaussian'")
    if( missing(setseed) ) stop("Since setseed is missed, please enter a value of setseed.")

    if( family=="binomial" ) y <- ifelse( as.numeric(factor(as.vector(as.matrix(y)))) == 1, 0, 1 )

    if(is.null(seq.alpha)) seq.alpha <- 1:9*0.1
    if(is.null(n.lambda)) n.lambda <- 10

    x <- as.matrix(x)
    y <- as.matrix(y)

    n <- nrow(x);
    p <- ncol(x);
    nsub <- n*psub;

    if( family=="binomial" ){
        wc <- which(y==1)
        wt <- which(y==0)
        nc <- floor(length(wc) * psub)
        nt <- floor(length(wt) * psub)
    }

    vector.lambda <- NULL

    for( i in 1:10 ){
        for( j in 1:length(seq.alpha) ){
            if( family=="binomial" ){
                wsub <- c(sample(wc, nc), sample(wt, nt))
            } else {
                wsub <- sample(n, nsub)
            }
            xsub <- x[wsub,]
            ysub <- y[wsub,]
            fitsub <- glmnet(x=xsub, y=ysub, alpha=seq.alpha[j], family=family)
            vector.lambda <- c( vector.lambda, fitsub$lambda )
        }
    }

    lambda.min <- quantile(vector.lambda, probs = lambda.min.quantile)
    lambda.max <- max(vector.lambda)
    seq.lambda <- seq(lambda.min, lambda.max, length.out=n.lambda);

    ncol.y <- length( append(fitsub$beta, NULL) )

    out <- array(0, c(ncol(x), length(seq.lambda), length(seq.alpha)) );
    for( j in 1:length(seq.alpha) ){
        for( i in 1:K ){
            if( i %% 10 == 0 ) print(paste0("iteration=", i, "  alpha=", seq.alpha[j]))

            set.seed( setseed*i )

            if( family=="binomial" ){
                wsub <- c(sample(wc, nc), sample(wt, nt))
            } else {
                wsub <- sample(n, nsub)
            }
            xsub <- x[wsub,];
            ysub <- y[wsub,];

            out.h <- array(0, c(ncol(x), length(seq.lambda), ncol.y) )
            for( h in seq_len(ncol.y) ){
                glmnet.fit <- glmnet(x=xsub, y=ysub, alpha=seq.alpha[j], lambda=seq.lambda, family=family, ... )
                out.h[,,h] <- as.numeric( append( NULL, glmnet.fit$beta )[[h]]!=0 ) ;
            }

                out[,,j] <- out[,,j] + apply(out.h, c(1,2), function(X) ifelse(any(X==1), 1, 0) )

        }
    }

    return( list( sp = out/K, alpha=seq.alpha, lambda=seq.lambda ) );
}
####################################################




sp.threshold <- function(sp, FD){
    if(max(sp)>1 | min(sp)<0) stop("'Sel.prob' should be in [0,1].")
    if(length( dim(sp) ) < 3) stop("'sp' should be a 3-dimensional array.")
    if( which.max(dim(sp))!=1 ) stop("'sp' array should have values for variables in the first dimension.")

    qhat <- sum(sp)/(dim(sp)[2]*dim(sp)[3])
    p <- dim(sp)[1]
    threshold <- qhat^2/(2*FD*p) + 0.5

    if(threshold>1) threshold <- 1
    return( threshold )
}
