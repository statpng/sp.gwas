# x <- rbind( mnormt::rmnorm(100, rep(0, 500), diag(500)), mnormt::rmnorm(100, rep(1, 500), diag(500)) )
# y <- data.frame(y1=rep(0:1, each=100), y2=x%*%rep(0.1,500)+rnorm(200), y3=rbinom(200, 5, prob=rep(c(0.4,0.6), each=100)))
# family <- rep("gaussian", ncol(y))
# type.multivariate <- "unet"
# seq.alpha <- 1*0.1
# seq.lambda <- replicate(5, seq(0.5, 6.0, length.out=10), simplify = F)
# seq.df = floor(seq(1, nrow(y), length.out=10))
# K <- 100
# psub <- 0.5
# 
# 
# fit <- png.cp_glmnet.df(x, y, family,
#                  type.multivariate,
#                  seq.alpha,
#                  seq.df, K=5
#                  )
# 
# fit.threshold <- png.sp_threshold(model = fit$model, out = fit, params = fit$params, nperm=10)

png.get_qhat <- function(array){
  Dim <- dim(array)
  sum(array)/(Dim[2]*Dim[3])
}


png.sp_threshold <- function(model,
                             FD=c(1,5,10),
                             nperm=100,
                             seed=1234,
                             ...){
  
  n <- nrow( params$x )
  p <- ncol( params$x )
  q <- ncol( params$y )
  n.methods <- length(out)
  
  out.perm <- array(NA, dim = c(length(FD), p, q, n.methods, nperm), dimnames = list(FD, 1:p, 1:q, names(out), 1:nperm) )
  Permuted_Theoretical_Threshold <- array(0, dim=c(length(FD), q, n.methods))
  
  Theoretical_Threshold <- array(NA, dim=c(length(FD), q, n.methods))
  
  Empirical_Threshold <- array(0, dim=c(length(FD), q, n.methods) )
  
  dimnames(Theoretical_Threshold) <- 
    dimnames(Permuted_Theoretical_Threshold) <- 
    dimnames(Empirical_Threshold) <- list( FD, paste0("y.", 1:q), names(out) )
  
  print("Permutation get started!")
  perm.params <- params
  pb <- txtProgressBar(min=0, max=nperm, style=3)
  for( perm.i in 1:nperm ){
    # set.seed((nperm+1)*j+perm.i)
    setTxtProgressBar(pb, perm.i)
    
    
    # permutation -------------------------------------------------------------
    set.seed( seed + perm.i )
    perm.params$y <- params$y[sample(n),]
    sp.perm <- do.call( get(model), perm.params )
    
    
    # Empirical threshold -----------------------------------------------------
    for( mm in 1:n.methods ){
      for( fdfd in 1:length(FD) ){
        Empirical_Threshold[fdfd, , mm] <- 
          Empirical_Threshold[fdfd, , mm] + 
          apply( sp.perm$sp[[mm]], 2, function(x) sort(x, decreasing=TRUE)[ FD[fdfd] ] ) / nperm
      }
    }
    
    # Permuted theoretical threshold ------------------------------------------
    # K <- sp.perm$params$K
    # 
    # for( mm in 1:n.methods ){
    #   threshold <- matrix(NA, length(FD), q)
    #   for( colcol in 1:q ){
    #     qhat <- mean( sp.perm$out[[mm]][,,,colcol] ) * p
    #     for( fdfd in 1:length(FD) ){
    #       threshold[fdfd, colcol] <- min( 1, qhat^2/(2*FD[fdfd]*p)+1/2 )
    #     }
    #   }
    #   Permuted_Theoretical_Threshold[,,mm] <- Permuted_Theoretical_Threshold[,,mm] + threshold/nperm
    # }
    
  }
  
  
  
  for( mm in 1:n.methods ){
    threshold <- matrix(NA, length(FD), q)
    for( colcol in 1:q ){
      qhat <- png.get_qhat( out[[mm]][,,,colcol,drop=F] )
      for( fdfd in 1:length(FD) ){
        threshold[fdfd, colcol] <- min( 1, qhat^2/(2*FD[fdfd]*p)+1/2 )
      }
    }
    Theoretical_Threshold[,,mm] <- threshold
  }
  
  
  return( list( theoretical = Theoretical_Threshold, 
                # permuted_theoretical = Permuted_Theoretical_Threshold, 
                empirical = Empirical_Threshold ) )
}
