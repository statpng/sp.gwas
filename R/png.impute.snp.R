png.impute.snp <- function(xx){
  # xx <- c("NN", "NA", "AA", "AT", "AA", "AA", "AA")
  # xx <- c("NN", "AA", "AA", "AA", "AA", "AA", "AA")
  # xx <- c("AA", "AA", "AA", "AA", "AA", "AA", "AA")
  # xx <- c("AT", "AT", "AT", "AT", "AT", "AT", "AT")
  # xx <- c("AT", "AT", "AT", "TA", "TA", "AT", "NN")
  # xx <- c("AA", "TT", "AA", "TT", "TT", "TT", "NN")
  # xx <- c("AA", "TT", "AA", "TT", "TT", "GG", "NN")
  NonGenotype <- c("-", "_2", "NN", "00", "--", "//", "++", "XX")
  
  if( any( class(xx) %in% c("data.frame", "list") ) ) xx <- as.character(unlist(xx))
  xx <- ifelse(xx %in% NonGenotype, NA, xx)
  
  if( all(is.na(xx)) ) return(xx)
  
  # xx <- sapply( xx, function(a) paste0(sort(strsplit(a, "")[[1]]), collapse="") )
  x.na <- xx[is.na(xx)]
  x.value <- xx[!is.na(xx)]
  
  if( sum(is.na(xx))==0 ) return(xx)
  if( length(x.value)==0 ) return(x.na)
  
  tb <- table(unlist(strsplit(x.value, "")))
  if(length(tb)>2) stop("There are three or more types of alleles. ")
  alleles <- names(tb)
  major.allele <- names( tb[which.max(tb)] )
  minor.allele <- alleles[ !alleles %in% major.allele ]
  
  if(length(minor.allele)==0){
    xx[is.na(xx)] <- unique(x.value)
    return(xx)
  }
  
  tf.hetero <- !sapply( strsplit( names(table(xx)), "" ), function(x) length(unique(x))==1 )
  if( any(tf.hetero) ){
    tmp.comb <- apply( expand.grid( alleles, alleles ), 1, function(x) paste0(x, collapse="" ) )
    rev.hetero <- paste0( rev( strsplit(names(table(xx))[tf.hetero], "")[[1]] ), collapse = "" )
    combs <- tmp.comb[ !tmp.comb %in% rev.hetero ]
  } else {
    combs <- unique( apply( expand.grid( alleles, alleles ), 1, function(x) paste0(sort(x), collapse="" ) ) )
  }
  
  
  tb <- table( factor( x.value, levels=combs ) )
  ord.tb <- order( sapply( strsplit( names(tb), "" ), function(x) sum(x==minor.allele) ) )
  tb.new <- tb[ord.tb]
  
  y <- sum( tb.new[2] + 2*tb.new[3] )
  n <- 2*sum(tb.new)
  # maf <- (tb[2]+2*tb[3])/(2*sum(tb))
  # pi(p) ~ beta(0.5, 0.5)
  # L(y|p) ~ B(p)
  # pi(p|y) \prop pi(p) * L(y|p)
  #         ~ Beta(y+0.5, n+0.5-y)
  maf <- rbeta(n=length(x.na), shape1=y+0.5, n+0.5-y )
  # curve(dbeta(x, shape1=y+2, n+2-y ))
  impute.value <- rbinom(n=length(x.na), size = 2, prob = maf)
  xx[is.na(xx)] <- names(tb.new)[ impute.value+1 ]
  xx
}
