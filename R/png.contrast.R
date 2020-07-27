#' @import pbapply
#' @export png.compare_impute
png.compare_impute <- function(before, after, equalize=TRUE){
  # before <- impute_before
  # after <- impute_after
  
  png.venn2 = function(x,y,duplicated=FALSE){
    inner = intersect(x,y)
    xnyc = x[!x%in%inner]
    ynxc = y[!y%in%inner]
    if(!duplicated){
      inner = unique(inner)
      xnyc = unique(xnyc)
      ynxc = unique(ynxc)
    }
    list( x = xnyc , y = ynxc , inner = inner )
  }
  
  if(equalize){
    id.before <- as.character( unlist( before[1,-(1:11)] ) )
    id.after <- as.character( unlist( after[1,-(1:11)] ) )
    
    snp.before <- before$rs[-1]
    snp.after <- after$rs[-1]
    
    snp.venn <- png.venn2(snp.before, snp.after)
    id.venn <- png.venn2(id.before, id.after)
    
    before <- before[ c(1, 1+which(before$rs[-1] %in% snp.venn$inner) ), 
                      c(1:11, 11+which(id.before %in% id.venn$inner)) ]
    after <- after[ c(1, 1+which(after$rs[-1] %in% snp.venn$inner) ), 
                    c(1:11, 11+which(id.after %in% id.venn$inner)) ]
    
    before <- before[c(1, order(before$rs[-1]) ),
                     c(1:11, 11+order( as.character( unlist( before[1,-(1:11)] ) ) ) )]
    
    after <- after[c(1, order(after$rs[-1]) ),
                     c(1:11, 11+order( as.character( unlist( after[1,-(1:11)] ) ) ) )]
  }
  
  
  
  NonGenotype <- c("-", "_2", "NN", "00", "--", "//", "++", "XX")
  wh.missing <- which( apply( before[-1,-(1:11)], 1, function(xj) any(xj %in% NonGenotype) ) )
  cat("The number of SNPs with missing values =", length(wh.missing), "\n")
  
  out_missing <- 
    pblapply( wh.missing, function(jj){
      Xjj <- as.character( unlist( before[-1,-(1:11)][jj,] ) )
      Xjj_imp <- as.character( unlist( after[-1,-(1:11)][jj,] ) )
      
      tab <- table( Xjj, Xjj_imp )
      
      rbind( REST=as.matrix(tab[rownames(tab)!="NN",,drop=F]),
             MISSING=tab[rownames(tab)=="NN",] )
    })
  
  names(out_missing) <- before$rs[-1][wh.missing]
  
  list( tab = out_missing,
        which = data.frame(before=which( snp.before %in% names(out_missing) ), 
                           after=which( snp.after %in% names(out_missing) ) ) 
        )
  
  
}





#' @export png.compare_numericalization
png.compare_numericalization <- function(before, after, equalize=TRUE){

  # before <- numericalization_before
  # after <- numericalization_after
  # rownames(after) <- after[,1];  after <- after[,-1]
  # equalize=TRUE
  
  before$rs <- gsub("[[:punct:]]", "_", before$rs)
  colnames(after) <- gsub("[[:punct:]]", "_", colnames(after))
  
  
  png.venn2 = function(x,y,duplicated=FALSE){
    inner = intersect(x,y)
    xnyc = x[!x%in%inner]
    ynxc = y[!y%in%inner]
    if(!duplicated){
      inner = unique(inner)
      xnyc = unique(xnyc)
      ynxc = unique(ynxc)
    }
    list( x = xnyc , y = ynxc , inner = inner )
  }
  
  if(equalize){
    
    id.before <- as.character( unlist( before[1,-(1:11)] ) )
    id.after <- as.character( unlist( rownames(after) ) )
    
    snp.before <- as.character( unlist( before$rs[-1] ) )
    snp.after <- as.character( unlist( colnames(after)[-1] ) )
    
    snp.venn <- png.venn2(snp.before, snp.after)
    id.venn <- png.venn2(id.before, id.after)
    
    before <- before[ c(1, 1+which(snp.before %in% snp.venn$inner) ), 
                      c(1:11, 11+which(id.before %in% id.venn$inner)) ]
    after <- after[ which(id.after %in% id.venn$inner), 
                    c(1, 1+which(snp.after %in% snp.venn$inner)) ]
    
    id.before2 <- as.character( unlist( before[1,-(1:11)] ) )
    id.after2 <- as.character( unlist( rownames(after) ) )
    snp.before2 <- as.character( unlist( before$rs[-1] ) )
    snp.after2 <- as.character( unlist( colnames(after)[-1] ) )
    
    before <- before[c(1, 1+order(snp.before2) ),
                     c(1:11, 11+order( id.before2 ) )]
    
    after <- after[ order(id.after2),
                    c(1, 1+order( snp.after2 ) ) ]
    
    
    id.before3 <- as.character( unlist( before[1,-(1:11)] ) )
    id.after3 <- as.character( unlist( rownames(after) ) )
    identical( id.before3, id.after3 )
    
    snp.before3 <- as.character( unlist( before$rs[-1] ) )
    snp.after3 <- as.character( unlist( colnames(after)[-1] ) )
    identical( snp.before3, snp.after3 )
    
  }
  
  
  
  
  print("Checking the numericalization function")
  tab_genotype_num <- pblapply( 1:(nrow(before)-1), function(jj){
    Xjj <- as.character(unlist(before[-1,-(1:11)][jj,]))
    Xjj <- ifelse( Xjj %in% NonGenotype, NA, Xjj )
    tab <- table(Geno=Xjj, Num=after[,-1][,jj] )
    tab_dimnames <- dimnames(tab)
    diag_rev <- (row(tab) + col(tab)) == (nrow(tab) + 1)
    
    if( length(tab) == 0 ) return(tab)
    
    if( all( tab[ which(!diag_rev, arr.ind=TRUE) ] == 0 ) ){
      diag_trans <- ifelse(diag_rev, 1, 0 )
      tab <- diag_trans %*% tab
      rownames(tab) <- rev(tab_dimnames[[1]])
      colnames(tab) <- tab_dimnames[[2]]
    }
    names(dimnames(tab)) <- NULL
    tab
  })
  
  
  print( head( tab_genotype_num ) )
  false_Numericalization <- sapply( tab_genotype_num, function(tab) sum( tab[row(tab)>col(tab)] ) != 0 )
  cat("The number of falsely converted SNPs in GAPIT.Numericalization =", sum(false_Numericalization), "\n")
  if( sum(false_Numericalization) > 0 ){
    print( head( tab_genotype_num[false_Numericalization] ) )
  }
  
  
  names(tab_genotype_num) <- snp.before3
  
  list( before = before,
        after = after,
        tab = tab_genotype_num,
        which = which(false_Numericalization) )
  
}






# Example
if(FALSE){
  library(dplyr)
  library(tidyr)
  library(pbapply)
  library(ggplot2)
  library(bestNormalize)
  library(genetics)
  library(HardyWeinberg)
  library(gridExtra)
  source("./R/png.impute.snp.R")
  source("./R/gapit.functions.R")

  genotype_without_NA <- R.utils::loadToEnv("./data/genotype.rda")$genotype
  phenotype <- R.utils::loadToEnv("./data/phenotype.rda")$phenotype
  # genotype_without_NA <- sp.gwas::genotype
  # phenotype <- sp.gwas::phenotype
  genotype_with_NA <- genotype_without_NA
  set.seed(123)
  for( j in 1:100 ){
    genotype_with_NA[-1,-(1:11)][j, sample(10, 2, replace = FALSE)] <- "NN"
  }
  input.type = c("object", "path")[1]
  imputation = TRUE
  QC = TRUE
  callrate.range = c(0.95, 1)
  maf.range = c(1e-3, 1)
  HWE.range = c(0, 1)
  heterozygosity.range = c(0, 1)
  remove.missingY = TRUE
  save.path = "./EXAMPLE_tmp"
  y.id.col = 1
  y.col = 2:4
  normalization = FALSE
  family="gaussian"



  
  # sp.gwas::import.hapmap(genotype = genotype_with_NA,
  #                        phenotype = phenotype,
  #                        input.type = c("object", "path")[1],
  #                        imputation = TRUE,
  #                        QC = TRUE,  # if TRUE, the following QC steps (callrate, maf, HWE, heterozygosity) are conducted.
  #                        callrate.range = c(0.95, 1),
  #                        maf.range = c(1e-3, 1),
  #                        HWE.range = c(0, 1),
  #                        heterozygosity.range = c(0, 1),
  #                        remove.missingY = TRUE,   # if TRUE, the samples with any missing phenotypes are filtered out in all data.
  #                        save.path = "./EXAMPLE_obj_imputation",
  #                        y.id.col = 1,
  #                        y.col = 2:4,
  #                        normalization = FALSE, #if family is not "gaussian", i.e. not continuous variable, normalization should be FALSE
  #                        family="gaussian")

  
  impute_before <- read.csv("./EXAMPLE_obj_imputation/[1]myX.csv")
  impute_after <- read.csv("./EXAMPLE_obj_imputation/[1]Impute_myX.csv")
  
  numericalization_before <- read.csv("./EXAMPLE_obj_imputation/[1]myX.csv")
  numericalization_after <- read.csv("./EXAMPLE_obj_imputation/[1]myGD.csv")
  
}
