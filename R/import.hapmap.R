#' @name import.hapmap
#' @title A function to import the hapmap formatted SNP data and the corresponding phenotype data
#' @description Input: Hapmap-formatted SNP data, phenotype data
#' 
#' Output: Matched data files (genotype, numerical, SNP information, QC information, and phenotype) with QC and/or imputation.
#' 
#' @param genotype Either R object or file path can be considered. A genotype data is not a data.frame but a matrix with dimension \code{p} by \code{(n+11)}. It is formatted by hapmap which has (rs, allele, chr, pos) in the first four(1-4) columns, (strand, assembly, center, protLSID, assayLSID, panel, Qcode) in the following seven(5-11) columns. If NULL, user can choose a path in interactive use.
#' @param phenotype Either R object or file path can be considered. A phenotype data is an \code{n} by \code{p} matrix. Since the first some columns can display attributes of the phenotypes, you should enter the arguments, y.col and y.id.col, which represent the columns of phenotypes to be analyzed and the column of sample ID. If NULL, user can choose a path in interactive use.
#' @param input.type Default is "object". If \code{input.type} is "object", obejects of genotype/phenotype will be entered, and if "path", paths of genotype/phenotype will be enterd. If you want to use an object, you have to make sure that the class of each column of genotype data is equal to "character".
#' @param imputation TRUE or FALSE for whether imputation will be conducted.
#' @param QC TRUE or FALSE for whether QC for SNPs will be conducted.
#' @param callrate.range A numeric vector indicating the range of non-missing proportion. Default is c(0, 1).
#' @param maf.range A numeric vector indicating the range of minor allele frequency (MAF) to be used. Default is c(0, 1).
#' @param HWE.range A numeric vector indicating the range of pvalue by Hardy-Weinberg Equillibrium to be used. Default is c(0, 1).
#' @param heterozygosity.range A numeric vector indicating the range of heterozygosity values to be used, because, in some cases, heterozygosity higher than expected indicates the low quality variants or sample contamination. Default is c(0, 1).
#' @param remove.missingY If TRUE, the samples with missing values in phenotype data are removed. Accordingly, the corresponding genotype samples are also filtered out. Default is TRUE.
#' @param y.col The columns of phenotypes. At most 4 phenotypes can be considered, because the plot of them will be fine. Default is 2.
#' @param y.id.col The column of sample ID in the phenotype data file. Default is 1.
#' @param normalization If TRUE. phenotypes are converted to be normal-shape using box-cox transformation when all phenotypes are positive.
#' @param family A family of response variable(phenotype). It is "gaussian" for continuous response variable, "binomial" for binary, "poisson" for count, etc. Now you can use only the same family for the multi phenotypes. For more details, see the function(\code{stats::glm}). Default is "gaussian".
#' @param save.path A save.path which has all output files. If there exists save.path, sp.gwas will check if there is an output file. Note that if there is an output RData file in "save.path", sp.gwas will just load the output files(.RData) in there, thereby not providing the results for new "genotype" and "phenotype".
#' 
#' @details Hardy-Weinberg Equillibrium test was derived from "genetics" package.
#' In imputation process, we first calculate the empirical allele frequencies.
#' If we use a beta distribution as a prior in order to estimate the posterior distribution of allele frequency, 
#' then the posterior distribution of allele frequecy is also beta distribution.
#' Accordingly, we impute the missing values with samples from the posterior distribution.
#'     
#'     
#' @author Kipoong Kim <kkp7700@gmail.com>
#' 
#' @examples
#' genotype <- sp.gwas::genotype # load("genotype.rda")
#' phenotype <- sp.gwas::phenotype # load("phenotype.rda")
#' 
#' # object
#' import.hapmap(genotype = genotype, 
#'               phenotype = phenotype, 
#'               input.type = c("object", "path")[1], 
#'               imputation = FALSE, 
#'               QC = TRUE,  # if TRUE, the following QC steps (callrate, maf, HWE, heterozygosity) are conducted.
#'               callrate.range = c(0.95, 1),
#'               maf.range = c(1e-3, 1),
#'               HWE.range = c(0, 1),
#'               heterozygosity.range = c(0, 1),
#'               remove.missingY = TRUE,   # if TRUE, the samples with any missing phenotypes are filtered out in all data.
#'               save.path = "./EXAMPLE_obj",
#'               y.id.col = 1, 
#'               y.col = 2:4, 
#'               normalization = FALSE, #if family is not "gaussian", i.e. not continuous variable, normalization should be FALSE
#'               family="gaussian")
#'
#'
#'
#' # path
#' 
#' write.table( x = sp.gwas::genotype, file = "./genotype.csv", row.names = FALSE, col.names = FALSE, sep=",")
#' write.table( x = sp.gwas::phenotype, file = "./phenotype.csv", row.names = FALSE, sep="," )
#' 
#' import.hapmap(genotype = "./genotype.csv", 
#'               phenotype = "./phenotype.csv", 
#'               input.type = c("object", "path")[2], 
#'               QC = TRUE,  # if TRUE, the following QC steps (callrate, maf, HWE, heterozygosity) are conducted.
#'               callrate.range = c(0.95, 1),
#'               maf.range = c(1e-3, 1),
#'               HWE.range = c(0, 1),
#'               heterozygosity.range = c(0, 0.5),
#'               remove.missingY = TRUE,   # if TRUE, the samples with any missing phenotypes are filtered out in all data.
#'               save.path = "./EXAMPLE_path",
#'               y.id.col = 1, 
#'               y.col = 2:4, 
#'               normalization = FALSE, #if family is not "gaussian", i.e. not continuous variable, normalization should be FALSE
#'               family="gaussian")
#'
#'
#' @return A folder containing a genomic data set in which the samples of genotype and phenotype data are matched, and that quality control steps can be conducted for genotype data
#' 
#' @import HardyWeinberg
#' @import pbapply
#' @import ggplot2
#' @importFrom data.table fread
#' @importFrom tidyr gather
#' @importFrom gridExtra grid.arrange
#' @importFrom bestNormalize boxcox
#' @importFrom genetics genotype
#' @importFrom genetics HWE.exact
#' @importFrom dplyr %>%
#' @export import.hapmap
import.hapmap <-
  function(genotype = NULL,
           phenotype = NULL,
           input.type = c("object", "path"),
           save.path,
           y.col = NULL,
           y.id.col = 2,
           family = "gaussian",
           normalization = TRUE,
           remove.missingY = TRUE,
           imputation = FALSE,
           QC = TRUE, 
           callrate.range = c(0, 1),
           maf.range = c(0, 1),
           HWE.range = c(0, 1),
           heterozygosity.range = c(0, 1)) {
    
    # library(pbapply)
    # library(genetics); library(HardyWeinberg)
    # library(tidyr); library(bestNormalize); library(ggplot2); library(gridExtra)
    
    # genotype.xlsx <- readxl::read_excel("./EXAMPLE/NEW_Genotype20190708.xlsx")
    # genotype.list <- R.utils::loadToEnv("./EXAMPLE/[data]TOTAL_GLOBAL_KOREA_final.RData")
    # phenotype.list <- R.utils::loadToEnv("./EXAMPLE/[data]PHENOTYPE_final.RData")
    # genotype <- genotype.xlsx
    # genotype <- genotype.list$TOTAL_final
    # phenotype <- phenotype.list$PHENOTYPE2018
    # class(phenotype$FloweringDate) <- "numeric"
    # input.type = "object"
    # save.path <- "EXAMPLE_gapit"
    # y.col = 7
    # y.id.col = 1
    # family = "gaussian"
    # normalization = TRUE
    # remove.missingY = F
    # QC = TRUE
    # imputation = TRUE
    # callrate.range = c(0.95, 1)
    # maf.range = c(0, 1)
    # HWE.range = c(0, 1)
    # heterozygosity.range = c(0, 1)
    NonGenotype <- c("-", "_2", "NN", "00", "--", "//", "++", "XX")
    
    if( !file.exists( paste0(save.path) ) ){
      dir.create(paste0(save.path))
    }
    
    if( any( class(genotype) %in% c("tbl", "tbl_df") ) ){
      genotype <- as.data.frame(genotype)
    }
    
    if( any( class(phenotype) %in% c("tbl", "tbl_df") ) ){
      phenotype <- as.data.frame(phenotype)
    }
    
    # Import a phenotype data -------------------------------------------------
    if( length( dim(phenotype) ) > 0 | input.type=="object" ){
      myY.init <- phenotype
    } else {
      # if(is.null(phenotype)) {
      #   print("Select the path of phenotype data which has sample IDs in its first column.")
      #   phenotype <- choose.files()
      # }
      if( grepl("*.xlsx", phenotype) ){
        myY.init <- read_excel(phenotype)
      } else if( grepl("*.csv|*.txt", phenotype) ){
        myY.init <- as.data.frame( fread(phenotype), stringsAsFactors=FALSE )
      } else {
        myY.init <- as.data.frame( fread(phenotype), stringsAsFactors=FALSE )
      }
    }
    
    
    # Restriction of the number of phenotypes ---------------------------------
    if( ncol(myY.init) > 5 & is.null(y.col) ){
      stop("The number of phenotypes must be equal or less than 4.")
    } else if (!is.null(y.col)){
      myY.init <- myY.init[,c(y.id.col, y.col)]
    }
    # myY.init <- myY.init[!duplicated(myY.init[,y.id.col]),]
    
    if( any(duplicated(myY.init[,1])) ){
      cat( "Duplicated phenotype IDs: {", myY.init[,1][ duplicated(myY.init[,1]) ], "}\n" )
      stop("There are duplicated sample IDs in phenotype data.")
    }
    
    print("Checking whether the class of phenotype corresponds to family")
    if( family == "gaussian" ){
      for( j in 1:(ncol(myY.init)-1) ){
        if( length( table( myY.init[,j+1] ) ) < 10 ) warning(paste0("The category of ", j, "-th phenotype is less than 5. Isn't it categorical variable?"))
        if( class( myY.init[,j+1] ) != "numeric" ) stop("Please confirm classes of your phenotypes.")
      }
    }
    
    
    
    
    # Import a genotype data --------------------------------------------------
    if( length( dim(genotype) ) > 0 | input.type=="object" ){
      myX.init <- genotype
    } else {
      if(is.null(genotype)) {
        print("Select genotype data with hapmap format.")
        genotype <- choose.files()
      }
      
      if( grepl("*.xlsx", genotype) ){
        myX.init <- read_excel(genotype, col_names=FALSE)
      } else if( grepl("*.csv|*.txt", genotype) ){
        myX.init <- as.data.frame( fread(genotype, header=FALSE), stringsAsFactors=FALSE )
      } else {
        myX.init <- as.data.frame( fread(genotype, header=FALSE), stringsAsFactors=FALSE )
      } 
    }
    
    if( any(duplicated(as.character( myX.init[1, 12:ncol(myX.init) ] ) )) ){
      cat( "Duplicated genotype IDs: {", 
           as.character( myX.init[1, 12:ncol(myX.init) ] )[
             duplicated(as.character( myX.init[1, 12:ncol(myX.init) ] ) )
             ],
           "}\n" )
      stop("There are duplicated sample IDs in genotype data.")
    }
    
    
    
    # Reorder the dataset according to the names ------------------------------
    myX.init <- myX.init[,c(1:11, order( myX.init[1, 12:ncol(myX.init) ] )+11)]
    myY.init <- myY.init[order( myY.init[, 1] ), ]
    
    
    # Remove missing values from Y ---------------------------------------------------
    
    if( remove.missingY ){
      print(paste0("Removing the missing values of phenotypes(",
                   sum(apply(myY.init, 1, function(x) any(is.na(x)))),
                   ")."
      ))
      
      myY.init <- myY.init[!apply(myY.init, 1, function(x) any(is.na(x))), ]
    }
    
    # myX.init[,c(1:11, order( myX.init[1, 12:ncol(myX.init) ] )+11)] %>% filter(rs %in% c("rs_1_0006", "rs_2_48660")) %>% .[,1:20]
    
    # Match the genotype with the phenotype -----------------------------------
    print("Matcing samples of genotpe and phenotype")
    taxa.snp <- as.character( myX.init[1,12:ncol(myX.init), drop=TRUE] )
    taxa.phenotype <- as.character( myY.init[, 1] )
    
    cat("IDs for genotype data are \n{", head(taxa.snp, 10), "}\n")
    cat("IDs for phenotype data are \n{", head(taxa.phenotype, 10), "}\n")
    
    ss.x <- which(taxa.snp %in% intersect( taxa.snp, taxa.phenotype ) )
    ss.y <- which(taxa.phenotype %in% intersect( taxa.snp, taxa.phenotype ) )
    
    if( length(ss.x) < length(taxa.snp)/2 ) warning("More half of samples does not match")
    if( length(ss.y) < length(taxa.phenotype)/2 ) warning("More half of samples does not match")
    
    if( length(ss.x) < 10 ) stop("Check your sample IDs")
    if( length(ss.y) < 10 ) stop("Check your sample IDs")
    
    
    if(QC){
    
    # Quality Control ---------------------------------------------------------
    png.maf <- function(genotype.hapmap, sep=""){
      
      # xx: vecor with c(GG, GC, CC, CC, CC)
      # This function will return 1e-22 for the SNP with no minor allele
      
      
      get.maf <- function(xj){
        xj <- as.character(xj)
        
        if(length(unique(xj))<=1) maf <- 0
        
        xj <- ifelse(xj %in% NonGenotype, NA, xj)
        x.na <- xj[is.na(xj)]
        x.value <- xj[!is.na(xj)]
        
        if(length(x.value)>1){
          tb <- table( unlist( unlist( strsplit( x.value, sep ) ) ) )
          maf <- min( prop.table(tb) )
        } else {
          maf <- 0
        }
        
        ifelse( maf>0.5, 1-maf, maf )
      }
      
      out <- pbapply( genotype.hapmap[-1,][,-(1:11)], 1, function(xj) get.maf(xj) )

      out
    }
    
    
    
    
    
    png.HWE <- function(genotype.hapmap){
      
      get.HWE <- function(xj){
        xj <- as.character(xj)
        xjNA <- replace(xj, list = (xj %in% NonGenotype ), NA)
        if( length(table(xjNA)) == 1 ) return("NA")
        input <- strsplit(xjNA,"")
        if( all(is.na(unlist(input))) ) return("NA")
        
        input.genotype <- genotype( sapply( input, paste0, collapse="/") )
        if(length(table(input.genotype))>1 & length(levels(input.genotype))<6){
          out <- HWE.exact(input.genotype)$p.value
        } else {
          out <- "NA"
        }
        out
      }
      
      png.get_allele.mat <- function(genotype.hapmap){
        allele.mat <- 
          genotype.hapmap[-1,-(1:11)] %>% 
          apply(1, function(xj) {
            allels.vec <- unlist( strsplit(names(table(ifelse(xj %in% NonGenotype, NA, xj))), "") ) %>% unique
            alleles <- allels.vec %>% paste0(collapse="/") 
            if( length(allels.vec) == 1 ){
              alleles <- paste0(c(alleles,"X"), collapse="/")
            }
            alleles
            # if( length(alleles) == 1 ){
            #   expand.grid(unlist(strsplit(alleles, "/")), unlist(strsplit(alleles, "/"))) %>% 
            #     apply(1, function(x) paste0(sort(x),collapse="")) %>% unique
            # }
            })
        allele.mat
      }
      
      
      png.get_allele <- function(xj){
        allels.vec <- unlist( strsplit(names(table(ifelse(xj %in% NonGenotype, NA, xj))), "") ) %>% unique
        alleles <- allels.vec %>% paste0(collapse="/") 
        if( length(allels.vec) == 1 ){
          alleles <- paste0(c(alleles,"X"), collapse="/")
        }
        alleles
      }
      
      # system.time({
      #   hwe.mat <- genotype.hapmap
      #   allele.vec <- png.get_allele.mat(hwe.mat)
      #   genotype.x.na.mat <- pbapply(hwe.mat[-1,-(1:11)], 2, function(xj) ifelse(xj %in% c("NN", "00", "--", "//", "++", "XX"), NA, xj) %>% unlist )
      #   out <- MakeCounts(t(genotype.x.na.mat), allele.vec)[,1:3] %>% HWExactMat(verbose=F) %>% .$pvalvec
      # })

      out <- pbapply(genotype.hapmap[-1,-(1:11)], 1, function(xj) {
        xj.na <- ifelse(xj %in% NonGenotype, NA, xj) %>% unlist
        if( all(is.na(xj.na)) | length(unique(xj.na))==1 ) return(0)
        suppressWarnings( MakeCounts(xj.na, png.get_allele(xj.na))[1:3] %>% HWExact(verbose=F) %>% .$pval %>% as.numeric )
      })
      
      # out <- pbapply( genotype.hapmap[-1,][,-(1:11)][1:20,], 1, function(xj) get.HWE(xj) )
      
      out
    }
    
    
    
    
    
    png.heterozygousCalls <- function(genotype.hapmap){
      
      get.heterozygousCalls <- function(xj){
        set.heterozygote <- apply( subset( expand.grid(c("T","C","A","G"), c("T","C","A","G")), Var1 != Var2 ), 1, paste0, collapse="")
        xj <- as.character(xj)
        xjna <- unlist( replace( xj, xj %in% NonGenotype, NA ) )
        if( all(is.na(xjna)) ) return("NA")
        out <- mean( xjna %in% set.heterozygote )
        out
      }
      
      out <- pbapply( genotype.hapmap[-1,][,-(1:11)], 1, function(xj) get.heterozygousCalls(xj) )
      
      out
    }
    
    
    print("Removing the SNPs with missing genotype in all samples")
    wh.allmissing <- pbapply( myX.init[-1,-(1:11)], 1, function(xj){
      all( unique(xj) %in% NonGenotype )
    } )
    cat( "# of SNPs with all missing =", sum( wh.allmissing ), "\n" )
    
    # myX.init <- myX.init[c(1, which(!wh.allmissing)+1),]
    print("Removing the SNPs with only one value")
    wh.onevalue <- pbapply( myX.init[-1,-(1:11)], 1, function(xj){
      ( length(unique(xj)) == 1 )
    } )
    cat( "# of SNPs with only one genotype =", sum( wh.onevalue ), "\n" )
    
    if( sum( wh.onevalue ) > 0 ){
      cat("Example of SNPs with only one value:", "\n")
      head( myX.init[-1,-(1:11)][which(wh.onevalue),] %>% {cbind(.[,1:10], "_"=".", "_"=".", "_"=".", .[,(ncol(.)-10+1):ncol(.)])} )
    }
    
    # myX.init <- myX.init[c(1, which(!wh.onevalue)+1),]
    print("Calculating call rate (the % of samples with a non-missing genotype")
    CallRate <- myX.init[-1,-(1:11)] %>% pbapply(1, function(xj){
      xjna <- unlist( replace( xj, xj %in% NonGenotype, NA ) )
      (length(xjna)-sum(is.na(xjna)))/length(xjna)
    })
    print("Calculating MAF")
    MAF <- png.maf( myX.init )
    print("Calculating pvalue by HWE")
    pvalue.HWE <- png.HWE( myX.init )
    print("Calculating heterozygosity")
    heterozygosity <- png.heterozygousCalls(myX.init)
    
    filter.allmissing <- which( !wh.allmissing )
    filter.onevalue <- which( !wh.onevalue )
    filter.CallRate <- which( CallRate <= max(callrate.range) & CallRate >= min(callrate.range) )
    filter.MAF <- which( MAF <= max(maf.range) & MAF >= min(maf.range) )
    filter.HWE <- which( sapply(pvalue.HWE, function(xx) min(xx, 1)) <= max(HWE.range) & pvalue.HWE >= min(HWE.range) )
    filter.heterozygosity <- which( heterozygosity <= max(heterozygosity.range) & heterozygosity >= min(heterozygosity.range) )
    
    filter.intersect <- intersect( intersect( intersect(filter.CallRate, filter.MAF), filter.HWE), filter.heterozygosity)
    
    # CallRate.final <- CallRate[filter.intersect]
    # MAF.final <- MAF[filter.intersect]
    # pvalue.HWE.final <- pvalue.HWE[filter.intersect]
    # heterozygosity.final <- heterozygosity[filter.intersect]
    
    
    
    write.csv(x=cbind.data.frame(myX.init[-1,1:4], 
                                 CallRate=CallRate,
                                 MAF=MAF,
                                 pvalue_HWE=pvalue.HWE,
                                 heterozygosity=heterozygosity), 
              file=paste0(save.path,"/[1]myQC.csv"), row.names=FALSE, fileEncoding = "UTF-8")
    
    
    }
    
    
    
    # myX.init <- myX.init[c(1, filter.intersect+1), ]
    
    
    # Numericalize the coding of genotypes ------------------------------------
    print("Performing the numericalization procedure for genotpe data.")
    myGD.init <- suppressWarnings( as.data.frame( apply(myX.init[-1,-(1:11)], 1, function(one) png.Numericalization(one, bit=2, Major.allele.zero=T)), stringsAsFactors = FALSE ) )
    
    print("Checking the numericalization function")
    tab_genotype_num <- 
      pbapply::pblapply( 1:floor(nrow(myGD.init)*0.01), function(jj){
      Xjj <- as.character(myX.init[-1,-(1:11)][jj,])
      Xjj <- ifelse( Xjj %in% NonGenotype, NA, Xjj )
      tab <- table(Geno=Xjj, Num=myGD.init[,jj] )
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
    
    
    
    myGM.init <- myX.init[,1:4]
    myGT.init <- myX.init[,c(12:ncol(myX.init))]
    colnames(myGD.init) <- myX.init[-1,1]
    rownames(myGD.init) <- myX.init[1,-c(1:11)]
    
    missing.snp <- apply( myX.init[-1,-(1:11)], 1, function(X) sum(is.na(X) | X %in% NonGenotype ) )
    missing.sample <- apply( myX.init[-1,-(1:11)], 2, function(X) sum(is.na(X) | X %in% NonGenotype ) )
    
    write.csv(x=data.frame(ID=myX.init[-1,1], missing.snp, stringsAsFactors = FALSE), file=paste0(save.path,"/[0]missing.snp.csv"), row.names=FALSE, fileEncoding = "UTF-8")
    write.csv(x=data.frame(ID=as.character(myX.init[1,-(1:11)]), missing.sample, stringsAsFactors = FALSE), file=paste0(save.path,"/[0]missing.sample.csv"), row.names=FALSE, fileEncoding = "UTF-8")
    
    # Saving Excel But its time is too long ----------------------------------
    # for( Excelobj in list(myX, myGD, myGM, myGT) ){
    #     write.xlsx(x = Excelobj,
    #                file = paste0("Data.GAPIT.",".xlsx"),
    #                sheetName=deparse(substitute(Excelobj)),
    #                append=T )
    # }
    # End ---------------------------------------------------------------------
    
    
    if(identical(as.character(myX.init[1,ss.x+11]), as.character(myY.init[ss.y, 1]) )){
      myX <- as.data.frame( myX.init[,c(1:11,ss.x+11)], stringsAsFactors=FALSE )
      myGD <- myGD.init[ss.x,]
      myGT <- myGT.init[,ss.x]
      myGM <- myGM.init
      myY <- myY.init[ss.y,]
    } else {
      xnyc <- setdiff(as.character(myX.init[1,ss.x+11]), as.character(myY.init[ss.y, 1]))
      ynxc <- setdiff(as.character(myY.init[ss.y, 1]), as.character(myX.init[1,ss.x+11]))
      xydifferences <- unique(c(xnyc, ynxc))
      if( length(xydifferences) == 0 ) xydifferences <- "NULL"
      cat("Sample IDs which differ between genotype and phenotype are {", xydifferences, "}\n")
      stop("Sample IDs between genotypes and phenotypes do not match")
    }
    
    
    # myX[ss.x+11,][1:10,1:10]
    # myY[ss.y,][1:10,]
    if( family=="gaussian" ){
      
      myY.original <- myY
      # for( j in 1:(ncol(myY)-1) ){
      #   myY.original[,j+1] <- as.numeric(as.character(myY.original[,j+1]))
      # }
      if( normalization ){
        print("Performing the normalization for phenotypes, but there is no evidence that they can be normalized.")
        
        for( j in 1:(ncol(myY)-1) ){
          myY[,j+1] <- boxcox(myY.original[,j+1], standardize = FALSE)$x.t
        }
      } else {
        myY <- myY.original
      }
      
      norm.pvalue.original <- norm.pvalue <- NULL
      for( j in 1:(ncol(myY.original)-1) ){
        norm.pvalue.original[j] <- round( shapiro.test(myY.original[,j+1])$p.value, 4)
        norm.pvalue[j] <- round( shapiro.test(myY[,j+1])$p.value, 4)
      }
      
      
      devAskNewPage(ask = FALSE)
      par(ask=F)
      if( normalization ){
        
        Hist.y.original <- myY.original %>%
          .[,-1,drop=F] %>%
          gather(PHENOTYPE) %>%
          data.frame(., pvalue=rep(norm.pvalue.original, each=nrow(.)/length(norm.pvalue.original))) %>%
          ggplot(aes(value)) +
          geom_histogram(fill="white", colour="black") +
          xlab("Original Phenotypes") +
          ggtitle(paste0("Histogram of original phenotypes")) +
          theme_bw() +
          facet_wrap(~PHENOTYPE, nrow=2, scales="free") +
          geom_label(x=Inf, y=Inf, aes(label=pvalue),
                     vjust = "inward", hjust = "inward")
        
        Hist.y <- myY %>%
          .[,-1,drop=FALSE] %>%
          gather(PHENOTYPE) %>%
          data.frame(., pvalue=rep(norm.pvalue, each=nrow(.)/length(norm.pvalue))) %>%
          ggplot(aes(value)) +
          geom_histogram(fill="white", colour="black") +
          ggtitle(paste0("Histogram of transformed phenotypes")) +
          theme_bw() +
          facet_wrap(~PHENOTYPE, nrow=2, scales="free") +
          geom_label(x=Inf, y=Inf, aes(label=pvalue),
                     vjust = "inward", hjust = "inward")
        
      } else {
        Hist.y <- myY %>%
          .[,-1,drop=FALSE] %>%
          gather(PHENOTYPE) %>%
          data.frame(., pvalue=rep(norm.pvalue, each=nrow(.)/length(norm.pvalue))) %>%
          ggplot(aes(value)) +
          geom_histogram(fill="white", colour="black") +
          ggtitle(paste0("Histogram of original phenotypes")) +
          theme_bw() +
          facet_wrap(~PHENOTYPE, nrow=2, scales="free") +
          geom_label(x=Inf, y=Inf, aes(label=pvalue),
                     vjust = "inward", hjust = "inward")
      }
      
      if( normalization ){
        Hist.norm <- grid.arrange(Hist.y.original, Hist.y, nrow=1, heights=15)
        ggsave(Hist.norm, filename = paste0(save.path,"/[1]Hist.phenotype.norm.pdf"), width = 10, height = 5)
      } else {
        Hist.norm <- grid.arrange(Hist.y, nrow=1, heights=15)
        ggsave(Hist.norm, filename = paste0(save.path,"/[1]Hist.phenotype.norm.pdf"), width = 10, height = 5)
      }
      
    }
    
    
    
    
    print("Imputation")
    if(imputation){
      myX.impute <- myX
      set.seed(2020)
      myX.impute[-1,-(1:11)] <- 
        pbapply::pbapply(myX[-1,-(1:11)], 1, function(xj){
          as.character(unlist(xj)) %>% png.impute.snp(., NonGenotype)
        }) %>% t
      
      print("Checking the difference between before and after imputation.")
      wh.missing <- which( apply( myX[-1,-(1:11)], 1, function(xj) any(xj %in% NonGenotype) ) )
      
      out_missing <- 
        pbapply::pblapply( wh.missing[1:floor(length(wh.missing)*0.01)], function(jj){
          Xjj <- as.character( myX[-1,-(1:11)][jj,] )
          Xjj_imp <- as.character( myX.impute[-1,-(1:11)][jj,] )
          
          tab <- table( Xjj, Xjj_imp )
          
          rbind( REST=as.matrix(tab[rownames(tab)!="NN",,drop=F]),
                 MISSING=tab[rownames(tab)=="NN",] )
        })
      
      if( length(out_missing)>0 ){
        sink(file=paste0(save.path,"/[0]Impute_BeforeAfter.csv"))
        out_missing
        sink()
      } else {
        print("There is no missing value")
      }
      
      myGD.impute <- suppressWarnings( as.data.frame( apply(myX.impute[-1,-(1:11)], 1, function(one) png.Numericalization(one, bit=2, Major.allele.zero=T)), stringsAsFactors = FALSE ) )
      myGT.impute <- myX.impute[,c(12:ncol(myX.impute))]
      
      print("Checking the Numericalization step after imputation")
      tab_genotype_num_impute <- 
        pbapply::pblapply( 1:floor(nrow(myGD.impute)*0.01), function(jj){
        Xjj <- as.character(myX.impute[-1,-(1:11)][jj,])
        Xjj <- ifelse( Xjj %in% NonGenotype, NA, Xjj )
        tab <- table(Geno=Xjj, Num=myGD.impute[,jj] )
        tab_dimnames <- dimnames(tab)
        diag_rev <- (row(tab) + col(tab)) == (nrow(tab) + 1)
        
        if( all( tab[ which(!diag_rev, arr.ind=TRUE) ] == 0 ) ){
          diag_trans <- ifelse(diag_rev, 1, 0 )
          tab <- diag_trans %*% tab
          rownames(tab) <- rev(tab_dimnames[[1]])
          colnames(tab) <- tab_dimnames[[2]]
        }
        names(dimnames(tab)) <- NULL
        tab
      })
      
      print( head( tab_genotype_num_impute ) )
      false_Numericalization_impute <- sapply( tab_genotype_num_impute, function(tab) sum( tab[row(tab)>col(tab)] ) != 0 )
      cat("The number of falsely converted SNPs in GAPIT.Numericalization after imputation =", sum(false_Numericalization), "\n")
      if( sum(false_Numericalization_impute) > 0 ){
        print( head( tab_genotype_num_impute[false_Numericalization_impute] ) )
      }
      
      
      colnames(myGD.impute) <- myX.impute[-1,1]
      rownames(myGD.impute) <- myX.impute[1,-c(1:11)]
      
      write.csv(x=myX.impute, file=paste0(save.path,"/[1]Impute_myX.csv"), row.names=FALSE, fileEncoding = "UTF-8")
      write.csv(x=data.frame(ID=rownames(myGD.impute), myGD.impute, stringsAsFactors = FALSE), file=paste0(save.path,"/[1]Impute_myGD.csv"), row.names=FALSE, fileEncoding = "UTF-8")
      write.csv(x=myGM, file=paste0(save.path,"/[1]Impute_myGM.csv"), row.names=FALSE, fileEncoding = "UTF-8")
      write.csv(x=myGT.impute, file=paste0(save.path,"/[1]Impute_myGT.csv"), row.names=FALSE, fileEncoding = "UTF-8")

    }
    
    print("Saving the data files")
    
    write.csv(x=myX, file=paste0(save.path,"/[1]myX.csv"), row.names=FALSE, fileEncoding = "UTF-8")
    if(family=="gaussian"){
      write.csv(x=myY.original, file=paste0(save.path,"/[1]myY.original.csv"), row.names=FALSE, fileEncoding = "UTF-8")
    }
    write.csv(x=myY, file=paste0(save.path,"/[1]myY.csv"), row.names=FALSE, fileEncoding = "UTF-8")
    write.csv(x=data.frame(ID=rownames(myGD), myGD, stringsAsFactors = FALSE), file=paste0(save.path,"/[1]myGD.csv"), row.names=FALSE, fileEncoding = "UTF-8")
    write.csv(x=myGM, file=paste0(save.path,"/[1]myGM.csv"), row.names=FALSE, fileEncoding = "UTF-8")
    write.csv(x=myGT, file=paste0(save.path,"/[1]myGT.csv"), row.names=FALSE, fileEncoding = "UTF-8")
    
    if(QC){
    write.csv(x=myX[c(1,1+filter.intersect),], file=paste0(save.path,"/[1]QC_myX.csv"), row.names=FALSE, fileEncoding = "UTF-8")
    write.csv(x=data.frame(ID=rownames(myGD), myGD[,filter.intersect], stringsAsFactors = FALSE), file=paste0(save.path,"/[1]QC_myGD.csv"), row.names=FALSE, fileEncoding = "UTF-8")
    write.csv(x=myGM[c(1, 1+filter.intersect),], file=paste0(save.path,"/[1]QC_myGM.csv"), row.names=FALSE, fileEncoding = "UTF-8")
    write.csv(x=myGT[c(1, 1+filter.intersect),], file=paste0(save.path,"/[1]QC_myGT.csv"), row.names=FALSE, fileEncoding = "UTF-8")
    }
    
    if(QC & imputation){
      write.csv(x=myX.impute[c(1,1+filter.intersect),], file=paste0(save.path,"/[1]ImputeQC_myX.csv"), row.names=FALSE, fileEncoding = "UTF-8")
      write.csv(x=data.frame(ID=rownames(myGD.impute), myGD.impute[,filter.intersect], stringsAsFactors = FALSE), file=paste0(save.path,"/[1]ImputeQC_myGD.csv"), row.names=FALSE, fileEncoding = "UTF-8")
      write.csv(x=myGM[c(1, 1+filter.intersect),], file=paste0(save.path,"/[1]ImputeQC_myGM.csv"), row.names=FALSE, fileEncoding = "UTF-8")
      write.csv(x=myGT.impute[c(1, 1+filter.intersect),], file=paste0(save.path,"/[1]ImputeQC_myGT.csv"), row.names=FALSE, fileEncoding = "UTF-8")
    }
    
 
    
    if(family=="gaussian"){
      if( !(QC|imputation) ){
        myDATA <- list(myX=myX, myY.original=myY.original, myY=myY, myGD=myGD, myGM=myGM, myGT=myGT)
      } else if( imputation & QC ){
        myDATA <- list(myX=myX.impute[c(1,1+filter.intersect),],
                       myY.original=myY.original,
                       myY=myY,
                       myGD=myGD.impute[,filter.intersect],
                       myGM=myGM[c(1,1+filter.intersect),],
                       myGT=myGT.impute[c(1,1+filter.intersect),])
      } else if( QC ){
        myDATA <- list(myX=myX[c(1,1+filter.intersect),],
                          myY.original=myY.original,
                          myY=myY,
                          myGD=myGD[,filter.intersect],
                          myGM=myGM[c(1,1+filter.intersect),],
                          myGT=myGT[c(1,1+filter.intersect),])
      } else if( imputation ){
        myDATA <- list(myX=myX.impute, myY.original=myY.original, myY=myY, myGD=myGD.impute, myGM=myGM, myGT=myGT.impute)
      }
    } else {
      if( (!QC) & (!imputation) ){
        myDATA <- list(myX=myX, myY=myY, myGD=myGD, myGM=myGM, myGT=myGT)  
      } else if( imputation & QC ){
        myDATA <- list(myX=myX.impute[c(1,1+filter.intersect),],
                       myY=myY,
                       myGD=myGD.impute[,filter.intersect],
                       myGM=myGM[c(1,1+filter.intersect),],
                       myGT=myGT.impute[c(1,1+filter.intersect),])
      } else if( QC ){
        myDATA <- list(myX=myX[c(1,1+filter.intersect),],
                          myY=myY,
                          myGD=myGD[,filter.intersect],
                          myGM=myGM[c(1,1+filter.intersect),],
                          myGT=myGT[c(1,1+filter.intersect),])
      } else if( imputation ){
        myDATA <- list(myX = myX.impute,
                              myY = myY,
                              myGD = myGD.impute,
                              myGM = myGM,
                              myGT = myGT.impute
                            )
      }
    }
    
    
    save( myDATA, file=paste0(save.path,"/[1]Data",".RData"))
    
    
    sink(file = paste0(save.path,"/[1]import.hapamp_information.txt"))
    print(
      list(
        genotype = deparse(substitute(genotype)),
        phenotype = deparse(substitute(phenotype)),
        input.type = c("object", "path"),
        save.path = save.path,
        y.col = y.col,
        y.id.col = y.id.col, 
        normalization = normalization,
        remove.missingY = remove.missingY,
        imputation = imputation,
        QC = QC,
        callrate.range = callrate.range,
        maf.range = maf.range,
        HWE.range = HWE.range,
        heterozygosity.range = heterozygosity.range
      )
    )
    sink()
    # End ---------------------------------------------------------------------
    
    
    invisible( myDATA )
  }


# library(dplyr)
# library(tidyr)
# library(pbapply)
# library(ggplot2)
# library(bestNormalize)
# library(genetics)
# library(HardyWeinberg)
# library(gridExtra)
# source("./R/png.impute.snp.R")
# source("./R/gapit.functions.R")
# 
# genotype_without_NA <- R.utils::loadToEnv("./data/genotype.rda")$genotype
# phenotype <- R.utils::loadToEnv("./data/phenotype.rda")$phenotype
# # genotype_without_NA <- sp.gwas::genotype
# # phenotype <- sp.gwas::phenotype
# genotype_with_NA <- genotype_without_NA
# set.seed(123)
# for( j in 1:100 ){
#   genotype_with_NA[-1,-(1:11)][j, sample(10, 2, replace = FALSE)] <- "NN"
# }
# genotype <- genotype_with_NA
# input.type = c("object", "path")[1]
# imputation = TRUE
# QC = TRUE
# callrate.range = c(0.95, 1)
# maf.range = c(1e-3, 1)
# HWE.range = c(0, 1)
# heterozygosity.range = c(0, 1)
# remove.missingY = TRUE
# save.path = "./EXAMPLE_tmp"
# y.id.col = 1
# y.col = 2:4
# normalization = FALSE
# family="gaussian"

