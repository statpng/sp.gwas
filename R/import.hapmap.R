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
import.hapmap <- function(genotype.path=NULL, phenotype.path=NULL, input.type=c("path", "object"), save.path, y.col=NULL, y.id.col=2, family="gaussian"){

  # Import a phenotype data -------------------------------------------------
  if( length( dim(phenotype.path) ) > 0 | input.type=="object" ){
    myY.init <- phenotype.path
  } else {
    if(is.null(phenotype.path)) {
      print("Select the path of phenotype data which has sample IDs in its first column.")
      phenotype.path <- choose.files()
    }
    
    if( grepl("*.xlsx", phenotype.path) ){
      myY.init <- read_excel(phenotype.path)
    } else if( grepl("*.csv|*.txt", phenotype.path) ){
      myY.init <- as.data.frame( fread(phenotype.path), stringsAsFactors=FALSE )
    } else {
      myY.init <- as.data.frame( fread(phenotype.path), stringsAsFactors=FALSE )
    }
  }
  
  
  # Restriction of the number of phenotypes ---------------------------------
  if( ncol(myY.init) > 5 & is.null(y.col) ){
    stop("There are too many phenotypes. Choose 4 or less.")
  } else if(ncol(myY.init) <= 5 & is.null(y.col)){
    myY.init <- myY.init[ , c(y.id.col, which( sapply( myY.init, is.numeric ) ) )]
  } else if (!is.null(y.col)){
    myY.init <- myY.init[,c(y.id.col, y.col)]
  }
  # myY.init <- myY.init[!duplicated(myY.init[,y.id.col]),]
  
  if( any(duplicated(myY.init[,y.id.col])) ){
    cat( "Duplicated phenotype IDs: {", myY.init[,y.id.col][ duplicated(myY.init[,y.id.col]) ], "}\n" )
    stop("There are duplicated sample IDs in phenotype data.")
  }
  
  
# Import a genotype data --------------------------------------------------
    if( length( dim(genotype.path) ) > 0 | input.type=="object" ){
      myX.init <- genotype.path
    } else {
      if(is.null(genotype.path)) {
          print("Select genotype data with hapmap format.")
          genotype.path <- choose.files()
      }
  
      if( grepl("*.xlsx", genotype.path) ){
          myX.init <- read_excel(genotype.path, col_names=FALSE)
      } else if( grepl("*.csv|*.txt", genotype.path) ){
          myX.init <- as.data.frame( fread(genotype.path, header=FALSE), stringsAsFactors=FALSE )
      } else {
        myX.init <- as.data.frame( fread(genotype.path, header=FALSE), stringsAsFactors=FALSE )
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


# Remove missing values ---------------------------------------------------
    print(paste0("Removing the missing values of phenotypes(",
                 sum(apply(myY.init, 1, function(x) any(is.na(x)))),
                 ")."
                 ))
    
    myY.init <- myY.init[!apply(myY.init, 1, function(x) any(is.na(x))), ]
    
# Match the genotype with the phenotype -----------------------------------
    print("Matcing procedure between genotpe and phenotype")
    taxa.snp <- as.character( myX.init[1,12:ncol(myX.init), drop=TRUE] )
    taxa.phenotype <- myY.init[, 1]
    
    cat("IDs for genotype data are \n{", head(taxa.snp, 10), "}\n")
    cat("IDs for phenotype data are \n{", head(taxa.phenotype, 10), "}\n")
    
    ss.x <- which(taxa.snp %in% intersect( taxa.snp, taxa.phenotype ) )
    ss.y <- which(taxa.phenotype %in% intersect( taxa.snp, taxa.phenotype ) )
    
    if( length(ss.x) < length(taxa.snp)/2 ) warning("More half of samples does not match")
    if( length(ss.y) < length(taxa.phenotype)/2 ) warning("More half of samples does not match")
    
    if( length(ss.x) < 10 ) stop("Check your sample IDs")
    if( length(ss.y) < 10 ) stop("Check your sample IDs")
    

# Numericalize the coding of genotypes ------------------------------------
    print("Performing the numericalization procedure for genotpe data.")
    myGD.init <- suppressWarnings( as.data.frame( apply(myX.init[-1,-(1:11)], 1, function(one) GAPIT.Numericalization(one, bit=2, impute="Middle", Major.allele.zero=TRUE)), stringsAsFactors = FALSE ) )
    myGM.init <- myX.init[,1:4]
    myGT.init <- myX.init[,c(12:ncol(myX.init))]
    colnames(myGD.init) <- myX.init[-1,1]
    rownames(myGD.init) <- myX.init[1,-c(1:11)]

    missing.snp <- apply( myX.init[-1,-(1:11)], 1, function(X) sum(is.na(X)) )
    missing.sample <- apply( myX.init[-1,-(1:11)], 2, function(X) sum(is.na(X)) )

    write_xlsx(x=data.frame(ID=myX.init[-1,1], missing.snp, stringsAsFactors = FALSE), path=paste0(save.path,"/[0]missing.snp.xlsx"))
    write_xlsx(x=data.frame(ID=as.character(myX.init[1,-(1:11)]), missing.sample, stringsAsFactors = FALSE), path=paste0(save.path,"/[0]missing.sample.xlsx"))

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
    print("Performing the normalization for phenotypes, but there is no evidence that they can be normalized.")

      
      myY.original <- myY
    for( j in 1:(ncol(myY)-1) ){
      myY[,j+1] <- boxcox(myY.original[,j+1], standardize = FALSE)$x.t
    }

    norm.pvalue.original <- norm.pvalue <- NULL
    for( j in 1:(ncol(myY.original)-1) ){
        norm.pvalue.original[j] <- round( shapiro.test(myY.original[,j+1])$p.value, 4)
        norm.pvalue[j] <- round( shapiro.test(myY[,j+1])$p.value, 4)
    }

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
        xlab("Transformed phenotypes") +
        ggtitle(paste0("Histogram of transformed phenotypes")) +
        theme_bw() +
        facet_wrap(~PHENOTYPE, nrow=2, scales="free") +
        geom_label(x=Inf, y=Inf, aes(label=pvalue),
                   vjust = "inward", hjust = "inward")

    Hist.norm <- grid.arrange(Hist.y.original, Hist.y, nrow=1, heights=15)

    ggsave(Hist.norm, filename = paste0(save.path,"/[1]Hist.phenotype.norm.pdf"), width = 10, height = 5)

    }
    
    write.csv(x=myX, file=paste0(save.path,"/[1]myX.csv"), row.names=FALSE, fileEncoding = "UTF-8")
    if(family=="gaussian"){
      write.csv(x=myY.original, file=paste0(save.path,"/[1]myY.original.csv"), row.names=FALSE, fileEncoding = "UTF-8")
    }
    write.csv(x=myY, file=paste0(save.path,"/[1]myY.csv"), row.names=FALSE, fileEncoding = "UTF-8")
    write.csv(x=data.frame(ID=rownames(myGD), myGD, stringsAsFactors = FALSE), file=paste0(save.path,"/[1]myGD.csv"), row.names=FALSE, fileEncoding = "UTF-8")
    write.csv(x=myGM, file=paste0(save.path,"/[1]myGM.csv"), row.names=FALSE, fileEncoding = "UTF-8")
    write.csv(x=myGT, file=paste0(save.path,"/[1]myGT.csv"), row.names=FALSE, fileEncoding = "UTF-8")
    
    if(family=="gaussian"){
      myDATA <- list(myX=myX, myY.original=myY.original, myY=myY, myGD=myGD, myGM=myGM, myGT=myGT)  
    } else {
      myDATA <- list(myX=myX, myY=myY, myGD=myGD, myGM=myGM, myGT=myGT)  
    }
    
    
    save( myDATA, file=paste0(save.path,"/[1]Data",".RData"))

    # End ---------------------------------------------------------------------

    return( myDATA )
}
