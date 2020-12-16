#' @name sp.gwas
#' @title Selection probabilities using generalized linear model with regularization for a SNP data in the hapmap format.
#' @description For analysis of high-dimensional genomic data, penalized regression can be a solution to accomodate correlations between predictors. Moreover, selection probabilities do not depend on tuning parameter selection so that it produces a stablity selection. Thresholds are also generated to control the false positives(errors). Thresholds vary with the expected number of false positives to be controlled by the user.
#' Input data is the hapmap formatted SNP data and phenotype data corresponding to SNP data.
#' Output includes files of three types: (1) Matched data files (genotype, numerical, snp info, and phenotype), (2) Results file (selection probabilities and thresholds), (3) Circular manhattan plot with (blue dotted) significant line corresponding to the largest value among user-defined false discoveries.
#' @param input.genotype Either R object or file path can be considered. A genotype data is not a data.frame but a matrix with dimension \code{p} by \code{(n+11)}. It is formatted by hapmap which has (rs, allele, chr, pos) in the first four(1-4) columns, (strand, assembly, center, protLSID, assayLSID, panel, Qcode) in the following seven(5-11) columns. If NULL, user can choose a path in interactive use.
#' @param input.phenotype Either R object or file path can be considered. A phenotype data is an \code{n} by \code{p} matrix. Since the first some columns can display attributes of the phenotypes, you should enter the arguments, y.col and y.id.col, which represent the columns of phenotypes to be analyzed and the column of sample ID. If NULL, user can choose a path in interactive use.
#' @param input.type Default is "object". If \code{input.type} is "object", obejects of genotype/phenotype will be entered, and if "path", paths of genotype/phenotype will be enterd. If you want to use an object, you have to make sure that the class of each column of genotype data is equal to "character".
#' @param input.y.col The columns of phenotypes. At most 4 phenotypes can be considered, because the plot of them will be fine. Default is 2.
#' @param input.y.id.col The column of sample ID in the phenotype data file. Default is 1.
#' @param input.y.norm If TRUE. phenotypes are converted to be normal-shape using box-cox transformation when all phenotypes are positive.
#' @param impute TRUE or FALSE for whether imputation will be conducted.
#' @param impute.type Two imputation methods are supported for (only) imputation=TRUE. Default is "distribution" which impute a genotype from allele distribution. The other is "mode" which indicates an imputation from the most frequent genotype. 
#' @param qc TRUE or FALSE for whether QC for SNPs will be conducted.
#' @param qc.callrate.range A numeric vector indicating the range of non-missing proportion. Default is c(0, 1).
#' @param qc.maf.range A numeric vector indicating the range of minor allele frequency (MAF) to be used. Default is c(0, 1).
#' @param qc.HWE.range A numeric vector indicating the range of pvalue by Hardy-Weinberg Equillibrium to be used. Default is c(0, 1).
#' @param qc.hetero.range A numeric vector indicating the range of heterozygosity values to be used, because, in some cases, heterozygosity higher than expected indicates the low quality variants or sample contamination. Default is c(0, 1).
#' @param qc.remove.missingY If TRUE, the samples with missing values in phenotype data are removed. Accordingly, the corresponding genotype samples are also filtered out. Default is TRUE.
#' @param pop Whether the population structure as a covariate into our model is incorporated or not. The estimation for subpopulation will be conducted using 'gap' statistic and the corresponding PCoA and tSNE plot will be provided.. Default is TRUE.
#' @param pop.tsne Whether the tSNE plot would be displayed or not. Defaults is FALSE.
#' @param pop.clustering The clustering method to classify the samples into subpopulations. Either "kmeans" or "pam" is available. Dafualt is "kmeans".
#' @param pop.gap The method of calculating the gap statistics. Default is "firstSEmax".
#' @param pop.gap.Kmax The maximum number of clusters to be taken into account. Default is 10\% of the sample size.
#' @param pop.gap.nboot The number of bootstraps when calculating the gap statistics. Default is 50. 
#' @param gwas.method A method of penalized regression. It includes "lasso" for the lasso and "enet" for the elastic-net.
#' @param gwas.family A family of response variable(phenotype). It is "gaussian" for continuous response variable, "binomial" for binary, "poisson" for count, etc. Now you can use only the same family for the multi phenotypes. For more details, see the function(\code{stats::glm}). Default is "gaussian".
#' @param gwas.alpha.seq A mixing proportion between lasso and ridge penalties. Default is from 0.1 to 0.9 by 0.2, i.e. elastic-net penalty. When you want to use the lasso penalty, set this value at 1.
#' @param gwas.nlambda The length of lambda sequence. The larger n.lambda, the more detailed lambda sequence will be.
#' @param gwas.lambda.min.quantile A range of lambda sequence. Default is 0.5 (median). If the range is so small that it can have many tied selection probabilities which is 1. To handle with this problem, you should increase the value of "lambda.min.quantile".
#' @param gwas.nrep The number of iterations in resampling when calculating the selection probabilities.
#' @param gwas.psub The subsampling proportion. For efficiency, default is 0.5.
#' @param gwas.seed An integer value to reproduce the analysis result. Default is current date (e.g. 20201111).
#' @param gwas.cpp Whether or not the code using Rcpp will be used.
#' @param threshold.false.discovery The expected number of false discovery to be controlled. The larger it is, the higher threshold becomes. Default is c(1, 5, 10).
#' @param threshold.perm Permutation-based empirical threshold can be provided for \code{false.discovery} if permutation is TRUE. Default is FALSE.
#' @param threshold.nperm The number of permutation replicates. Note that the calculation time increases as \code{nperm} increases. Default is 100.
#' @param plot.display Figures will be displayed or not
#' @param plot.manhattan.type A type of manhattan plot to be drawn includes circular('c') and rectangular('m').
#' @param plot.name A name of plot file.
#' @param plot.ylim A range of the y-axis. If NULL, automatic range in the y-axis will be provided. For plot.ylim=c(0,1), the y-axis has a range of 0 and 1.
#' @param plot.extension A type of plot file which includes "jpg", "pdf", "tiff", etc.
#' @param plot.dpi A resolution of plot. If you want to get a high-resolution image, plot.dpi should be large.
#' @param save.path A save.path which has all output files. If there exists save.path, sp.gwas will check if there is an output file. Note that if there is an output RData file in "save.path", sp.gwas will just load the output files(.RData) in there, thereby not providing the results for new "genotype" and "phenotype".
#' @details The penalty function of \code{elastic-net} is defined as \deqn{\alpha||\beta||_1+(1-\alpha)||\beta||_2/2,} where \eqn{\alpha} is a mixing proportion of ridge and the lasso, and \eqn{\beta} is regression coefficients. This penalty is equivalent to the Lasso penalty if \code{alpha=1}. \cr \cr
#' \code{An algorithm of selection probabilities with elastic-net.} \cr \cr
#' 0 : Let us assume that a genomic data has \code{n} samples and \code{p} variables. \cr
#' 1 : For all \eqn{\Lambda=(\alpha, \lambda)}, where \eqn{\alpha in [0,1]}, \eqn{\lambda>0}. \cr
#' 2 : \code{for} k=1 to K \code{do}. \cr
#' 3 : ------- Subsample \eqn{I_k} with size \eqn{[n/2].} \cr
#' 4 : ------- Compute \eqn{\hat{\beta}_j^{\Lambda}(I_k)} with regularization model. \cr
#' 5 : \code{end for.} \cr
#' 6 : \eqn{SP_j^\Lambda = \frac{1}{K}\#\{k\le K: \hat{\beta}_j^\Lambda(I_k) \ne 0 \}. } \cr
#' 7 : \eqn{SP_j =  \max_\Lambda SP_j^\Lambda, ~~j=1,\cdots, p.} \cr
#' 8 : \code{return} \code{SP}=\eqn{(SP_1, \cdots, SP_p).} \cr
#' @return
#'     \item{Histogram of original and transformed phenotypes}{Histogram of phenotypes with p-value by Shapiro-test on the top right corner.}
#'     \item{myDATA}{A list of myX, myGD, myGM, myGT, myY, and myY.original(for "gaussian").}
#'     \item{sp.res}{A list of sp.df and threshold.}
#'     \item{Circular Manhattan plot}{Manhattan plot for the first phenotype is the innermost circle. Colors for chromosome is fixed, so that if you want to change colors, you would edit the R code of sp.manhattan function.}
#'     
#'     Note that, before your code, you have to specify the setseed value to get reproducible results, because it uses the resampling approach when calculating the selection probability.
#' 
#' @author Kipoong Kim <kkp7700@gmail.com>
#' @references 
#' Zou, H., & Hastie, T. (2005). Regularization and variable selection via the elastic net. Journal of the royal statistical society: series B (statistical methodology), 67(2), 301-320.
#' Meinshausen, N., & Bühlmann, P. (2010). Stability selection. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 72(4), 417-473.
#' Kim, K., Koo, J., & Sun, H. (2020). An empirical threshold of selection probability for analysis of high-dimensional correlated data. Journal of Statistical Computation and Simulation, 1-12.
#' 
#' @examples
#' if(FALSE){
#' 
#' # Note that you have to change the parameters used in this example when applying this code into your works.
#' genotype <- sp.gwas::genotype # load("genotype.rda")
#' phenotype <- sp.gwas::phenotype # load("phenotype.rda")
#' 
#' # elastic-net model
#' sp.gwas(input.genotype = genotype,
#'         input.phenotype = phenotype,
#'         input.type = "object",
#'         input.y.col = 2:3,
#'         input.y.id.col = 1,
#'         input.y.norm = FALSE, # due to the negative values of phenotypes
#' 
#'         impute = TRUE,
#'         impute.type = "distribution", # based on the allele frequency distribution
#' 
#'         qc = TRUE,
#'         qc.callrate.range = c(0.9, 1),
#'         qc.maf.range = c(0.01, 1),
#'         qc.HWE.range = c(0, 1),
#'         qc.hetero.range = c(0, 1),
#'         qc.remove.missingY = TRUE,
#' 
#'         pop = TRUE,
#'         pop.tsne = FALSE,
#'         pop.clustering = "kmeans",
#'         pop.gap = "firstSEmax",
#'         pop.gap.Kmax = NULL,
#'         pop.gap.nboot = 50,
#' 
#'         gwas.method = "enet",
#'         gwas.family = "gaussian",
#'         gwas.alpha.seq = 0.1, # for a fixed alpha value
#'         gwas.nlambda = 10,
#'         gwas.lambda.min.quantile = 0.5,
#'         gwas.nrep = 10, # the number of resampling replicates
#'         gwas.psub = 0.5,
#'         gwas.seed = NULL,
#'         gwas.cpp = FALSE,
#' 
#'         threshold.false.discovery = c(1,5,10),
#'         threshold.perm = TRUE,
#'         threshold.nperm = 10,
#' 
#'         plot = FALSE,
#'         plot.manhattan.type = "c",
#'         plot.name = "",
#'         plot.ylim = NULL,
#'         plot.extension = "jpg",
#'         plot.dpi = 300,
#' 
#'         save.path = "./EXAMPLE_enet")
#'         
#' # Manhattan plot for the first phenotype (Y1)
#' results <- read.csv("./EXAMPLE_enet/[2]sp.results.csv")
#' class(results$chr) <- "numeric"
#' class(results$pos) <- "numeric"
#' thresholds <- read.csv("./EXAMPLE_enet/[2]sp.thresholds.csv")
#' threshold_Y1_FD1 <- subset( subset(thresholds, Method=="Theoretical"), FD==1)$Y1
#' threshold_Y1_FD10 <- subset( subset(thresholds, Method=="Theoretical"), FD==10)$Y1
#' highlight1 <- results$rs[results$Y1 > threshold_Y1_FD1]
#' highlight10 <- setdiff( results$rs[results$Y1 > threshold_Y1_FD10], highlight1 )
#' 
#' jpeg("./EXAMPLE_enet/Manhattan_from_qqman.jpeg", width=12, height=5, unit="in", res=600)
#' qqman_manhattan(results,
#'                 chr="chr", bp="pos", snp="rs", p="Y1", logp=FALSE, 
#'                 suggestiveline = threshold_Y1_FD1,
#'                 genomewideline = threshold_Y1_FD10,
#'                 highlight = list(highlight1, highlight10),
#'                 col.highlight = c("blue", "red"),
#'                 ylab="Selection probabilities", ylim=c(0,1))
#' dev.off()
#' 
#' 
#' # Lasso model with permuted threshold for multinomial phenotypes
#' 
#' genotype <- sp.gwas::genotype # load("genotype.rda")
#' phenotype <- sp.gwas::phenotype # load("phenotype.rda")
#' 
#' phenotype[,2] <- gtools::quantcut(phenotype[,2], 3)
#' phenotype[,3] <- gtools::quantcut(phenotype[,3], 3)
#' phenotype[,4] <- gtools::quantcut(phenotype[,4], 3)
#' 
#' 
#' sp.gwas(input.genotype = genotype,
#'         input.phenotype = phenotype,
#'         input.type = "object",
#'         input.y.col = 2:4,
#'         input.y.id.col = 1,
#'         input.y.norm = FALSE, # due to the negative values of phenotypes
#'         
#'         impute = TRUE,
#'         impute.type = "distribution", # based on the allele frequency distribution
#'         
#'         qc = TRUE,
#'         qc.callrate.range = c(0.9, 1),
#'         qc.maf.range = c(0.01, 1),
#'         qc.HWE.range = c(0, 1),
#'         qc.hetero.range = c(0, 1),
#'         qc.remove.missingY = TRUE,
#'         
#'         pop = TRUE,
#'         pop.tsne = FALSE,
#'         pop.clustering = "kmeans",
#'         pop.gap = "firstSEmax",
#'         pop.gap.Kmax = NULL,
#'         pop.gap.nboot = 50,
#'         
#'         gwas.method = "lasso",
#'         gwas.family = "multinomial",
#'         gwas.alpha.seq = 1, # for a fixed alpha value
#'         gwas.nlambda = 10,
#'         gwas.lambda.min.quantile = 0.5,
#'         gwas.nrep = 10, # the number of resampling replicates
#'         gwas.psub = 0.5,
#'         gwas.seed = NULL,
#'         gwas.cpp = FALSE,
#'         
#'         threshold.false.discovery = c(1,5,10),
#'         threshold.perm = TRUE,
#'         threshold.nperm = 10,
#'         
#'         plot = FALSE,
#'         plot.manhattan.type = "c",
#'         plot.name = "",
#'         plot.ylim = NULL,
#'         plot.extension = "jpg",
#'         plot.dpi = 300,
#'         
#'         save.path = "./EXAMPLE_lasso_multinomial")
#' 
#' 
#' png.manhattan_from_dir("./EXAMPLE_lasso_multinomial", "Theoretical", FD=c(1,5,10))
#' png.manhattan_from_dir("./EXAMPLE_lasso_multinomial", "Permuted", FD=c(1,5,10))
#' 
#' 
#'
#' }
#'
#'
#' @return The dataframe with new mean and sum columns
#' @import CMplot
#' @importFrom dplyr %>%
#' @import ggplot2
#' @import glmnet
#' @import utils
#' @importFrom genetics HWE.exact
#' @importFrom genetics genotype
#' @importFrom bestNormalize boxcox
#' @importFrom compiler cmpfun
#' @importFrom data.table fread
#' @importFrom gridExtra grid.arrange
#' @importFrom gtools mixedorder
#' @importFrom readxl read_xlsx
#' @importFrom writexl write_xlsx
#' @importFrom tidyr gather
#' @importFrom tidyr spread
#' @importFrom tidyr unite
#' @importFrom stats quantile
#' @export sp.gwas
sp.gwas <- function( input.genotype = NULL,
                     input.phenotype = NULL,
                     input.type = c("object", "path"),
                     input.y.col = 2,
                     input.y.id.col = 1,
                     input.y.norm = TRUE,
                     
                     impute = FALSE,
                     impute.type = c("distribution", "mode"),
                     
                     qc = TRUE,
                     qc.callrate.range = c(0, 1),
                     qc.maf.range = c(0, 1),
                     qc.HWE.range = c(0, 1),
                     qc.hetero.range = c(0, 1),
                     qc.remove.missingY = TRUE,
                     
                     pop = TRUE,
                     pop.tsne = FALSE,
                     pop.clustering = "kmeans",
                     pop.gap = "firstSEmax",
                     pop.gap.Kmax = NULL,
                     pop.gap.nboot = 50,
                     
                     gwas.method = "lasso",
                     gwas.family = "gaussian",
                     gwas.alpha.seq = seq(0.1,0.9,by=0.2),
                     gwas.nlambda = 10,
                     gwas.lambda.min.quantile = 0.5,
                     gwas.nrep = 100,
                     gwas.psub = 0.5,
                     gwas.seed = NULL,
                     gwas.cpp = FALSE,
                     
                     threshold.false.discovery = c(1,5,10),
                     threshold.perm = FALSE,
                     threshold.nperm = 100,
                     
                     plot.display = FALSE,
                     plot.manhattan.type = "c",
                     plot.name = "",
                     plot.ylim = NULL,
                     plot.extension = "jpg",
                     plot.dpi = 300,
                     
                     save.path = "./"
){
    
    if( gwas.method == "lasso" & gwas.alpha.seq!=1 ) stop("For the lasso, gwas.alpha.seq have to be equal to 1")
    if( gwas.method == "enet" & gwas.alpha.seq==1 ) stop("For the elastic-net, gwas.alpha.seq should not be equal to 1")
    
    print("Removing the SNPs with only one value")
    wh.onevalue <- pbapply( input.genotype[-1,-(1:11)], 1, function(xj){
        ( length(unique(xj)) == 1 )
    } )
    
    genotype.onevalue <- input.genotype[ c(1, which(wh.onevalue)+1), ]
    input.genotype <- input.genotype[ c(1, which(!wh.onevalue)+1), ]
    
    cat( "# of SNPs with only one genotype =", sum( wh.onevalue ), "\n" )
    
    
    if( is.null(gwas.seed) ) gwas.seed <- gsub("-", "", Sys.Date())
    
    # system(paste0("mkdir ", save.path))
    
    
    # if( !file.exists( paste0(save.path,"/[1]Data",".RData") ) ){
    
    if( ifelse( is.null(save.path), TRUE, save.path %in% c("./", ".", "/", "" ) ) ){
        save.path <- paste0(save.path, "/method=", gwas.method, "__seed=", gwas.seed, "__time=", gsub("|\\:", "", gsub("([0-9]{2}:[0-9]{2}):[0-9]{2}", "\\1", gsub("\\ ", "-", Sys.time())) ) )
        save.path <- gsub("//", "/", save.path)
    } else {
        if( !file.exists( save.path ) ){
            dir.create(paste0(save.path))
        } else {
            stop("Change your save.path that already exists.")
        }
    }
        
        set.seed(gwas.seed)
        myDATA <- import.hapmap(genotype=input.genotype,
                                phenotype=input.phenotype,
                                input.type=input.type,
                                y.col = input.y.col,
                                y.id.col = input.y.id.col,
                                normalization = input.y.norm,
                                
                                imputation=impute,
                                impute.type=impute.type,
                                
                                QC = qc,
                                maf.range = qc.maf.range,
                                remove.missingY = qc.remove.missingY, 
                                HWE.range = qc.HWE.range, 
                                callrate.range = qc.callrate.range,
                                heterozygosity.range = qc.hetero.range,
                                
                                family=gwas.family,
                                
                                save.path = save.path
                                )
    # }
    
    # else {
    #     print("The dataset file of genotype and phenotype already exists.")
    #     load(paste0(save.path,"/[1]Data",".RData"))
    #     
    #     if( qc.remove.missingY ){
    #         
    #         if( any(apply(myDATA$myY,1,function(x) any(is.na(x)))) ){
    #             print(paste0("Removing the missing values of phenotypes(",
    #                          sum(apply(myData$myY, 1, function(x) any(is.na(x)))),
    #                          ")."
    #             ))
    #             
    #             wh.nonmissing <- which(!apply(myDATA$myY,1,function(x) any(is.na(x))))
    #             id.nonmissing <- myDATA$myY[,1][which(!apply(myDATA$myY,1,function(x) any(is.na(x))))]
    #             
    #             
    #             myDATA$myX <- myDATA$myX[, c(1:11, 11+which(unlist(myDATA$myX[1,-(1:11)])%in%id.nonmissing))]
    #             myDATA$myY <- myDATA$myY[ which(myDATA$myY[,1] %in% id.nonmissing), ]
    #             if( family=="gaussian" ){
    #                 myDATA$myY.original <- myDATA$myY.original[ which(myDATA$myY.original[,1] %in% id.nonmissing), ]
    #             }
    #             myDATA$myGD=myDATA$myGD[ which(rownames(myDATA$myGD) %in% id.nonmissing), ]
    #             myDATA$myGT=myDATA$myGT[ , which(myDATA$myGT[1,] %in% id.nonmissing)]
    #             
    #             save( myDATA, file=paste0(save.path,"/[1]Data_missingY_removed",".RData"))
    #             
    #             # myDATA$myX <- myDATA$myX[, c(1:11, 11+wh.nonmissing)]
    #             # myDATA$myY <- myDATA$myY[ wh.nonmissing, ]
    #             # if( family=="gaussian" ){
    #             #     myDATA$myY.original <- myDATA$myY.original[ wh.nonmissing, ]
    #             # }
    #             # myDATA$myGD=myDATA$myGD[ wh.nonmissing, ]
    #             # myDATA$myGT=myDATA$myGT[ , wh.nonmissing]
    #         }
    #     }
    # }
    
    
    # Estimating population structure -----------------------------------------
    
    # png.population()
    
    if( pop ){
        
        if(is.null(pop.clustering)) pop.clustering <- "kmeans"
        if(is.null(pop.gap)) pop.gap <- "firstSEmax"
        
        set.seed(gwas.seed)
        res.population <- png.population(snp = myDATA$myGD,
                                         cluster.method = pop.clustering,
                                         gap.method = pop.gap,
                                         pca.scale = FALSE,
                                         K.max = pop.gap.Kmax,
                                         nboot = pop.gap.nboot,
                                         tsne = pop.tsne,
                                         save.path = save.path)
    } else {
        res.population <- list(group=NULL)
    }
    
    # Beginning of Selection Probabilities ------------------------------------
    
    # x=myDATA$myGD;
    # y=myDATA$myY;
    # pop=res.population$group;
    # save.path=save.path;
    # snp.info=myDATA$myGM;
    # method=method;
    # family=family;
    # false.discovery=false.discovery;
    # permutation=permutation;
    # nperm=nperm;
    # lambda.min.quantile=lambda.min.quantile;
    # alpha.seq=alpha.seq;
    # n.lambda=n.lambda;
    # K=K;
    # psub=psub
    # ... <- NULL

    # if( !file.exists( paste0(save.path,"/[2]sp.res", ".RData") ) ){


    
    sp.res <- selection.prob(x=myDATA$myGD,
                             y=myDATA$myY,
                             pop=res.population$group,
                             save.path=save.path,
                             snp.info=myDATA$myGM,
                             
                             method=gwas.method,
                             family=gwas.family,
                             lambda.min.quantile=gwas.lambda.min.quantile,
                             alpha.seq=gwas.alpha.seq,
                             n.lambda=gwas.nlambda,
                             K=gwas.nrep,
                             psub=gwas.psub,
                             seed=gwas.seed,
                             
                             false.discovery=threshold.false.discovery,
                             permutation=threshold.perm,
                             nperm=threshold.nperm
                             )
            
    
    

# Time comparison ---------------------------------------------------------
    
    # microbenchmark::microbenchmark(
    #     selection.prob_cpp(x=myDATA$myGD,
    #                                  y=myDATA$myY,
    #                                  pop=res.population$group,
    #                                  save.path=save.path,
    #                                  snp.info=myDATA$myGM,
    #                                  method=method,
    #                                  family=family,
    #                                  false.discovery=false.discovery,
    #                                  permutation=permutation,
    #                                  nperm=nperm,
    #                                  lambda.min.quantile=lambda.min.quantile,
    #                                  alpha.seq=alpha.seq,
    #                                  n.lambda=n.lambda,
    #                                  K=K,
    #                                  psub=psub,
    #                                  seed=as.double(seed)),
    #     
    #     selection.prob(x=myDATA$myGD,
    #                              y=myDATA$myY,
    #                              pop=res.population$group,
    #                              save.path=save.path,
    #                              snp.info=myDATA$myGM,
    #                              method=method,
    #                              family=family,
    #                              false.discovery=false.discovery,
    #                              permutation=permutation,
    #                              nperm=nperm,
    #                              lambda.min.quantile=lambda.min.quantile,
    #                              alpha.seq=alpha.seq,
    #                              n.lambda=n.lambda,
    #                              K=K,
    #                              psub=psub,
    #                              seed=seed)
    # )
    

# -------------------------------------------------------------------------
    
    
    
    # sp.res$threshold
    # sp.res2$threshold
    
    # } else {
    
    
    #     # Load the saved results file ---------------------------------------------
    #     print("The results file of selection probabilities already exists.")
    #     load( paste0(save.path,"/[2]sp.res", ".RData") )
    # }
    
    
    if( plot.display ){
        sp.manhattan( sp.df=sp.res$sp.df,
                      threshold=sp.res$threshold,
                      save.path=save.path,
                      plot.manhattan.type=plot.manhattan.type,
                      plot.ylim=plot.ylim,
                      plot.extension=plot.extension,
                      dpi=plot.dpi )
        
        plot.manhattan.type <- NULL
        plot.name <- NULL
        plot.extension <- NULL
        plot.dpi <- NULL
    }
    
    
    
    sink(file = paste0(save.path,"/[4]information.txt"))
    print(
        list(
            input.genotype = deparse(substitute(input.genotype)),
            input.phenotype = deparse(substitute(input.phenotype)),
            input.y.col = input.y.col,
            input.y.id.col = input.y.id.col,
            input.y.norm = input.y.norm,
            
            impute = impute,
            impute.type = impute.type,
            
            qc = qc,
            qc.callrate.range = qc.callrate.range,
            qc.maf.range = qc.maf.range,
            qc.HWE.range = qc.HWE.range,
            qc.hetero.range = qc.hetero.range,
            qc.remove.missingY = qc.remove.missingY,
            
            pop = pop,
            pop.tsne = pop.tsne,
            pop.gap.Kmax = pop.gap.Kmax,
            
            gwas.method = gwas.method,
            gwas.family = gwas.family,
            gwas.alpha.seq = gwas.alpha.seq,
            gwas.nlambda = gwas.nlambda,
            gwas.lambda.min.quantile = gwas.lambda.min.quantile,
            gwas.nrep = gwas.nrep,
            gwas.psub = gwas.psub,
            gwas.seed = gwas.seed,
            gwas.cpp = gwas.cpp,
            
            threshold.false.discovery = threshold.false.discovery,
            threshold.perm = threshold.perm,
            threshold.nperm = threshold.nperm,
            
            plot = plot,
            plot.manhattan.type = plot.manhattan.type,
            plot.name = plot.name,
            plot.ylim = plot.ylim,
            plot.extension = plot.extension,
            plot.dpi = plot.dpi,
            
            save.path = save.path
        )
        
    )
    sink()
    
    cat("\n", "The results are saved in ", getwd(), "\n")
    
    out.final <- list( myDATA = myDATA,
                       pop = res.population,
                       sp = sp.res )
    
    save(out.final, file = paste0(save.path,"/[4]output.RData"))
    
    
    return( out.final )
    
}


