#' @name sp.gwas
#' @title Selection probabilities using generalized linear model with regularization for a SNP data in the hapmap format.
#' @description For analysis of high-dimensional genomic data, penalized regression can be a solution to accomodate correlations between predictors. Moreover, selection probabilities do not depend on tuning parameter selection so that it produces a stablity selection. Thresholds are also generated to control the false positives(errors). Thresholds vary with the expected number of false positives to be controlled by the user.
#' Input data is the hapmap formatted SNP data and phenotype data corresponding to SNP data.
#' Output includes files of three types: (1) Matched data files (genotype, numerical, snp info, and phenotype), (2) Results file (selection probabilities and thresholds), (3) Circular manhattan plot with (blue dotted) significant line corresponding to the largest value among user-defined false discoveries.

#' @param genotype.path A path of a snp data which is a \code{p} by \code{(n+11)} matrix of genotypes. It is formatted by hapmap which has (rs, allele, chr, pos) in the first four(1-4) columns, (strand, assembly, center, protLSID, assayLSID, panel, Qcode) in the following seven(5-11) columns. If NULL, user can choose a path in interactive use.
#' @param phenotype.path A path of a phenotype data which is an \code{n} by \code{p} matrix of phenotypes. Since the first some columns can display attributes of the phenotypes, you should enter the arguments, y.col and y.id.col, which represent the columns of phenotypes to be analyzed and the column of sample ID. If NULL, user can choose a path in interactive use.
#' @param input.type If \code{input.type} is "object", obejects of genotype/phenotype.path will be entered. If "path", paths of genotype and phenotype will do.
#' @param y.col The columns of phenotypes. At most 4 phenotypes can be considered, because the plot of them will be fine. Default is 4.
#' @param y.id.col The column of sample ID in the phenotype data file. Default is 2.
#' @param method A method of penalized regression. It includes "lasso" for the lasso and "enet" for the elastic-net.
#' @param family A family of response variable(phenotype). It is "gaussian" for continuous response variable, "binomial" for binary, "poisson" for count, etc. Now you can use only the same family for the multi phenotypes. For more details, see the function(\code{stats::glm}). Default is "gaussian".
#' @param Falsediscovery The expected number of false discovery to be controlled. The larger it is, the higher threshold becomes. Default is c(1, 5, 10).
#' @param plot.ylim A range of the y-axis. If NULL, automatic range in the y-axis will be provided. For plot.ylim=c(0,1), the y-axis has a range of 0 and 1.
#' @param save.path A save.path which has all output files. If there exists save.path, sp.gwas will check if there is an output file. If then, sp.gwas is not going to generate the results but just loading the output files(.RData).
#' @param lambda.min.quantile A range of lambda sequence. Default is 0.5 (median). If the range is so small that it can have many tied selection probabilities which is 1.
#' @param n.lambda The length of lambda sequence. The larger n.lambda, the more detailed lambda sequence will be.
#' @param K	The number of iterations in resampling when calculating the selection probabilities.
#' @param psub The subsampling proportion. For efficiency, default is 0.5.
#' @param manhattan.type A type of manhattan plot to be drawn includes circular('c') and rectangular('m').
#' @param plot.name A name of plot file.
#' @param plot.type A type of plot file which includes "jpg", "pdf", "tiff", etc.
#' @param plot.dpi A resolution of plot. If you want to get a high-resolution image, plot.dpi should be large.

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

#' @author Kipoong Kim <kkp7700@gmail.com>

#' @examples
#' # Not run
#' # sp.gwas(genotype.path = "./input.snp.csv",
#' #         phenotype.path = "./input.phenotype.csv",
#' #         input.type = c("path", "object")[1],
#' #         save.path = "./Test",
#' #         y.col=5:8,
#' #         y.id.col=2,
#' #         method="lasso",
#' #         family="gaussian",
#' #         Falsediscovery = c(1,5,10),
#' #         plot.ylim = NULL,
#' #         lambda.min.quantile = 0.5,
#' #         n.lambda = 10,
#' #         K = 100,
#' #         psub = 0.5,
#' #         manhattan.type = c("c", "r")[2],
#' #         plot.name = "Test",
#' #         plot.type = "jpg",
#' #         plot.dpi = 300)
#' #' A function
#'
#'
#' @return The dataframe with new mean and sum columns
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
#' @importFrom stats quantile
#' @export sp.gwas
sp.gwas <- function( genotype.path = NULL,
                     phenotype.path = NULL,
                     input.type = c("object", "path"),
                     save.path = "./sp.folder",
                     y.col = 4,
                     y.id.col = 2,
                     method = "lasso",
                     family = "gaussian",
                     Falsediscovery = c(1,5,10),
                     plot.ylim = NULL,
                     lambda.min.quantile = 0.5,
                     n.lambda = 10,
                     K = 100,
                     psub = 0.5,
                     manhattan.type = "c",
                     plot.name = "",
                     plot.type = "jpg",
                     plot.dpi = 300
){


    if( !file.exists( paste0(save.path) ) ){
        # system(paste0("mkdir ", save.path))
        dir.create(paste0(save.path))
    }

    if( !file.exists( paste0(save.path,"/[1]Data",".RData") ) ){
        myDATA <- import.hapmap(genotype.path=genotype.path,
                                phenotype.path=phenotype.path,
                                input.type=input.type,
                                save.path=save.path,
                                y.col=y.col,
                                y.id.col=y.id.col,
                                family=family)
    } else {
        print("The data file already exists.")
        load(paste0(save.path,"/[1]Data",".RData"))
    }

    # Beginning of Selection Probabilities ------------------------------------

    if( !file.exists( paste0(save.path,"/[2]sp.res", ".RData") ) ){
        sp.res <- selection.prob(x=myDATA$myGD,
                                 y=myDATA$myY,
                                 save.path=save.path,
                                 snp.info=myDATA$myGM,
                                 method=method,
                                 family=family,
                                 Falsediscovery=Falsediscovery,
                                 lambda.min.quantile=lambda.min.quantile,
                                 n.lambda=n.lambda,
                                 K=K,
                                 psub=psub)
    } else {
        # Load the saved results file ---------------------------------------------
        print("The results file of selection probabilities already exists.")
        load( paste0(save.path,"/[2]sp.res", ".RData") )
    }

    sp.manhattan(sp.df=sp.res$sp.df,
                 threshold=sp.res$threshold,
                 save.path=save.path,
                 manhattan.type=manhattan.type,
                 plot.ylim=plot.ylim,
                 plot.type=plot.type,
                 dpi=plot.dpi
    )

    sink(file = paste0(save.path,"/[4]information.txt"))
    print(
        list(
            genotype.path = genotype.path,
            phenotype.path = phenotype.path,
            save.path = save.path,
            y.col = y.col,
            method = method,
            family = family,
            Falsediscovery = Falsediscovery,
            plot.ylim = plot.ylim,
            lambda.min.quantile = lambda.min.quantile,
            n.lambda = n.lambda,
            K = K,
            psub = psub,
            plot.name = plot.name,
            plot.type = plot.type
        )
    )
    sink()

}


