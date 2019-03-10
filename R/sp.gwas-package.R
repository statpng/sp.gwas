#' sp.gwas
#'
#' @name sp.gwas-package
#' @docType package
#' @title Selection probabilities using generalized linear model with regularization for a SNP data in the hapmap format.
#' @details The penalty function of \code{elastic-net} is defined as \deqn{\alpha||\beta||_1+(1-\alpha)||\beta||_2/2,} where \eqn{\alpha} is a mixing proportion of ridge and the lasso, and \eqn{\beta} is regression coefficients. This penalty is equivalent to the Lasso penalty if \code{alpha=1}. \cr
#' @return A list of data files(genotype, phenotype, etc.), results for selection probabilities, and manhattan plot for multiple traits.

#' @import glmnet CMplot stats ggplot2 utils
#' @importFrom tidyr gather
#' @importFrom gridExtra grid.arrange
#' @importFrom compiler cmpfun
#' @importFrom bestNormalize boxcox
#' @importFrom writexl write_xlsx
#' @importFrom readxl read_xlsx
#' @importFrom data.table fread
#' @importFrom dplyr "%>%"
#' @importFrom dplyr select

#' @references Zou, H., & Hastie, T. (2005). Regularization and variable selection via the elastic net. Journal of the royal statistical society: series B (statistical methodology), 67(2), 301-320.

#' @keywords hapmap gwas regularization "selection probabilities"

#' @examples
#' # Not run
#' # sp.gwas(genotype.path = "./input.snp.csv",
#' #         phenotype.path = "./input.phenotype.csv",
#' #         y.col=5:8,
#' #         y.id.col=2,
#' #         method="lasso",
#' #         family="gaussian",
#' #         Falsediscovery = c(1,5,10),
#' #         plot.ylim = NULL,
#' #         save.path = "./Test2",
#' #         lambda.min.quantile = 0.5,
#' #         n.lambda = 10,
#' #         K = 100,
#' #         psub = 0.5,
#' #         plot.name = "Test2",
#' #         plot.type = "jpg",
#' #         plot.dpi = 300)
NULL
