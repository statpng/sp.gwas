#' sp.gwas
#'
#' @name sp.gwas-package
#' @docType package
#' @title Selection probabilities using generalized linear model with regularization for a SNP data in the hapmap format.
#' @details The penalty function of \code{elastic-net} is defined as \deqn{\alpha||\beta||_1+(1-\alpha)||\beta||_2/2,} where \eqn{\alpha} is a mixing proportion of ridge and the lasso, and \eqn{\beta} is regression coefficients. This penalty is equivalent to the Lasso penalty if \code{alpha=1}. \cr
#' @return A list of data files(genotype, phenotype, etc.), results for selection probabilities, and manhattan plot for multiple traits.

#' @import glmnet CMplot stats ggplot2 utils dplyr 
#' @importFrom gtools mixedorder
#' @importFrom tidyr spread gather unite
#' @importFrom gridExtra grid.arrange
#' @importFrom compiler cmpfun
#' @importFrom bestNormalize boxcox
#' @importFrom writexl write_xlsx
#' @importFrom readxl read_xlsx
#' @importFrom data.table fread

#' @references Zou, H., & Hastie, T. (2005). Regularization and variable selection via the elastic net. Journal of the royal statistical society: series B (statistical methodology), 67(2), 301-320.

#' @keywords hapmap gwas regularization "selection probabilities"


NULL
