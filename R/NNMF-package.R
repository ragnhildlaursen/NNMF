#' NNMF: Neighborhood Non-negative Matrix Factorization
#'
#' A tool to analyze spatial single cell to recover soft clusters together with corresponding gene programs. 
#' The package provides a fast function 'nnmf' to calculate the factors of non-negative matrix factorization (NMF) for spatial data. 
#' It makes neighboring observations share similar weights for the different factors using gaussian smoothing 
#' in the multiplicative updates of NMF.
#'
#' @details
#' For a guide and overview of the package, please see the
#' vignette:
#' \code{vignette("NNMFguide", package="NNMF")}
#'
#' @name NNMF
#'
#' @import Rcpp
#' @import RcppArmadillo
#' 
#' @keywords internal
#'
"_PACKAGE"