# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @name gibbs
#' @title gibbs sample by Rcpp
#' @description Generate the gibbs sample by Rcpp,using the method of gibbs sampling.
#' @param N the sample size
#' @return a Numeric Matrix,sp stands for random sample of multivariate distribution.
#' @examples
#' \dontrun{
#' res <- gibbs_chain(N)
#' }
#' @export
gibbs_chain <- function(N) {
    .Call('_StatComp_gibbs_chain', PACKAGE = 'StatComp', N)
}

