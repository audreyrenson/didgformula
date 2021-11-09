

#' Bootstrap standard error estimates for didgformula estimators
#'
#' @param data Wide format data frame with one row per individual, and columns Yt for t = 0,1,...,Tt.
#' @param nboots integer. Number of boostrap iterations
#' @param estimator function. One of `*_pipeline()` functions.
#' @param ... Arguments passed to `estimator`
#' @param tibble logical. Return results as a tibble (TRUE) or vector (FALSE)?
#'
#' @return
#' @export
#'
#' @examples
bootstrap_se <- function(data, nboots, estimator, ..., tibble=TRUE) {

  one_bootstrap <- function(data)   return (data[ sample(1:nrow(data), size=nrow(data), replace=TRUE), ])
  many_bootstraps <- function(data, nboots) replicate(nboots, one_bootstrap(data), simplify=FALSE)

  boots = many_bootstraps(data, nboots)
  estimates = lapply(boots, estimator, ..., tibble=FALSE)
  boot_se =  matrixStats::colSds( Reduce(rbind, estimates)  )

    if(tibble) {
    return (tibble::tibble(t=1:length(boot_se), se=boot_se) )
  } else {
    return (boot_se)
  }
}
