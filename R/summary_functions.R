#' Combine estimates from different pipelines for comparison
#'
#' @param truth tibble. Return value of `estimate_truth()`
#' @param iptw tibble. Return value of `iptw_pipeline()`
#' @param ice tibble. Return value of `ice_pipeline()`
#' @param or tibble. Return value of `or_pipeline()`
#' @param Tt max periods.
#'
#' @return tibble.
#' @export
#'
#' @examples
combine_estimates <- function(truth, iptw, ice, or, Tt) {

  result = tibble::tibble(t=1:Tt)

  if(!missing(truth)) result = dplyr::left_join(result, dplyr::rename(truth, truth=estimate), by='t')
  if(!missing(iptw)) result = dplyr::left_join(result, dplyr::rename(iptw, iptw=estimate), by='t')
  if(!missing(ice)) result = dplyr::left_join(result, dplyr::rename(ice, ice=estimate), by='t')
  if(!missing(or)) result = dplyr::left_join(result, dplyr::rename(or, or=estimate), by='t')

  return(result)
}


#' Estimate true values of the parameters from the simulated potential outcomes
#'
#' @param df_po data frame. Return value of `generate_data(potential_outcomes=TRUE)`
#' @param Tt max periods.
#'
#' @return tibble.
#' @export
#'
#' @examples
estimate_truth <- function(df_po, Tt, link_fun=NULL) {

  if(is.null(link_fun)) {
    return( tibble::tibble(t=1:Tt, estimate = colMeans(calc_ydiffs(df_po, Tt))) )
  } else {
    stopifnot(is.function(link_fun))
    ymeans = colMeans(df_po[ , sapply(0:Tt, function(t) glue('Y{t}')) ])
    return( tibble::tibble(t=1:Tt, estimate = diff(link_fun(ymeans)) ) )
  }


}


#' Estimate confidence interval coverage for a normality-assuming interval based on the simulation variance.
#'
#' @param results_df Data frame with columns 'rep' and 'estimates' where 'estimates' is a list of tibble output results from either iptw_pipeline(), ice_pipeline(), or or_pipeline().
#' @param estimates_truth Data frame with columns 't' and 'estimate', with the 'true' value of E(Y_t(a)-Y_t-1(a)) for t=1,...,T
#' @param conf_level num. level of confidence, e.g. .95 is a 95% confidence interval
#'
#' @return data frame with coverage estimates
#' @export
#'
#' @examples
estimate_ci_coverage = function(results_df, estimates_truth, conf_level=0.95) {
  z = -qnorm((1 - conf_level)/2)
  results_df %>%
    tidyr::unnest(estimates) %>%
    dplyr::left_join(estimates_truth %>% dplyr::rename(truth=estimate)) %>%
    dplyr::group_by(t) %>%
    dplyr::mutate(lwr_normal = estimate - z*sd(estimate),
           upr_normal = estimate + z*sd(estimate),
           cover = 1*(lwr_normal < truth)*(truth < upr_normal))  %>%
    dplyr::summarise(coverage = mean(cover))
}


#' Estimate the bias and variance of DID g-formula estimators
#'
#' @param results_df Data frame with columns 'rep' and 'estimates' where 'estimates' is a list of tibble output results from either iptw_pipeline(), ice_pipeline(), or or_pipeline().
#' @param estimates_truth Data frame with columns 't' and 'estimate', with the 'true' value of E(Y_t(a)-Y_t-1(a)) for t=1,...,T
#'
#' @return data frame with bias and variance estimates
#' @export
#'
#' @examples
estimate_bias_variance = function(results_df, estimates_truth) {
  results_df %>%
    tidyr::unnest(estimates) %>%
    dplyr::left_join(estimates_truth %>% dplyr::rename(truth = estimate)) %>%
    dplyr::group_by(t) %>%
    dplyr::summarise(bias = mean(estimate) - mean(truth),
              variance = var(estimate))
}


