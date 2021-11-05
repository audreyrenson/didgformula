#' Get histograms of DID g-formula results
#'
#' @param results_df Data frame with columns 'rep' and 'estimates' where 'estimates' is a list of tibble output results from either iptw_pipeline(), ice_pipeline(), or or_pipeline().
#' @param estimates_truth Data frame with columns 't' and 'estimate', with the 'true' value of E(Y_t(a)-Y_t-1(a)) for t=1,...,T
#' @param show_bias Logical. Should the histograms show deviations from 'truth' (TRUE) or the estimates themselves (FALSE)?
#'
#' @return
#' @export
#'
#' @examples
results_histograms = function(results_df, estimates_truth, show_bias=TRUE) {
  #results_df should have a column t and a column estimates, with multiple reps of t=1,2,...,Tt
  results_df %>%
    unnest(estimates) %>%
    left_join(estimates_truth %>% rename(truth = estimate)) %>%
    mutate(bias = estimate - truth) %>%
    ggplot(aes(x=if(show_bias) bias else estimate)) +
    geom_histogram() +
    geom_vline(data=estimates_truth, mapping=aes(xintercept= if(show_bias) 0 else estimate), color='red') +
    facet_wrap(~t, scales='free') +
    labs(x=if(show_bias) 'bias' else 'estimate')

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
    unnest(estimates) %>%
    left_join(estimates_truth %>% rename(truth=estimate)) %>%
    group_by(t) %>%
    mutate(lwr_normal = estimate - z*sd(estimate),
           upr_normal = estimate + z*sd(estimate),
           cover = 1*(lwr_normal < truth)*(truth < upr_normal))  %>%
    summarise(coverage = mean(cover))
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
    unnest(estimates) %>%
    left_join(estimates_truth %>% rename(truth = estimate)) %>%
    group_by(t) %>%
    summarise(bias = mean(estimate) - mean(truth),
              variance = var(estimate))
}
