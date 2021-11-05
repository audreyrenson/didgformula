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
