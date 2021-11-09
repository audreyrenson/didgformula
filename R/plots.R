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


#' View histograms of deviations from parallel trends in simulated data
#'
#' @param df_po Dataset generated using generate_data(potential_outcomes=TRUE)
#' @param Tt int. max period.
#' @param k int. period index in conditioning event. If null, all periods are show.
#' @param link function. Corresponds to the scale of parallel trends. Typically either I (identity-scale) or qlogis (logit-scale)
#'
#' @return
#' @export
#'
#' @examples
pt_deviation_histograms <- function(df_po, Tt, k=NULL, link_fun=NULL) {

  if(is.null(k)) {

    all_trend_diffs <- lapply(1:Tt, function(k) trend_diffs_k(df_po, k, Tt, link_fun))

    all_trend_diffs %>% #why is it just a tiny bit off? these really should be precisely zero
      purrr::map(tidyr::pivot_longer, dplyr::starts_with("Ydiff")) %>%
      purrr::map(dplyr::ungroup) %>%
      purrr::map(dplyr::select, -dplyr::starts_with('L')) %>%
      dplyr::bind_rows(.id='k') %>%
      dplyr::mutate(k=as.numeric(k), t=as.numeric(stringr::str_extract(name, '[0-9]+'))) %>%
      ggplot2::ggplot(ggplot2::aes(x=value)) +
      ggplot2::geom_histogram() +
      ggplot2::geom_vline(xintercept=0, color='red') +
      ggplot2::facet_grid(k~t, scales="free")
  } else {
    trend_histograms_k(df_po, k, Tt, link_fun)
  }
}

trend_diffs_k <- function(df_po, k, Tt, link_fun) {
  # Used by pt_deviation_histograms
  #
  #checks that parallel trends holds.
  #Showing estimates of:
  #
  # Ydifft = E[Y_t(0) - Y_{t-1}(0)|A_k=0, \bar A_{k-1}=0, \bar L_k ]
  #       - E[Y_t(0) - Y_{t-1}(0)|A_k=1,  bar A_{k-1}=0, \bar L_k ]
  # for t = k, k+1, ..., Tt.
  #
  #These should look like unbiased estimates of zero (since df_po is a finite sample)

  Lbar_k = sapply(0:k, function(m) glue('L{m}'))
  A_kminus1 = glue('A{k-1}')
  A_k = glue('A{k}')
  trends_t = sapply(k:Tt, function(t) glue('Ydiff{t}'))
  Y_t = sapply(0:Tt, function(t) glue('Y{t}'))
  mean_fun = if(is.null(link_fun)) function(...) mean(...) else function(...) link_fun(mean(...))


  # restrict to \bar A_{k-1} = \bar a_{k-1}
  trend_diffs = dplyr::filter(df_po, df_po[[A_kminus1]]==0)
  #condition on \bar L_k
  trend_diffs = dplyr::group_by(trend_diffs, dplyr::across(c({{Lbar_k}}, {{A_k}})))
  #mean Y_t according to A_k, L_k
  trend_diffs = dplyr::summarise(trend_diffs, dplyr::across({{Y_t}}, mean_fun), .groups='keep')
  #take the difference in (possibly transformed) means
  trend_diffs = add_ydiffs(trend_diffs, Tt)
  #ungroup by A_k
  trend_diffs = dplyr::group_by(trend_diffs, dplyr::across({{Lbar_k}}))
  #and take the differences, which should be zero
  trend_diffs = dplyr::summarise(trend_diffs, dplyr::across({{trends_t}}, diff), .groups='keep')

  return(trend_diffs)
}


trend_histograms_k <- function(df_po, k, Tt, link_fun) {
  trend_diffs = trend_diffs_k(df_po, k, Tt, link_fun)
  ydiffs_names = grep('Ydiff', names(trend_diffs), value=TRUE)
  trend_diffs %>%
    tidyr::pivot_longer( {{ydiffs_names}} , names_to="t", names_prefix = "Ydiff") %>%
    dplyr::mutate(t=as.numeric(t)) %>%
    ggplot2::ggplot(ggplot2::aes(x=value)) +
    ggplot2::geom_histogram() +
    ggplot2::geom_vline(xintercept=0, color='red') +
    ggplot2::facet_wrap(~t)
}


add_ydiffs <- function(data, Tt) {
  for(t in 1:Tt) {
    data[[glue('Ydiff{t}')]] = data[[glue('Y{t}')]] - data[[glue('Y{t-1}')]]
  }
  return (data)
}
