#' @importFrom rlang .data
calc_inclusion_indicators_long <- function(df_obs, df_interv, idvar='uid', timevar = 't', exposure='A') {

  df_obs %>%
    dplyr::left_join(df_interv %>% dplyr::select(!!idvar, !!timevar, A_interv = !!exposure), by=c('uid','t')) %>%
    dplyr::group_by( .data[[idvar]] ) %>%
    dplyr::mutate(incl = cumprod(.data[[exposure]] == .data[['A_interv']])) %>% #this works regardless of whether A is monotonic.
    dplyr::pull(incl)

}
#' @importFrom rlang .data
calc_weights_long <- function(df_obs, den_preds, num_preds, idvar='uid') {
  #assume data are sorted by time in asending order
  df_obs %>%
    dplyr::mutate(den=den_preds, num=num_preds) %>%
    dplyr::group_by(.data[[idvar]]) %>%
    dplyr::mutate(w = cumprod(.data$num/ .data$den)) %>%
    dplyr::pull(w)
}


estimate_iptw_long <- function(df_obs, yvar, ylagvar, ip_weights, inclusion_indicators, freq_weights, timevar, link_fun = NULL) {

  link_fun = check_link_function(link_fun)

  df_obs[['ipw']] = ip_weights
  df_obs[['incl']] = inclusion_indicators
  df_obs[['fw']] = freq_weights

  t_list = split(df_obs, df_obs[[timevar]])

  estimate_for_one_timepoint = function(incl, y, ipw, fw) sum(incl*y*ipw*fw)/sum(incl*ipw*fw)

  Eyat = sapply(t_list, function(df) estimate_for_one_timepoint(incl = df$incl, y = df[[yvar]], ipw = df$ipw, fw = df$fw))
  Eyas = sapply(t_list, function(df) estimate_for_one_timepoint(incl = df$incl, y = df[[ylagvar]], ipw = df$ipw, fw = df$fw))

  link_fun(Eyat) - link_fun(Eyas)
}

iptw_pipeline_long <- function(df_obs,
                               df_interv,
                               den_formula,
                               num_formula,
                               family,
                               yvar='Y',
                               ylagvar = 'Y_lag',
                               idvar = 'uid',
                               timevar = 't',
                               tibble=TRUE,
                               pt_link_fun=NULL,
                               binomial_n=NULL,
                               models=TRUE) {

  if(!is.null(binomial_n)) {
    stopifnot(length(binomial_n) != nrow(data))
    df_obs$freq_w = binomial_n / sum(binomial_n)
  } else {
    df_obs$freq_w = 1
  }

  exposure_variable = terms(as.formula(den_formula))[[2]]

  den_mod = glm(den_formula, family, df_obs, weights = freq_w)
  den_preds = dens(df_interv[[exposure_variable]], model=den_mod, newdata=df_interv)

  if(!missing(num_formula)) {
    num_mod = glm(num_formula, family, df_obs, weights = freq_w)
    num_preds = dens(df_interv[[exposure_variable]], model=num_mod, newdata=df_interv)
  } else {
    num_mod = NULL
    num_preds = 1
  }

  weights = calc_weights_long(df_obs=df_obs, den_preds=den_preds, num_preds=num_preds, idvar=idvar)

  inclusion_indicators = calc_inclusion_indicators_long(df_obs,
                                                        df_interv,
                                                        idvar = idvar,
                                                        timevar=timevar,
                                                        exposure=exposure_variable)

  estimates = estimate_iptw_long(df_obs = df_obs,
                                 yvar = yvar,
                                 ylagvar = ylagvar,
                                 ip_weights = weights,
                                 inclusion_indicators = inclusion_indicators,
                                 freq_weights = df_obs$freq_w,
                                 timevar = timevar,
                                 link_fun = pt_link_fun)
  if(!tibble) {
    result = estimates
  } else {
    result = tibble::tibble(t = 1:length(estimates),
                            estimate = estimates)
  }

  if(models) attr(result, 'models') = list(den=den_mod, num=num_mod)

  return(result)
}



