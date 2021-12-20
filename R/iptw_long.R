#' @importFrom rlang .data
calc_inclusion_indicators_long <- function(df_obs, df_interv, idvar, timevar, exposure) {

  df_obs[['A_interv']] = df_interv[[exposure]]

  df_obs %>%
    group_by(.data[[idvar]]) %>%
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


estimate_iptw_long <- function(df_obs,
                               yvar,
                               ylagvar,
                               ip_weights,
                               inclusion_indicators,
                               timevar,
                               link_fun = NULL,
                               binomial_n = NULL) {

  link_fun = check_link_function(link_fun)

  df_obs[['IPW__']] = ip_weights
  df_obs[['INCL__']] = inclusion_indicators
  df_obs[['BINN__']] = binomial_n

  t_list = split(df_obs, df_obs[[timevar]])

  estimate_for_one_timepoint = function(incl, y, ipw, binn) sum(incl*y*ipw)/sum(incl*ipw*binn)

  Eyat = sapply(t_list, function(df) estimate_for_one_timepoint(incl = df$INCL__, y = df[[yvar]], ipw = df$IPW__, binn = df$BINN__))
  Eyas = sapply(t_list, function(df) estimate_for_one_timepoint(incl = df$INCL__, y = df[[ylagvar]], ipw = df$IPW__, binn = df$BINN__))

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
    stopifnot(length(binomial_n) == nrow(df_obs))
    df_obs$freq_w = binomial_n / sum(binomial_n)
  } else {
    binomial_n = df_obs$freq_w = 1
  }

  exposure_variable = terms(as.formula(den_formula))[[2]]

  den_mod = withCallingHandlers({
    glm(den_formula, family, df_obs, weights = freq_w)
  }, warning = function(w) {
    #the non-integer successes warning happens anytime you weight a binomial or other discrete glm, and isn't actually a problem
    if (startsWith(conditionMessage(w), "non-integer #successes"))
      invokeRestart("muffleWarning")})

  den_preds = dens(df_interv[[exposure_variable]], model=den_mod, newdata=df_interv)

  if(!missing(num_formula)) {

    num_mod = withCallingHandlers({
      glm(num_formula, family, df_obs, weights = freq_w)
    }, warning = function(w) {
      #the non-integer successes warning happens anytime you weight a binomial or other discrete glm, and isn't actually a problem
      if (startsWith(conditionMessage(w), "non-integer #successes"))
        invokeRestart("muffleWarning")})

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
                                 timevar = timevar,
                                 link_fun = pt_link_fun,
                                 binomial_n = binomial_n)
  if(!tibble) {
    result = estimates
  } else {
    result = tibble::tibble(t = 1:length(estimates),
                            estimate = estimates)
  }

  if(models) attr(result, 'models') = list(den=den_mod, num=num_mod)

  return(result)
}



