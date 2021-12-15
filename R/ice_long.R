#' Iterated conditional DID g-formula estimator
#'
#' @param data Wide format data frame with one row per individual, and columns Yt for t = 0,1,...,Tt.
#' @param inside_formula_t chr, right-hand-side formula for inside model for Yt
#' @param inside_formula_tmin1 chr, right-hand-side formula for inside model for Yt-1
#' @param outside_formula chr, right-hand-side formula for outside models
#' @param Tt int. max periods
#' @param inside_family stats::family object or string referring to one, as in `glm`.
#' @param binomial_n int length nrow(data). Group sizes for binomial aggregate data.
#' @param pt_link_fun function. The scale on which parallel trends is assumed (e.g., `qlogis` for logit scale). Default `NULL` for untransformed scale.
#' @param tibble logical. return results as a tibble (TRUE) or vector (FALSE)?
#'
#' @return
#' @export
#'
#' @examples
ice_pipeline_long <- function(df_obs,
                              df_interv,
                              inside_formula_t,
                              inside_formula_tmin1,
                              outside_formula,
                              Tt,
                              t_col,
                              inside_family='gaussian',
                              pt_link_fun=NULL,
                              binomial_n=NULL,
                              tibble=TRUE) {

  ice_preds = recursive_ice_long(Tt=Tt,
                                 n_nested=Tt,
                                 df_obs = df_obs,
                                 df_interv=df_interv,
                                 inside_formula_t = inside_formula_t,
                                 inside_formula_tmin1 = inside_formula_tmin1,
                                 outside_formula = outside_formula,
                                 inside_family=inside_family,
                                 binomial_n = binomial_n)

  ice_estimates = estimate_ice_long(ice_preds, t_col, pt_link_fun, binomial_n)
  ice_estimates = sapply(ice_preds, estimate_ice, link_fun=pt_link_fun, binomial_n)

  if(!tibble) {
    return (ice_estimates)
  } else {
    return ( tibble::tibble(t=1:Tt,
                            estimate = ice_estimates))
  }

}

#'  Recursively apply the ICE algorithm for a long dataset
#'
#' @param Tt int. Max time period (total of TT+1) periods
#' @param n_nested int. Index of the current nested expectation being estimated. 0=innermost (non-nested) expectation, Tt=outermost. Generally always Tt when called by the user.
#' @param df_obs long data frame of the observed data with columns uid, t (0,1,...,Tt), Y, Y_lag, exposure and covariates.
#' @param df_interv long data frame with all columns in df_obs, but exposure variables set to their intervened values.
#' @param inside_formula_t chr. full formula for innermost expectation model for Y_t. Will be `glue()`ed internally, so can include {tvars}. (See details)
#' @param inside_formula_tmin1 chr. full formula for innermost expectation model for Y_t-1. Will be `glue()`ed internally, so can include {tvars}. (See details)
#' @param outside_formula chr. RHS-only formula for outer (nested) expectation models for n_nested=0,1,...,Tt. Will be `glue()`ed internally, so can include {tvars} and {n}. (See details)
#' @param inside_family stats::families object or character naming one, for inside model fit using `glm`
#' @param binomial_n int length nrow(df_obs). Group sizes for aggregate binomial data.
#' @param models lgl. Return all models as an attribute?
#'
#' @details Formulas in `recursive_ice_long()` are glue-style and admit 2 special references.
#'
#'         {n} is shorthand for `n_nested`. For example, we typically want exposure and covariate values
#'         lagged exactly `n_nested` periods in each outside model, so we might pass `'~A_lag{n}*L_lag{n}'`
#'         as `outside_formula`. This is simply a more concise version of `'~A_lag{n_nested}*L_lag{n_nested}`.
#'
#'         {tvars} is an internally generated string variable that looks like `"(t{n} + t{n+1} + ... + t{Tt})"`,
#'         where again {n} is shorthand for `{n_nested}`. Specifically, `recursive_ice_long()` retrieves the `timevars`
#'         attribute of `df_obs`, takes the indices `n_nested:length(tvars)`, and pastes these together with `collapse='+'`.
#'         This is useful for including time indicators in saturated models, since the `n`th outside model will only have
#'         observations with `t>=n` in it.
#'
#' @return
#' @export
#'
#' @examples
recursive_ice_long <- function(Tt,
                               n_nested,
                               df_obs,
                               df_interv,
                               inside_formula_t,
                               inside_formula_tmin1,
                               outside_formula,
                               inside_family,
                               binomial_n=NULL,
                               models=TRUE) {

  tvars = attr(df_obs, 'timevars')
  tvars = paste0('(', paste( tvars[max(1, n_nested) : length(tvars)] , collapse='+'), ')') #this is like (t{n} + t{n+1} + ... + t{Tt}) for referencing formulas

  if(n_nested == 0) {
    two_models = fit_two_outcome_models_long(data=df_obs,
                                             formula_t = glue(inside_formula_t),
                                             formula_tmin1 = glue(inside_formula_tmin1),
                                             family=inside_family)
    predictions = sapply(two_models, predict, newdata=df_interv, type='response', simplify=TRUE)

    if(models) attr(predictions, 'models') = list(two_models)

    return (  predictions )
  } else {
    preds = recursive_ice_long(Tt,
                               n_nested - 1,
                               df_obs,
                               df_interv,
                               inside_formula_t,
                               inside_formula_tmin1,
                               outside_formula,
                               inside_family,
                               binomial_n,
                               models)

    obs_keep = df_obs$t >= n_nested #we can only go outward n_nested expectations if that doesn't take us before time 0
    n = n_nested #for convenience, allow user to specify formulas with x{n} where {n} is shorthand for {n_nested}
    formula_n = glue(outside_formula)


    two_models = list()
    two_models[['tmin1']] = fit_outer_exp_model_long(n_nested, df_obs, preds[,1], formula_n, binomial_n, obs_keep)
    two_models[['t']] = fit_outer_exp_model_long(n_nested, df_obs, preds[,2], formula_n, binomial_n, obs_keep)

    predictions = preds #carry outward for observations where t < n_nested
    predictions[obs_keep, ] = sapply(two_models, predict, newdata=df_interv[obs_keep, ], simplify=TRUE)

    if(models) {
      attr(predictions, 'models') = c(attr(predictions, 'models'), list(two_models))
      names(attr(predictions, 'models')) = paste0('n_nested', 0:n_nested)
    }

    return ( predictions )
  }
}

fit_two_outcome_models_long <- function(data, formula_t, formula_tmin1, family) {
  list(tmin1 = glm(formula_tmin1, family, data),
       t     = glm(formula_t,     family, data))

}

fit_outer_exp_model_long <- function(n_nested, data, preds, formula_n, binomial_n=NULL,
                                     obs_keep=NULL #this variable is used to select observations where t >= k, so we don't get into pre-time-zero lags. Is there a way to automate this?
) {

  if(is.null(binomial_n)) {
    data$lm_weights = 1 #these are 1's unless binomial aggregate data
  } else {
    data$lm_weights = binomial_n / sum(binomial_n)
  }

  data$preds = preds

  return (
    lm( formula = formula_n,
        data=data[obs_keep, ],
        weights = lm_weights)
  )
}
estimate_ice_long <- function(ice_preds, #(NxTt)x2 matrix
                              t_col, #(NxTt)x1 vector indicating which times the ice_preds correspond to
                              link_fun = NULL,
                              binomial_n = NULL) {
  if(is.null(link_fun)) {
    link_fun = function(x) x
  } else {
    stopifnot(is.function(link_fun))
  }

  if(is.null(binomial_n)) binomial_n = rep(1, nrow(ice_preds))
  freq_w = binomial_n / sum(binomial_n)

  df = tibble::tibble(t = t_col, ice_t = ice_preds[,2], ice_s = ice_preds[,1], freq_w=freq_w) %>%
    dplyr::group_by(t) %>%
    dplyr::summarise(estimate = link_fun(sum(ice_t * freq_w) / sum(freq_w)) - link_fun(sum(ice_s * freq_w) / sum(freq_w)))

  return ( df )
}
