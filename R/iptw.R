#' Fit treatment models for IPTW
#'
#' @param data Wide format data frame with one row per individual, and columns Lt, At, Yt for t = 0,1,...,Tt.
#' @param rhs_formula chr. Formulas to use Format '~L{t}+L{t-1}' etc.
#' @param Tt int. Final period in dataset (t=0,1,...,Tt)
#'
#' @return list of models for period t=0,1,...,Tt, in that order.
#' @export
#'
#' @examples
fit_treatment_models <- function(data, rhs_formula, Tt) {
  #returns a list of models

  fit_one_model <- function(t) glm(formula = as.formula(glue::glue('A{t}', rhs_formula)),
                                   family=binomial,
                                   data=data,
                                   subset=data[[glue::glue('A{t-1}')]] == 0)

  return (sapply(1:Tt, fit_one_model, simplify=FALSE))

}


pred_treatment_models <- function(data, treatment_models) {
  #returns a matrix of predicted Pr[A{t}=0|\bar L{t}, \bar A{t-1} = 0]
  return (1 - sapply(treatment_models, predict, newdata=data, type='response', simplify=TRUE))

}

calc_weights <- function(denominator, numerator = 1) {
  #denominator is a NxTt matrix, numerator is a scalar or NxTt matrix
  matrixStats::rowCumprods( numerator / denominator )

}

calc_ydiffs <- function(data, Tt) {
  #returns an NxTt matrix
  sapply(1:Tt, function(t) data[[glue::glue('Y{t}')]] - data[[glue::glue('Y{t-1}')]], simplify = TRUE)
}

calc_inclusion_indicators <- function(data, Tt) {
  #returns an NxTt matrix of 1s and 0s
  sapply(1:Tt, function(t) 1 * (data[[glue::glue('A{t}')]] == 0), simplify = TRUE)
}

estimate_iptw <- function(ydiffs, weights,inclusion_indicators) {
  #returns a Tx1 vector of \hat\E[Y_t(\bar a) - Y_{t-1}(\bar a)], t=1,2,...,T
  return (colMeans(inclusion_indicators * ydiffs * weights))
}

iptw_pipeline <- function(data, rhs_formula, Tt, tibble=TRUE) {
  den_mods = fit_treatment_models(data, rhs_formula, Tt)
  den_preds = pred_treatment_models(data, den_mods)
  weights = calc_weights(den_preds)
  estimates = estimate_iptw(ydiffs = calc_ydiffs(data, Tt),
                            weights = weights,
                            inclusion = calc_inclusion_indicators(data, Tt))
  if(!tibble) {
    return (estimates)
  } else {
    return ( tibble::tibble(t = 1:length(estimates),
                            estimate = estimates) )
  }
}





