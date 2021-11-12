#' Fit treatment models for IPTW
#'
#' @param data Wide format data frame with one row per individual, and columns Lt, At, Yt for t = 0,1,...,Tt.
#' @param rhs_formula chr. Formulas to use Format: '~L\{t\}+L\{t-1\}' etc.
#' @param Tt int. Final period in dataset (t=0,1,...,Tt)
#' @param ... Parameters passed to `glm`.
#'
#' @return list of models for period t=0,1,...,Tt, in that order.
#' @export
#'
#' @examples
fit_treatment_models <- function(data, rhs_formula, Tt, ...) {
  #returns a list of models


  fit_one_model <- function(t) {
    withCallingHandlers({
      glm(formula = as.formula(glue::glue('A{t}', rhs_formula)),
          family=binomial,
          data=data,
          subset=data[[glue::glue('A{t-1}')]] == 0,
          ...)
    }, warning = function(w) {
      #the non-integer successes warning happens anytime you weight a binomial or other discrete glm, and isn't actually a problem
      if (startsWith(conditionMessage(w), "non-integer #successes"))
        invokeRestart("muffleWarning")
    })
  }

  return (lapply(1:Tt, fit_one_model))

}


pred_treatment_models <- function(data, treatment_models) {
  #returns a matrix of predicted Pr[A{t}=0|\bar L{t}, \bar A{t-1} = 0]
  return (1 - sapply(treatment_models, predict, newdata=data, type='response', simplify=TRUE))

}

calc_weights <- function(denominator, numerator = 1) {
  #denominator is a NxTt matrix, numerator is a scalar or NxTt matrix
  matrixStats::rowCumprods( numerator / denominator )

}

#' Calculate Y_t - Y_t-1
#'
#' @param data Wide format data frame with one row per individual, and columns Yt for t = 0,1,...,Tt.
#' @param Tt int. Final period in dataset (t=0,1,...,Tt)
#'
#' @return matrix with nrow(data) rows and Tt columns, each with values of Y_t - Y_t-1 for t=1,2,...,Tt.
#' @export
#'
#' @examples
calc_ydiffs <- function(data, Tt) {
  #returns an NxTt matrix
  sapply(1:Tt, function(t) data[[glue::glue('Y{t}')]] - data[[glue::glue('Y{t-1}')]], simplify = TRUE)
}

calc_inclusion_indicators <- function(data, Tt) {
  #returns an NxTt matrix of 1s and 0s
  sapply(1:Tt, function(t) 1 * (data[[glue::glue('A{t}')]] == 0), simplify = TRUE)
}

estimate_iptw <- function(ydiffs, weights, inclusion_indicators, binomial_n = 1) {
  binomial_n = rep(1, nrow(ydiffs)) * binomial_n #get a column of 1's if this is not binomial aggregate data

  #returns a Tx1 vector of \hat\E[Y_t(\bar a) - Y_{t-1}(\bar a)], t=1,2,...,T
  return (colSums(inclusion_indicators * ydiffs * weights) / sum(binomial_n))
}

#' Estimate the DID g-formula using inverse-probability-of-treatment-weights.
#'
#' @param data Wide format data frame with one row per individual, and columns Lt, At, Yt for t = 0,1,...,Tt.
#' @param rhs_formula chr. Right-hand-side formula to use in predicting the treatment. Format: '~L\{t\}+L\{t-1\}' etc.
#' @param Tt int. Final period in dataset (t=0,1,...,Tt)
#' @param tibble logical. Should the results be returned as a tibble with columns (t, estimates) (TRUE) or a vector of just the estimates  (FALSE)?
#' @param binomial_n int length nrow(data). Group sizes for aggregate binomial data.
#'
#' @return Estimates of E(Yt(a) - Yt-1(a)), in the form of a tibble or vector (dependening on argument tibble), for times t=1,2,...,Tt (in that order).
#' @export
#'
#' @examples
iptw_pipeline <- function(data, rhs_formula, Tt, tibble=TRUE, binomial_n=1) {

  if(!length(binomial_n) %in% c(1, nrow(data))) stop('binomial_n must be lenth 1 or nrow(data)')
  binomial_n = binomial_n * rep(1, nrow(data)) # get a column of 1's if not aggregate binomial data

  den_mods = fit_treatment_models(data, rhs_formula, Tt, weights = binomial_n / sum(binomial_n))
  den_preds = pred_treatment_models(data, den_mods)
  weights = calc_weights(den_preds)
  estimates = estimate_iptw(ydiffs = calc_ydiffs(data, Tt),
                            weights = weights,
                            inclusion = calc_inclusion_indicators(data, Tt),
                            binomial_n = binomial_n)
  if(!tibble) {
    return (estimates)
  } else {
    return ( tibble::tibble(t = 1:length(estimates),
                            estimate = estimates) )
  }
}





