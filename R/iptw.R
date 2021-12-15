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
fit_treatment_models <- function(data, rhs_formula, Tt, freq_w) {
  #returns a list of models

  fit_one_model <- function(t, data, rhs_formula, freq_w) {
    data$freq_w = freq_w #to avoid glm() being confused about its calling environment and allow us to manually subset the data (since 'subset' argument runs into problems)

    withCallingHandlers({

      glm(formula = glue_formula(paste0('A{t}', rhs_formula), t=t),
          family=binomial,
          data=data[data[[glue("A{t-1}")]]==0, ],
          weights=freq_w)

    }, warning = function(w) {
      #the non-integer successes warning happens anytime you weight a binomial or other discrete glm, and isn't actually a problem
      if (startsWith(conditionMessage(w), "non-integer #successes"))
        invokeRestart("muffleWarning")
    })
  }

  return (lapply(1:Tt, fit_one_model, data, rhs_formula, freq_w))

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

estimate_iptw <- function(data, Tt, weights, inclusion_indicators, link_fun=NULL, binomial_n = 1) {

  if(!length(binomial_n) %in% c(1, nrow(data))) stop('binomial_n must be lenth 1 or nrow(data)')
  binomial_n = rep(1, nrow(data)) * binomial_n #get a column of 1's if this is not binomial aggregate data

  if(is.null(link_fun)) {
    link_fun = function(x) x
  } else {
    stopifnot(is.function(link_fun))
  }

  yt = data %>% dplyr::select(Y1:glue::glue('Y{Tt}')) %>% as.matrix()
  ytmin1 = data %>% dplyr::select(Y0:glue::glue('Y{Tt-1}')) %>% as.matrix()
  gEyt = link_fun(colSums(inclusion_indicators * yt * weights) / colSums(inclusion_indicators * weights))
  gEytmin1 = link_fun(colSums(inclusion_indicators * ytmin1 * weights) /colSums(inclusion_indicators * weights))

  #returns a Tx1 vector of \hat\E[Y_t(\bar a) - Y_{t-1}(\bar a)], t=1,2,...,T
  return (gEyt - gEytmin1)
}

#' Estimate the DID g-formula using inverse-probability-of-treatment-weights.
#'
#' @param data Wide format data frame with one row per individual, and columns Lt, At, Yt for t = 0,1,...,Tt.
#' @param rhs_formula chr. Right-hand-side formula to use in predicting the treatment. Format: '~L\{t\}+L\{t-1\}' etc.
#' @param Tt int. Final period in dataset (t=0,1,...,Tt)
#' @param tibble logical. Should the results be returned as a tibble with columns (t, estimates) (TRUE) or a vector of just the estimates  (FALSE)?
#' @param pt_link_fun function. The scale on which parallel trends is assumed (e.g., `qlogis` for logit scale). Default `NULL` for untransformed scale.
#' @param binomial_n int length nrow(data). Group sizes for aggregate binomial data.
#' @param models logical. Return all models as an attribute?
#'
#' @return Estimates of E(Yt(a) - Yt-1(a)), in the form of a tibble or vector (dependening on argument tibble), for times t=1,2,...,Tt (in that order).
#' @export
#'
#' @examples
iptw_pipeline <- function(data, Tt, den_formula, num_formula=NULL, tibble=TRUE, pt_link_fun=NULL, binomial_n=1, models=TRUE) {

  if(!length(binomial_n) %in% c(1, nrow(data))) stop('binomial_n must be lenth 1 or nrow(data)')
  binomial_n = binomial_n * rep(1, nrow(data)) # get a column of 1's if not aggregate binomial data

  den_mods = fit_treatment_models(data, den_formula, Tt, freq_w = binomial_n / sum(binomial_n))
  den_preds = pred_treatment_models(data, den_mods)
  if(!is.null(num_formula)) {
    num_mods = fit_treatment_models(data, num_formula, Tt, freq_w = binomial_n / sum(binomial_n))
    num_preds = pred_treatment_models(data, num_mods)
    weights = calc_weights(denominator = den_preds, numerator = num_preds)
  } else {
    num_mods = NULL
    num_preds = 1
    weights = calc_weights(denominator = den_preds, numerator = 1)
  }

  inclusion_indicators = calc_inclusion_indicators(data, Tt)

  estimates = estimate_iptw(data=data,
                            Tt=Tt,
                            weights = weights,
                            inclusion =inclusion_indicators,
                            link_fun = pt_link_fun,
                            binomial_n = binomial_n)
  if(!tibble) {
    result = estimates
  } else {
    result = tibble::tibble(t = 1:length(estimates),
                            estimate = estimates)
  }

  if(models) attr(result, 'models') = list(den=den_mods, num=num_mods)

  return(result)
}




