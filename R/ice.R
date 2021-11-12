fit_two_outcome_models <- function(t, data, rhs_formula, family, binomial_n=NULL) {
  #this is for estimating the innermost expectation of ICE, or for the outcome models in monte carlo
  #two at once because we typically want E[Y_t|L_t, A_t] and E[Y_{t-1}|L_t, A_t]

  models = list()

  if(is.null(binomial_n)) { #y is not aggregate binomial
    formulas =     c(glue('Y{t-1}', rhs_formula),
                     glue('Y{t}', rhs_formula))
  } else {#y is aggregate binomial
    formulas =     c(glue('cbind(Y{t-1}, binomial_n - Y{t-1})', rhs_formula),
                     glue('cbind(Y{t}, binomial_n - Y{t})', rhs_formula))
    data$binomial_n = binomial_n #and we don't want any ambiguity in case the 'data' already has a column 'binomial_n' (!!)
  }

  models = lapply(1:2, function(s) glm(formula = as.formula(formulas[s]),
                                       family,
                                       data=data,
                                       subset=data[[glue::glue('A{t}')]] == 0))

  return (models)
}

fit_outer_exp_model <- function(k, data, preds, rhs_formula, binomial_n=NULL) {
  #t is the timepoint for which we are estimating \hat\E[Y_t(\bar a) - Y_{t-1}(\bar a)]
  #k is the the timepoint the outer expectation applies to
  #preds are the k+1th outer expectation

  lm_weights = if(is.null(binomial_n)) rep(1, nrow(data)) else  binomial_n / sum(binomial_n) #these are 1's unless binomial aggregate data

  return (
    lm( formula = as.formula(glue::glue('preds', rhs_formula)),
        data=data,
        subset=data[[glue::glue('A{k}')]] == 0,
        weights = lm_weights)
  )
}

recursive_ice <- function(t, k, data, inside_formula, outside_formula, inside_family, binomial_n=NULL) {
  if(k==t) {
    two_models = fit_two_outcome_models(t=k, data=data, rhs_formula = inside_formula, family=inside_family, binomial_n)
    predictions = sapply(two_models, predict, newdata=data, type='link', simplify=TRUE)
    return (  predictions[,2] - predictions[,1] )
  } else {
    preds = recursive_ice(t, k+1, data, inside_formula, outside_formula, inside_family, binomial_n)
    one_model = fit_outer_exp_model(k, data=data, preds=preds, rhs_formula=outside_formula, binomial_n)
    return (predict(one_model, newdata=data))
  }
}

estimate_ice <- function(ice_diffs, data, binomial_n = 1) {
  freq_w = binomial_n / ( binomial_n * sum(rep(1,nrow(data))) )

  return ( colSums(ice_diffs * freq_w) )
}

#' Iterated conditional DID g-formula estimator
#'
#' @param data Wide format data frame with one row per individual, and columns Yt for t = 0,1,...,Tt.
#' @param inside_formula chr, right-hand-side formula for inside models
#' @param outside_formula chr, right-hand-side formula for outside models
#' @param Tt int. max periods
#' @param tibble logical. return results as a tibble (TRUE) or vector (FALSE)?
#'
#' @return
#' @export
#'
#' @examples
ice_pipeline <- function(data, inside_formula, outside_formula, Tt, inside_family='gaussian', binomial_n=NULL, tibble=TRUE) {

  ice_diffs <- sapply(1:Tt,
                      recursive_ice,
                      k=0,
                      data=data,
                      inside_formula=inside_formula,
                      outside_formula=outside_formula,
                      inside_family=inside_family,
                      binomial_n=binomial_n)

  ice_estimates <- estimate_ice(ice_diffs, data=data)

  if(!tibble) {
    return (ice_estimates)
  } else {
    return ( tibble::tibble(t=1:Tt,
                            estimate = ice_estimates))
  }

}
