#' Fit outcome models using `glm` at t and t-1 for ICE and outcome regression pipelines
#'
#' @param t int.
#' @param data wide data frame
#' @param rhs_formula_t chr. RHS glue-style formula for Yt model
#' @param rhs_formula_tmin chr. RHS glue-style formula for Yt-1 model
#' @param family stats::families object or string referring to one
#' @param binomial_n vector of group sizes of length(nrow(data))
#'
#' @return
#' @export
#'
#' @examples
fit_two_outcome_models <- function(t, data, rhs_formula_t, rhs_formula_tmin1, family, binomial_n=NULL) {
  #this is for estimating the innermost expectation of ICE, or for the outcome models in monte carlo
  #two at once because we typically want E[Y_t|L_t, A_t] and E[Y_{t-1}|L_t, A_t]

  if(is.null(binomial_n)) { #y is not aggregate binomial
    formulas =     c(glue_formula(paste0('Y{t-1}', rhs_formula_tmin1), t=t),
                     glue_formula(paste('Y{t}',rhs_formula_t ), t=t))
  } else {#y is aggregate binomial
    formulas =     c(glue_formula(paste0('cbind(Y{t-1}, binomial_n - Y{t-1})', rhs_formula_tmin1), t=t),
                     glue_formula(paste('cbind(Y{t}, binomial_n - Y{t})',rhs_formula_t ), t=t))
    data$binomial_n = binomial_n #and we don't want any ambiguity in case the 'data' already has a column 'binomial_n' (!!)
  }

  subst = data[[glue('A{t}')]] == 0

  models = lapply(1:2, function(s) glm(formula = as.formula(formulas[[s]]),
                                       family,
                                       data=data[subst, ]))

  names(models) = c('tmin1', 't')

  return (models)
}

#' Fit outer moderls for ICE using `lm`
#'
#' @param k int. Timepoint to which the outer expectation, i.e. E(mu_k^t|L_k, A_k), applies.
#' @param data wide data frame
#' @param preds vector of length(nrow(data)) of values predicted by the k+1 model.
#' @param rhs_formula_t chr. RHS glue-style formula for Yt model
#' @param rhs_formula_tmin chr. RHS glue-style formula for Yt-1 model
#' @param binomial_n vector of group sizes if Y refers to binomial aggregate data. Used to construct weights.
#'
#' @return
#' @export
#'
#' @examples
fit_outer_exp_model <- function(k, data, preds, rhs_formula, binomial_n=NULL) {
  #t is the timepoint for which we are estimating \hat\E[Y_t(\bar a) - Y_{t-1}(\bar a)]
  #k is the the timepoint the outer expectation applies to
  #preds are the k+1th outer expectation

  if(is.null(binomial_n)) {
    data$lm_weights = 1 #these are 1's unless binomial aggregate data
  } else {
    data$lm_weights = binomial_n / sum(binomial_n)
  }

  data$preds = preds

  subst = data[[glue('A{k}')]] == 0

  return (
    lm( formula = glue_formula(paste0('preds', rhs_formula), k=k),
        data=data[subst, ],
        weights = lm_weights)
  )
}

#' Recursively apply the ICE algorithm
#'
#' @param t int. Timepoint the estimand applies to
#' @param k int. Always 0 if called at the user level. Timepoint the outer expectation applies to
#' @param data wide data frame
#' @param inside_formula_t chr. glue-style formula for inside model for Yt (times indexed by t, e.g. '~L{t}')
#' @param inside_formula_tmin1 chr. glue-style formula for inside model for Yt (times indexed by t, e.g. '~L{t-1}')
#' @param outside_formula chr. glue-style formula for outside models (times indexed by k, e.g. '~L{k})
#' @param inside_family stats::families object (or chr) for glm inside models
#' @param binomial_n vector of length nrow(data) of group sizes if Y is binomial aggregate data
#' @param models lgl. Return all models as an attribute?
#'
#' @return
#' @export
#'
#' @examples
recursive_ice <- function(t, k, data, inside_formula_t, inside_formula_tmin1,
                          outside_formula, inside_family, binomial_n=NULL, models=TRUE) {

  if(k==t) {
    two_models = fit_two_outcome_models(t=k, data=data,
                                        rhs_formula_t = inside_formula_t,
                                        rhs_formula_tmin1 = inside_formula_tmin1,
                                        family=inside_family,
                                        binomial_n)
    predictions = sapply(two_models, predict, newdata=data, type='response', simplify=TRUE)

    if(models)  attr(predictions, 'models') = list(two_models)

    return (  predictions )
  } else {
    preds = recursive_ice(t,
                          k+1,
                          data,
                          inside_formula_t,
                          inside_formula_tmin1,
                          outside_formula,
                          inside_family,
                          binomial_n)

    two_models = lapply(1:2, function(h)  fit_outer_exp_model(k, data, preds[,h], outside_formula, binomial_n))
    predictions = sapply(two_models, predict, newdata=data, simplify=TRUE)

    if(models) {
      names(two_models) = c('tmin1','t')
      attr(predictions, 'models') = c(attr(preds,'models'), list(two_models))
      names(attr(predictions, 'models')) = paste0('n_nested', 0:(t-k))
    }


    return ( predictions )
  }
}

estimate_ice <- function(ice_preds,  # Nx2 matrix
                         link_fun,
                         binomial_n = NULL) {

  if(is.null(link_fun)) {
    link_fun = function(x) x
  } else {
    stopifnot(is.function(link_fun))
  }

  if(is.null(binomial_n)) binomial_n = rep(1, nrow(ice_preds))
  freq_w = binomial_n / sum(binomial_n)

  return ( diff(link_fun(colSums(ice_preds * freq_w) )) )
}

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
#' @param models logical. Return all models as an attribute?
#'
#' @return
#' @export
#'
#' @examples
ice_pipeline <- function(data,
                         inside_formula_t,
                         inside_formula_tmin1,
                         outside_formula,
                         Tt,
                         inside_family='gaussian',
                         pt_link_fun=NULL,
                         binomial_n=NULL,
                         tibble=TRUE,
                         models=TRUE) {

  ice_preds = lapply(1:Tt, function(t) recursive_ice(t, k=0,
                                                     data=data,
                                                     inside_formula_t=inside_formula_t,
                                                     inside_formula_tmin1=inside_formula_tmin1,
                                                     outside_formula=outside_formula,
                                                     inside_family=inside_family,
                                                     binomial_n=binomial_n))

  ice_estimates = sapply(ice_preds, estimate_ice, link_fun=pt_link_fun, binomial_n)

  if(!tibble) {
    result = ice_estimates
  } else {
    result = tibble::tibble(t=1:Tt,
                            estimate = ice_estimates)
  }

  if(models) {
    attr(result, 'models') = lapply(ice_preds, function(x) attr(x, 'models'))
    names(attr(result, 'models')) = paste0('t', 1:Tt)
  }

  return(result)

}


