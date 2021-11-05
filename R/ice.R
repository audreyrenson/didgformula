fit_two_outcome_models <- function(t, data, rhs_formula) {
  #this is for estimating the innermost expectation of ICE, or for the outcome models in monte carlo
  #two at once because we typically want E[Y_t|L_t, A_t] and E[Y_{t-1}|L_t, A_t]

  models = list(2)

  models[[1]] = glm(formula = as.formula(glue::glue('Y{t-1}', rhs_formula)),
                    data=data,
                    subset=data[[glue::glue('A{t}')]] == 0)

  models[[2]] = glm(formula = as.formula(glue::glue('Y{t}', rhs_formula)),
                    data=data,
                    subset=data[[glue::glue('A{t}')]] == 0)
  return (models)
}

fit_outer_exp_model <- function(k, data, preds, rhs_formula) {
  #t is the timepoint for which we are estimating \hat\E[Y_t(\bar a) - Y_{t-1}(\bar a)]
  #k is the the timepoint the outer expectation applies to
  #preds are the k+1th outer expectation
  lm( formula = as.formula(glue::glue('preds', rhs_formula)),
      data=data,
      subset=data[[glue::glue('A{k}')]] == 0 )
}

recursive_ice <- function(t, k, data, inside_formula, outside_formula) {
  if(k==t) {
    two_models = fit_two_outcome_models(t=k, data=data, rhs_formula = inside_formula)
    predictions = sapply(two_models, predict, newdata=data, type='response', simplify=TRUE)
    return (  predictions[,2] - predictions[,1] )
  } else {
    preds = recursive_ice(t, k+1, data, inside_formula, outside_formula)
    one_model = fit_outer_exp_model(k, data=data, preds=preds, rhs_formula=outside_formula)
    return (predict(one_model, newdata=data))
  }
}

estimate_ice <- function(ice_diffs, data) {
  return ( colMeans(ice_diffs) )
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
ice_pipeline <- function(data, inside_formula, outside_formula, Tt, tibble=TRUE) {

  ice_diffs <- sapply(1:Tt,
                      recursive_ice,
                      k=0,
                      data=data,
                      inside_formula=inside_formula,
                      outside_formula=outside_formula)

  ice_estimates <- estimate_ice(ice_diffs, data=data)

  if(!tibble) {
    return (ice_estimates)
  } else {
    return ( tibble::tibble(t=1:Tt,
                            estimate = ice_estimates))
  }

}

ice_bootstrap_se <- function(data, inside_formula, outside_formula, Tt, nboots, tibble=TRUE) {

  boot_data <- many_bootstraps(data, nboots, Tt)
  boot_estimates <- sapply(boot_data,
                           ice_pipeline,
                           inside_formula=inside_formula,
                           outside_formula=outside_formula,
                           Tt=Tt,
                           tibble=FALSE,
                           simplify=FALSE)
  boot_se <- matrixStats::colSds(  Reduce(rbind, boot_estimates) )

  if(!tibble) {
    return (boot_se)
  } else {
    return ( tibble::tibble(t = 1:length(boot_se),
                            se = boot_se) )
  }

}
