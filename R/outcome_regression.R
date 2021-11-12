fit_covariate_models <- function(data, rhs_formula, family, Tt, binomial_n=NULL) {

  cov_weights = if(is.null(binomial_n)) rep(1, nrow(data)) else  binomial_n / sum(binomial_n) #these are 1's unless binomial aggregate data

  fit_one_model <- function(t) {
    return ( glm(formula = as.formula(glue('L{t}', rhs_formula)),
                 family = family,
                 data = data,
                 subset = data[[glue('A{t-1}')]] == 0,
                 weights = cov_weights) )
  }

  return (sapply(1:Tt, fit_one_model, simplify=FALSE))
}

fit_outcome_models <- function(data, rhs_formula, family, Tt, binomial_n=NULL) {
  return (sapply(1:Tt,
                 fit_two_outcome_models,
                 data=data,
                 rhs_formula=rhs_formula,
                 family=family,
                 binomial_n=binomial_n,
                 simplify=FALSE))
}

get_replicate_data <- function(data, nreps, binomial_n=NULL) {

  if(is.null(binomial_n)) binomial_n = rep(1, nrow(data))

  row_probs       = binomial_n/sum(binomial_n)   #to sample in proportion to binomial_n
  data$binomial_n = binomial_n                   #include this column in the resampled data for later summarizing
  indices         = sample(1:nrow(data), size=nreps, replace=TRUE, prob=row_probs)

  return ( data [ indices, c('uid','binomial_n', 'A0','L0','Y0')] )

}

simulate_Lt <- function(replicate_data, cov_model, t) {

  cov_preds = predict(cov_model, newdata=replicate_data, type='response')
  return ( rbinom(n=nrow(replicate_data), size=1, prob = cov_preds) )
}

simulate_Ydifft <- function(replicate_data,
                            outcome_models, #list of 2 glm objects, for t-1 and t
                            t) {
  outcome_preds = sapply(outcome_models, predict, newdata=replicate_data, type="link", simplify = TRUE)
  return ( outcome_preds[,2] - outcome_preds[,1] )
}

simulate_fulldata <- function(replicate_data,
                              cov_models, #list ot Tt glm objects
                              outcome_models, #list of Tt lists of 2 glm objects, each for t-1 and t, t=1,2,...,Tt
                              Tt) {
  for (t in 1:Tt) {
    replicate_data[[glue('L{t}')]] = simulate_Lt(replicate_data, cov_models[[t]], t)
    replicate_data[[glue('Ydiff{t}')]] = simulate_Ydifft(replicate_data, outcome_models[[t]], t)
  }

  return ( replicate_data )
}

estimate_or <- function(simulated_data, Tt, binomial_n=NULL) {

  if(is.null(binomial_n)) binomial_n = rep(1, nrow(simulated_data))
  freq_w = binomial_n / sum(binomial_n)  #so that weighted sum = weighted mean

  ydiff_names = sapply(1:Tt, function(t) glue('Ydiff{t}'))

  return (colSums (freq_w * simulated_data[, ydiff_names]) )
}

#' Outcome regression DID g-formula estimator
#'
#' @param data Wide format data frame with one row per individual, and columns Yt for t = 0,1,...,Tt.
#' @param y_formula chr. right-hand-side formula for outcome models
#' @param l_formula chr. right hand side formula for covariate models
#' @param Tt int. Max periods
#' @param nreps int. Number of reps for monte-carlo simulation
#' @param tibble logical. Return results as tibble (TRUE) or vector (FALSE)?
#'
#' @return
#' @export
#'
#' @examples
or_pipeline <- function(data, y_formula, l_formula, y_family='gaussian', l_family='binomial', Tt, nreps, binomial_n=NULL, tibble=TRUE) {

  cov_models = fit_covariate_models(data, rhs_formula = l_formula, family=l_family, Tt=Tt, binomial_n = binomial_n)
  out_models = fit_outcome_models(data, rhs_formula = y_formula, family=y_family, Tt=Tt, binomial_n = binomial_n)
  rep_data = get_replicate_data(data, nreps = nreps, binomial_n = binomial_n)
  sim_data = simulate_fulldata(rep_data, cov_models, out_models, Tt)
  or_estimates = estimate_or(sim_data, Tt, binomial_n=sim_data$binomial_n)

  if(!tibble) {
    return (or_estimates)
  } else {
    return ( tibble::tibble(t=1:Tt,
                            estimate = or_estimates))
  }
}

