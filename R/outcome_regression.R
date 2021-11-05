fit_covariate_models <- function(data, rhs_formula, Tt) {
  fit_one_model <- function(t) {
    return ( glm(formula = as.formula(glue::glue('L{t}', rhs_formula)),
                 data = data,
                 subset = data[[glue::glue('A{t-1}')]] == 0) )
  }

  return (sapply(1:Tt, fit_one_model, simplify=FALSE))
}

fit_outcome_models <- function(data, rhs_formula, Tt) {
  return (sapply(1:Tt, fit_two_outcome_models, data=data, rhs_formula=rhs_formula, simplify=FALSE))
}

get_replicate_data <- function(data, nreps) {
  return ( data [ sample(1:nrow(data), size=nreps, replace=TRUE), c('uid', 'A0','L0','Y0')] )
}

simulate_Lt <- function(replicate_data, cov_model, t) {
  cov_preds = predict(cov_model, newdata=replicate_data, type='response')
  return ( rbinom(n=nrow(replicate_data), size=1, prob = cov_preds) )
}

simulate_Ydifft <- function(replicate_data,
                            outcome_models, #list of 2 glm objects, for t-1 and t
                            t) {
  outcome_preds = sapply(outcome_models, predict, newdata=replicate_data, type="response", simplify = TRUE)
  return ( outcome_preds[,2] - outcome_preds[,1] )
}

simulate_fulldata <- function(replicate_data,
                              cov_models, #list ot Tt glm objects
                              outcome_models, #list of Tt lists of 2 glm objects, each for t-1 and t, t=1,2,...,Tt
                              Tt) {
  for (t in 1:Tt) {
    replicate_data[[glue::glue('L{t}')]] = simulate_Lt(replicate_data, cov_models[[t]], t)
    replicate_data[[glue::glue('Ydiff{t}')]] = simulate_Ydifft(replicate_data, outcome_models[[t]], t)
  }

  return ( replicate_data )
}

estimate_or <- function(simulated_data, Tt) {
  ydiff_names = sapply(1:Tt, function(t) glue::glue('Ydiff{t}'))

  return (colMeans (simulated_data[, ydiff_names]) )
}

or_pipeline <- function(data, y_formula, l_formula, Tt, nreps, tibble=TRUE) {
  rep_data = get_replicate_data(data, nreps)
  cov_models = fit_covariate_models(data, rhs_formula = l_formula, Tt=Tt)
  out_models = fit_outcome_models(data, rhs_formula = y_formula, Tt=Tt)
  rep_data = get_replicate_data(data, nreps = nreps)
  sim_data = simulate_fulldata(rep_data, cov_models, out_models, Tt)
  or_estimates = estimate_or(sim_data, Tt)

  if(!tibble) {
    return (or_estimates)
  } else {
    return ( tibble::tibble(t=1:Tt,
                            estimate = or_estimates))
  }
}

or_bootstrap_se <- function(data, y_formula, l_formula, Tt, nreps, nboots, tibble=TRUE) {

  boot_data <- many_bootstraps(data, nboots, Tt)
  boot_estimates <- sapply(boot_data,
                           or_pipeline,
                           y_formula=y_formula,
                           l_formula=l_formula,
                           Tt=Tt,
                           nreps=nreps,
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

