#' Simulate from a fitted `glm` object based on new data, with automatic detection of family and link.
#'
#' @param model a fitted `glm` object
#' @param newdata data frame with values of the covariates based on which to simulate the response, similar to `predict.glm`
#' @param binomial_n int length nrow(newdata). If `model` was fit with `link=binomial`, you can optionally pass a vector of group sizes to simulate aggregate binomial data.
#'
#' @return A vector of length nrow(newdata) with simulated responses.
#' @export
#'
#' @examples
#' n  = 100
#' df = tibble::tibble(x = rnorm(n), y = rpois(n, lambda = exp(x)))
#' mod = glm(y~x, poisson('log'), df)
#' sims = sim(mod, newdata=df)
#' qqplot(df$y, sims)
#'
sim <- function(model, newdata, binomial_n=1) {

  if(!length(binomial_n) %in% c(1, nrow(newdata))) stop('binomial n must be length 1 or nrow(newdata)')

  N = nrow(newdata)
  n = binomial_n
  eta = predict(model, newdata=newdata, type='link')
  dispersion = summary(model)$dispersion
  simfun = get_simfun(family(model))

  simfun(eta, dispersion, n)

}

#' Generate a function that simulates data in accordance with the specification in a `stats::family` object
#'
#' @description This function is useful for simulating data that will follow a certain `glm` model
#'
#' @param family an object of type `stats::family`
#'
#' @return A function that can receive a linear predictor as input and return simulated samples based on the family and link specification in the `stats::family` object. The arguments of this function are: `eta`, vector of linear predictors; `dispersion`, the return value of `summary(model)$dispersion`, where `model` is a fitted `glm` object; and `n`, vector of group sizes (optional unless family=binomial)
#' @export
#'
#' @examples
#'
#' fam = Gamma('log')
#' simfun = get_simfun(fam)
#' sims = simfun(eta = 1, dispersion=1)
#'
get_simfun <- function(family) {
  function(eta, dispersion, n) {
    mu = family$linkinv(eta)
    switch ( family$family,
             'binomial' = rbinom(length(mu), n, mu),
             'gaussian' = rnorm(length(mu), mu, sqrt(dispersion)),
             'Gamma' = rgamma(length(mu), shape = 1/dispersion, scale = dispersion*mu),
             'poisson' = rpois(length(mu), lambda = mu) )
  }

}

