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
#' Evaluate data against the implied density from a fitted `glm` object based on possibly new covariate values, with automatic detection of family and link.
#'
#' @param x vector of values at which to evaluate the density
#' @param model a fitted `glm` object
#' @param newdata data frame with values of the covariates based on which to simulate the response, similar to `predict.glm`. nrow(newdata) must equal length(x).
#' @param binomial_n int length nrow(newdata). If `model` was fit with `link=binomial`, you can optionally pass a vector of group sizes to simulate aggregate binomial data.
#'
#' @return A vector of length(x) of density estimates.
#' @export
#'
#' @examples
#' n  = 100
#' df = tibble::tibble(x = rnorm(n), y = rpois(n, lambda = exp(x)))
#' mod = glm(y~x, poisson('log'), df)
#' dens_y = dens(df$y, mod, newdata=df)
#'
#' # equivalent to
#' dpois(df$y, lambda = predict(mod, type='response'))
#'
dens <- function(x, model, newdata, binomial_n=1) {

  if(length(x) != nrow(newdata)) stop('length(x) must equal nrow(newdata)')
  if(!length(binomial_n) %in% c(1, length(x))) stop('binomial n must be length 1 or length(x)')

  N = length(x)
  n = binomial_n
  eta = predict(model, newdata=newdata, type='link')
  dispersion = summary(model)$dispersion
  densfun = get_densfun(family(model))

  densfun(x, eta, dispersion, n)
}
#' Generate a function that looks up a density in accordance with the specification in a `stats::family` object
#'
#' @description This function is useful for evaluating data against the implied density of a fitted `glm` object
#'
#' @param family an object of type `stats::family`
#'
#' @return A function with arguments `x` (data), `eta`, (linear predictor from `glm` object), `dispersion` (`summary(m)$dispersion` where `m` is a fitted `glm`), and  `n`, vector of group sizes (optional unless family=binomial)
#' @export
#'
#' @examples
#'
#' fam = Gamma('log')
#' densfun = get_densfun(fam)
#' densfun(0.2, eta = 2, dispersion=0.5)
#'
#' # the above is equivalent to:
#' dgamma(0.2, scale = exp(2), shape = 1/0.5)
get_densfun <- function(family) {
  function(x, eta, dispersion, n) {
    mu = family$linkinv(eta)
    stopifnot(length(x) == length(eta))

    switch ( family$family,
             'binomial' = dbinom(x, size=n, prob=mu),
             'gaussian' = dnorm(x, mu, sqrt(dispersion)),
             'Gamma' = dgamma(x, shape = 1/dispersion, scale = dispersion*mu),
             'poisson' = dpois(x, lambda = mu),
             'inverse.gaussian' = SuppDists::dinvGauss(x, nu=mu, lambda = 1/dispersion))
  }
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
             'poisson' = rpois(length(mu), lambda = mu),
             'inverse.gaussian' = SuppDists::rinvGauss(length(mu), nu=mu, lambda = 1/dispersion))
  }

}

