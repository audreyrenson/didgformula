test_that("sim example produces the right number of sims", {
  n  = 100
  df = tibble::tibble(x = rnorm(n), y = rpois(n, lambda = exp(x)))
  mod = glm(y~x, poisson('log'), df)
  sims = sim(mod, newdata=df[1:10, ])

  expect_equal(length(sims), 10)
})
test_that('sim example produces no sims outside the support of response',{
  n  = 100
  df = tibble::tibble(x = rnorm(n), y = rpois(n, lambda = exp(x)))
  mod = glm(y~x, poisson('log'), df)
  sims = sim(mod, newdata=df[1:10, ])

  expect_true(all(sims) >= 0)
})
