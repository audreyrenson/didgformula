test_that("parameter dimensions are correct", {
  Tt = 10
  Beta = generate_parameters(Tt, range_ymeans = qlogis(c(0.01, 0.05)), check_CDF=TRUE)

  check_dimensions <- function(Beta, Tt) {
    expect_equal(length(Beta$U), 1)
    expect_equal(dim(Beta$A), c(Tt + 1, 3))
    expect_equal(dim(Beta$L), c(Tt + 1, 2))
    expect_equal(dim(Beta$Y), c(Tt + 1, 5))
  }

  check_dimensions(Beta, Tt)
})

test_that('dydu is constant across t=1,2,...,Tt', {

  deriv_Beta <- function(Beta) {
    b1 = Beta[,1]
    b2 = Beta[,2]
    b3 = Beta[,3]
    b4 = Beta[,4]
    b5 = Beta[,5]
    u=1

    #expr = if(ylink == 'identity') expression(b1 + u*b2 + b3 + u*b5) else expression(1/ (1+exp(-b1-u*b2-b3-u*b5)))
    expr = expression(b1 + u*b2 + b3 + u*b5)

    dydu = D(expr, 'u')
    eval(dydu)
  }

  Beta_identity <- generate_parameters(Tt = 10)$Y
  derivs_identity <- deriv_Beta(Beta_identity)
  expect_equal( length(unique(derivs_identity)), 1)

  # Beta_logit <- generate_parameters(Tt = 10, ylink='logit')$Y
  # derivs_logit <- deriv_Beta(Beta_logit, ylink='logit')
  # expect_equal( length(unique(derivs_logit)), 1)

})

test_that('can generate binomial identity parameters iteratively', {
  set.seed(2)
  expect_true( length(generate_parameters_until_noerror(Tt=5, range_ymeans=c(0.4, 0.6))) > 0 )
})
