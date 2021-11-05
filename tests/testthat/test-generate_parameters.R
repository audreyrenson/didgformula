test_that("parameter dimensions are correct", {
  Tt = 10
  Beta_id = generate_parameters(Tt)
  Beta_logit = generate_parameters(Tt, ylink = 'logit')

  check_dimensions <- function(Beta, Tt) {
    expect_equal(length(Beta$U), 1)
    expect_equal(dim(Beta$A), c(Tt + 1, 3))
    expect_equal(dim(Beta$L), c(Tt + 1, 2))
    expect_equal(dim(Beta$Y), c(Tt + 1, 5))
  }

  check_dimensions(Beta_id, Tt)
  check_dimensions(Beta_logit, Tt)
})

test_that('dydu is constant across t=1,2,...,Tt', {

  deriv_Beta <- function(Beta, ylink) {
    b1 = Beta[,1]
    b2 = Beta[,2]
    b3 = Beta[,3]
    b4 = Beta[,4]
    b5 = Beta[,5]
    u=1

    expr = if(ylink == 'identity') expression(b1 + u*b2 + b3 + u*b5) else expression(1/ (1+exp(-b1-u*b2-b3-u*b5)))

    dydu = D(expr, 'u')
    eval(dydu)
  }

  Beta_identity <- generate_parameters(Tt = 10, ylink='identity')$Y
  derivs_identity <- deriv_Beta(Beta_identity, ylink='identity')
  expect_equal( length(unique(derivs_identity)), 1)

  Beta_logit <- generate_parameters(Tt = 10, ylink='logit')$Y
  derivs_logit <- deriv_Beta(Beta_logit, ylink='logit')
  expect_equal( length(unique(derivs_logit)), 1)

})
