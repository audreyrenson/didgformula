test_that("data dimensions are correct", {

  N = 1e3
  Tt = 10
  Beta = generate_parameters(Tt=Tt, range_ymeans=qlogis(c(0.01, 0.05)))

  correct_dimensions <- c(N, 4*(Tt+1) + 2) # Tt+1 periods of 4 variables (Lt,Wt,At,Yt), plus U0 and uid

  df_normal <- generate_data(N=N, Tt=Tt, Beta=Beta)
  df_binomial <- generate_data(N=N, Tt=Tt, Beta=Beta, ylink='rbinom_logit', binomial_n = round(runif(N, min=20, max=100)))
  df_hazard  <- generate_data(N=N, Tt=Tt, Beta=Beta, ylink='rbinom_logit_hazard')

  expect_equal(dim(df_normal), correct_dimensions)
  expect_equal(dim(df_binomial), correct_dimensions + c(0, 1)) #has additional binomial_n column, but same # rows
  expect_equal(dim(df_hazard), correct_dimensions)
})
