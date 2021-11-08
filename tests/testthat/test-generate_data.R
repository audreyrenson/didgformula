test_that("data dimensions are correct", {

  N = 1e3
  Tt = 10
  Beta = generate_parameters(Tt=Tt, range_ymeans=qlogis(c(0.01, 0.05)))

  correct_dimensions <- c(N, 3*(Tt+1) + 2) # Tt+1 periods of 3 variables (Lt,At,Yt), plus U0 and uid

  df_normal <- generate_data(N=N, Tt=Tt, Beta=Beta)
  df_binomial <- generate_data(N=N, Tt=Tt, Beta=Beta, ylink='rbinom_logit')
  df_hazard  <- generate_data(N=N, Tt=Tt, Beta=Beta, ylink='rbinom_logit_hazard')

  expect_equal(dim(df_normal), correct_dimensions)
  expect_equal(dim(df_binomial), correct_dimensions)
  expect_equal(dim(df_hazard), correct_dimensions)
})
