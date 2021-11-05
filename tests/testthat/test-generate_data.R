test_that("data dimensions are correct", {

  N = 1e3
  Tt = 10
  df <- generate_data(N=N, Tt=Tt, Beta = generate_parameters(Tt=Tt))

  expect_equal(dim(df), c(N, 3*(Tt+1) + 2)) # 2 additional columns are U0 and uid

})
