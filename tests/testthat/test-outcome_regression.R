test_that("outcome regression pipeline returns results", {
  Tt=5
  N=1e3
  df <- generate_data(N=N, Tt=Tt, Beta = generate_parameters(Tt=Tt))
  or_result = or_pipeline(df, y_formula = '~L{t}', l_formula =  '~1', Tt=5, nreps=1e3)
  expect_false( any(is.na(or_result)) )
})
