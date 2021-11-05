test_that("ice pipeline returns results", {

  Tt=5
  N=1e3
  df <- generate_data(N=N, Tt=Tt, Beta = generate_parameters(Tt=Tt))

  ice_result = ice_pipeline(df, inside_formula = '~L{t}', outside_formula = '~L{k}', Tt=5)

  expect_false( any(is.na(ice_result)) )
})
