test_that("ice pipeline returns results", {

  Tt=5
  N=1e3
  Beta=generate_parameters(Tt)
  df <- generate_data(N=N, Tt=Tt, Beta=Beta)

  ice_result = ice_pipeline(df, inside_formula_t = '~L{t}', inside_formula_tmin1 = '~L{t-1}',
                            outside_formula = '~L{k}', Tt=5)

  expect_false( any(is.na(ice_result)) )
})
