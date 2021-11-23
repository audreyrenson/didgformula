test_that("outcome regression pipeline returns results", {
  Tt=5
  N=1e3
  Beta = generate_parameters(Tt=Tt)
  df <- generate_data(N=N, Tt=Tt, Beta = Beta)
  or_result = or_pipeline(df,
                          yt_formula = '~L{t}*L{t-1}*L{t-2}*L{t-3}',
                          ytmin1_formula='~L{t}*L{t-1}*L{t-2}*L{t-3}',
                          l_formula =  '~1',
                          Tt=5,
                          nreps=1e3)
  expect_false( any(is.na(or_result)) )
})
