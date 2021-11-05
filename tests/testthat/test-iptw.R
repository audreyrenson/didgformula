test_that("iptw pipeline returns results", {
  Tt=5
  N=1e3
  df <- generate_data(N=N, Tt=Tt, Beta = generate_parameters(Tt=Tt))

  iptw_result = iptw_pipeline(df, '~L{t}', Tt=5)

  expect_false( any(is.na(iptw_result)) )
})
