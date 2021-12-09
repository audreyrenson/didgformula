test_that("iptw pipeline returns previous results", {
  set.seed(5)
  Tt=5
  N=1e3
  Beta = generate_parameters(Tt=Tt)
  df <- generate_data(N=N, Tt=Tt, Beta = Beta)

  iptw_result = iptw_pipeline(df, Tt=5, den_formula='~L{t}')
  expect_snapshot_value(iptw_result, style='serialize')
})
