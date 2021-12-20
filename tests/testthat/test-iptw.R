test_that("iptw pipeline returns previous results", {
  set.seed(5)
  Tt=5
  N=1e3
  Beta = generate_parameters(Tt=Tt)
  df <- generate_data(N=N, Tt=Tt, Beta = Beta, ylink = 'rbinom_logit', binomial_n = round(runif(N, min=1, max=100)))

  iptw_result = iptw_pipeline(df,
                              Tt=5,
                              den_formula='~L{t}',
                              num_formula='~1',
                              pt_link_fun = qlogis,
                              binomial_n = df$binomial_n,
                              models=FALSE)

  expect_snapshot_value(iptw_result, style='serialize')
})
