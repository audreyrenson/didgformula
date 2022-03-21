test_that("iptw pipeline returns previous results", {
  set.seed(5)
  Tt=5
  N=1e3
  Beta = generate_parameters(Tt=Tt)
  df <- generate_data(N=N, Tt=Tt, Beta = Beta, ylink = 'rbinom_logit', binomial_n = round(runif(N, min=1, max=100)))
  iptw_result = iptw_pipeline(df,
                              Tt=5,
                              den_formula='~L{t}+W{t}+I(W{t}^2)',
                              num_formula='~1',
                              pt_link_fun = qlogis,
                              binomial_n = df$binomial_n,
                              models=FALSE)
  expect_equal(unname(iptw_result$estimate[3]), 3.5754744)

})
