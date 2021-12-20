test_that("iptw_long returns same results as iptw", {

  set.seed(87)
  Tt=5
  Beta = generate_parameters(Tt=Tt)

  df_wide = generate_data(N=1e2, Tt=Tt, Beta=Beta, ylink = 'rbinom_logit', binomial_n = round(runif(1e2, min=1, max=50)))
  df_long = pivot_longer_gendata(df_wide) %>%
    dplyr::group_by(uid) %>%
    dplyr::mutate(A_lag = dplyr::lag(A),
                  Y_lag = dplyr::lag(Y)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(t>0) %>%
    make_time_dummies('t')

  df_interv <- df_long %>% dplyr::mutate(A=0, A_lag=0)

  tvars = paste(attr(df_long, 'timevars')[-1], collapse='+')

  den_formula = glue('A ~ ({tvars})*L*A_lag - (L*A_lag) + L')
  num_formula = glue('A ~ ({tvars})*A_lag - A_lag')

  iptw_long_result = iptw_pipeline_long(df_obs = df_long,
                                        df_interv = df_interv,
                                        den_formula = den_formula,
                                        family = binomial,
                                        yvar = 'Y',
                                        ylagvar = 'Y_lag',
                                        idvar = 'uid',
                                        timevar = 't',
                                        binomial_n = df_long$binomial_n,
                                        pt_link_fun = qlogis,
                                        models=FALSE)

  iptw_wide_result = iptw_pipeline(data = df_wide,
                                   Tt=Tt,
                                   den_formula = '~L{t}',
                                   pt_link_fun = qlogis,
                                   binomial_n=df_wide$binomial_n,
                                   models=FALSE)

  expect_equal( unname(iptw_long_result$estimate), unname(iptw_wide_result$estimate) )
})
