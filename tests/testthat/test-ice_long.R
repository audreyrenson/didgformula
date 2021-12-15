test_that("ice_long returns same result as ice", {
  Tt=3
  Beta = generate_parameters(Tt=Tt)


  set.seed(87)
  df_wide = generate_data(N=1e3, Tt=Tt, Beta=Beta)

  df_long = pivot_longer_gendata(df_wide) %>%
    dplyr::group_by(uid) %>%
    append_lags(n_lags = Tt, lag_vars = c('A','L'), default=NA) %>%
    dplyr::mutate(Y_lag = dplyr::lag(Y)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(t>0) %>%
    make_time_dummies('t')

  df_interv <- df_long %>% dplyr::mutate(dplyr::across(dplyr::starts_with('A'), ~0))
  tvars = attr(df_long, 'timevars')

  ice_preds = recursive_ice_long(Tt=Tt,
                     n_nested=Tt,
                     df_obs=df_long,
                     df_interv=df_interv,
                     inside_formula_t = 'Y~-1 + {tvars}*L*A - A*L',
                     inside_formula_tmin1 ='Y_lag~ -1 + {tvars}*L_lag1*A - L_lag1*A',
                     outside_formula = 'preds ~ -1 + {tvars}*L_lag{n}*A_lag{n} - L_lag{n}*A_lag{n} - t{n}:A_lag{n} - t{n}:L_lag{n}:A_lag{n} ',
                     inside_family=gaussian,
                     models=TRUE)


  ice_long = estimate_ice_long(ice_preds, t_col = df_long$t) #why are these slightly different?
  ice_wide = ice_pipeline(df_wide, '~L{t}', '~L{t-1}', '~L{k}', Tt)

  expect_equal(ice_long$estimate, ice_wide$estimate)

})
