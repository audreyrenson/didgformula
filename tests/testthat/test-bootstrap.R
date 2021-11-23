test_that("bootstrap_se returns results", {
  set.seed(143) #always show the computer some love

  Tt = 4
  N = 100
  df = generate_data(N = N,
                     Tt = Tt,
                     Beta = generate_parameters(Tt = 4))

  boot_iptw = bootstrap_se(
    data = df,
    nboots = 3,
    estimator = iptw_pipeline,
    rhs_formula = '~L{t}',
    Tt = Tt
  )
  boot_ice = bootstrap_se(
    data = df,
    nboots = 3,
    estimator = ice_pipeline,
    inside_formula_t = '~L{t}',
    inside_formula_tmin1 = '~L{t-1}',
    outside_formula = '~L{k}',
    Tt = Tt
  )
  boot_or = bootstrap_se(
    data = df,
    nboots = 3,
    estimator = or_pipeline,
    yt_formula = '~L{t}',
    ytmin1_formula = '~L{t-1}',
    l_formula = '~1',
    nreps = 100,
    Tt = Tt
  )

  expect_equal(dim(boot_iptw), c(Tt, 2))
  expect_equal(dim(boot_ice), c(Tt, 2))
  expect_equal(dim(boot_or), c(Tt, 2))
})
