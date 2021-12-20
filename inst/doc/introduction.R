## ---- include = FALSE----------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----parameters----------------------------------------------------------------------------------------------------------
library(didgformula)

set.seed(10)
N_obs = 1e4
time_periods = 5
parameters = generate_parameters(Tt=time_periods)
parameters


## ----generate_data-------------------------------------------------------------------------------------------------------
df_linear <- generate_data(N=N_obs, Tt=time_periods, Beta = parameters)
head(df_linear)

## ------------------------------------------------------------------------------------------------------------------------
df_po = generate_data(N=N_obs*10, Tt=time_periods, Beta=parameters,  potential_outcomes = TRUE)

truth = estimate_truth(df_po, Tt=time_periods)
truth

## ----histograms----------------------------------------------------------------------------------------------------------
pt_deviation_histograms(df_po, Tt=time_periods)


## ----histograms_k--------------------------------------------------------------------------------------------------------
pt_deviation_histograms(df_po, Tt=time_periods, k=2)


## ------------------------------------------------------------------------------------------------------------------------
estimates_iptw = iptw_pipeline(data = df_linear, den_formula = '~L{t}', num_formula = '~1', Tt=time_periods)
estimates_iptw

## ------------------------------------------------------------------------------------------------------------------------
estimates_ice = ice_pipeline(data = df_linear, inside_formula_t = '~L{t}', inside_formula_tmin1 = '~L{t-1}', outside_formula = '~L{k}', Tt=time_periods)
estimates_ice

## ------------------------------------------------------------------------------------------------------------------------
estimates_or = or_pipeline(data = df_linear, yt_formula = '~L{t}', ytmin1_formula = '~L{t}', l_formula = '~1', Tt=time_periods, nreps=N_obs) #usually nreps should be much larger but in this example it appears fine
estimates_or

## ------------------------------------------------------------------------------------------------------------------------
combine_estimates(truth, estimates_iptw, estimates_ice, estimates_or, Tt=time_periods)

## ------------------------------------------------------------------------------------------------------------------------
EY0_truth = mean(df_po$Y0)
EY0_data  = mean(df_linear$Y0)

EYat_truth = tibble::tibble(t=1:time_periods, estimate = EY0_truth + cumsum(truth$estimate))
EYat_iptw = tibble::tibble(t=1:time_periods, estimate = EY0_data + cumsum(estimates_iptw$estimate))
EYat_ice = tibble::tibble(t=1:time_periods, estimate = EY0_data + cumsum(estimates_ice$estimate))
EYat_or = tibble::tibble(t=1:time_periods, estimate = EY0_data + cumsum(estimates_or$estimate))

combine_estimates(EYat_truth, EYat_iptw, EYat_ice, EYat_or, Tt=time_periods)

## ------------------------------------------------------------------------------------------------------------------------

se_iptw = bootstrap_se(data = df_linear, 
                       nboots = 10, #obviously you should use many more than 10
                       estimator = iptw_pipeline,
                       den_formula = '~L{t}',
                       num_formula = '~1',
                       Tt = time_periods)

se_iptw                       


