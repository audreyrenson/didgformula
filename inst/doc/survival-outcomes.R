## ---- include = FALSE------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup-----------------------------------------------------------------------------------------------------------------
library(didgformula)
set.seed(777)

## --------------------------------------------------------------------------------------------------------------------------
N_obs = 1e4
time_periods = 10
parameters = generate_parameters(Tt=time_periods, range_ymeans = qlogis(c(0.01, 0.05)), check_CDF=TRUE)
parameters


## --------------------------------------------------------------------------------------------------------------------------
df_survival = generate_data(N=N_obs, Tt=time_periods, Beta=parameters, ylink="rbinom_logit_hazard")
head(df_survival)

## --------------------------------------------------------------------------------------------------------------------------
df_survival_po = generate_data(N=N_obs*10, Tt=time_periods, Beta=parameters, ylink="rbinom_logit_hazard", potential_outcomes = TRUE)
head(df_survival_po)

## --------------------------------------------------------------------------------------------------------------------------
pt_deviation_histograms(df_survival_po, Tt=time_periods, link_fun = qlogis)

## --------------------------------------------------------------------------------------------------------------------------
pt_deviation_histograms(df_survival_po, Tt=time_periods, k=5, link_fun = qlogis)

## --------------------------------------------------------------------------------------------------------------------------
estimates_ice   = ice_pipeline(df_survival, inside_formula_t = '~L{t}', inside_formula_tmin1='~L{t}', outside_formula = '~L{k}', inside_family = 'binomial', Tt=time_periods)
estimates_or    = or_pipeline(df_survival, yt_formula = '~L{t}', ytmin1_formula = '~L{t-1}', l_formula = '~1', y_family = 'binomial', Tt=time_periods, nreps=N_obs)
estimates_truth = estimate_truth(df_survival_po, Tt=time_periods, link_fun=qlogis)

combine_estimates(truth=estimates_truth, ice=estimates_ice, or=estimates_or, Tt=time_periods)

## --------------------------------------------------------------------------------------------------------------------------
gEY0_truth = qlogis(mean(df_survival_po$Y0))
gEY0_data  = qlogis(mean(df_survival$Y0))

gEYat_truth = tibble::tibble(t=1:time_periods, estimate = plogis(gEY0_truth + cumsum(estimates_truth$estimate)))
gEYat_ice = tibble::tibble(t=1:time_periods, estimate =  plogis(gEY0_data + cumsum(estimates_ice$estimate)))
gEYat_or = tibble::tibble(t=1:time_periods, estimate =  plogis(gEY0_data + cumsum(estimates_or$estimate)))

combined = combine_estimates(gEYat_truth, ice=gEYat_ice, or=gEYat_or, Tt=time_periods)
combined

## --------------------------------------------------------------------------------------------------------------------------
combined %>%
  dplyr::mutate(dplyr::across(truth:or, cumsum))

