## ---- include = FALSE------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup-----------------------------------------------------------------------------------------------------------------
library(didgformula)
library(tidyverse)
library(splines)
theme_set(theme_minimal())


## ----data------------------------------------------------------------------------------------------------------------------
data("stayathome2020")
stayathome2020

## ----recoding--------------------------------------------------------------------------------------------------------------

sah = stayathome2020
#the analysis is simpler if the exposure variable is 0 under the intervention in question
sah$open = 1-sah$stayathome

#Week 3 was the first time a stay-at-home ordder was lifted, so sah$open = 0 for everyone before that. 
#Thus, we code a variable that is the interaction of week*open that drops terms were sah$open is constant.
sah$week_by_open = ifelse(sah$week - 2 < 0, 0, sah$week - 2) * sah$open  

#coding the outcome as an Nx2 matrix for binomial likelihoods
sah$mort_bin = cbind(sah$mort, sah$pop - sah$mort)

#coding lags
sah$open_lag = sah %>% group_by(state) %>% mutate(x=dplyr::lag(open)) %>% pull()
sah$mort_lag = sah %>% group_by(state) %>% mutate(x=dplyr::lag(mort)) %>% pull()
sah$mort_bin_lag = cbind(sah$mort_lag, sah$pop - sah$mort_lag)
sah = sah %>% group_by(state) %>% append_lags(n_lags = 11, lag_vars = c('week_by_open', 'case_gr'), default = NA) %>% ungroup()

#most modeling uses a dataset starting with week=1, so create a separate dataset to simplify
sah_mod = sah[sah$week >0, ]

## ----estimation------------------------------------------------------------------------------------------------------------

# a function to set the exposure variable equal to the intervened status
intervene = function(df_obs) {
  df_obs %>% mutate(across(contains('open'), ~ifelse(is.na(.x), NA, 0)))
}

# convenience function for "positive part" of a variable, useful for keeping non-degenerative parts of time*exposure interaction terms
pos <- function(x) ifelse(x<0, 0, x)

ice_estimate = ice_pipeline_long(df_obs = sah_mod,
                                 df_interv = intervene(sah_mod),
                                 inside_formula_t = 'mort_bin ~ ns(week, df=3) + log(case_gr) + week_by_open ',
                                 inside_formula_tmin1 = 'mort_bin_lag ~ ns(week, df=3) + lag(case_gr) + week_by_open ',
                                 outside_formula = 'preds ~  ns(week, df={pos(Tt-n-5)+1}) + log(case_gr_lag{n}) + week_by_open_lag{n}',
                                 Tt = 11,
                                 t_col = sah_mod$week,
                                 n_nested = 7,
                                 inside_family = binomial,
                                 pt_link_fun = qlogis,
                                 binomial_n = sah_mod$pop)

iptw_estimate = iptw_pipeline_long(df_obs = sah_mod,
                                   df_interv = intervene(sah_mod),
                                   den_formula = 'open ~ ns(week, df=3) + log(case_gr)*open_lag',
                                   num_formula = 'open ~ ns(week, df=3) + open_lag',
                                   family = binomial,
                                   yvar = 'mort',
                                   ylagvar = 'mort_lag',
                                   idvar='state',
                                   timevar = 'week',
                                   pt_link_fun = qlogis,
                                   binomial_n = sah_mod$pop)
ice_estimate

iptw_estimate


## ----combining-------------------------------------------------------------------------------------------------------------


gEY0 = sah %>%
  filter(week==0) %>%
  summarise(gEY0 = qlogis(sum(mort)/sum(pop))) %>%
  pull()


make_rates = function(estimates, gEY0) {
  estimates %>% tibble::add_row(t=0, estimate=0, .before=1) %>%
    mutate(estimate = plogis( gEY0 + cumsum(estimate) ))
}

combine_estimates(iptw = iptw_estimate %>% make_rates(gEY0),
                  ice  = ice_estimate %>% make_rates(gEY0), Tt=11)



## ----compare_counts--------------------------------------------------------------------------------------------------------


make_death_counts = function(estimates, gEY0, pop) {
  estimates %>% make_rates(gEY0) %>% mutate(deaths = estimate * pop) %>% select(t, deaths)
}

make_lives_saved = function(estimates, data) {
  gEY0 = data %>%
    filter(week==0) %>%
    summarise(gEY0 = qlogis(sum(mort)/sum(pop))) %>%
    pull()

  counterfactual_deaths = make_death_counts(estimates, gEY0, pop=sum(data$pop[data$week==0])) %>% pull(deaths) %>% sum()
  factual_deaths = sum(data$mort)

  factual_deaths - counterfactual_deaths
}

make_lives_saved(ice_estimate, sah)
make_lives_saved(iptw_estimate, sah)

## --------------------------------------------------------------------------------------------------------------------------

R=5 #number of resamples

#an equivalent function to *_pipeline()s to estimate rates under the observed treatments.
obs_pipeline <- function(data) {
  data %>%
    filter(week>0) %>%
    group_by(t=week) %>%
    summarise(estimate = qlogis( sum(mort)/sum(pop)) - qlogis( sum(mort_lag) / sum(pop) ))
}

# bootstrap based on iid sampling -----------------------------------------

#add 'censored' observations - this is used in bootstrapping the aggregate data
boot_template = sah %>%
  bind_rows(sah %>%
              group_by(state, pop) %>%
              summarise(mort =  max(pop) - sum(mort)) %>%
              mutate(week = 12)) %>% #just have to remember that 12=never)
  arrange(state, week)


boot_agg = tibble(rep = 1:R) %>%
  mutate(data = map(rep, ~boot_template %>% #bootstrap the data with the 'censored' row
                      mutate(mort = rmultinom(1, sum(mort), prob=mort/sum(mort))[,1]) %>% 
                      group_by(state) %>%
                      mutate(pop = sum(mort),
                             mort_lag = dplyr::lag(mort),
                             mort_bin = cbind(mort, pop),
                             mort_bin_lag = cbind(mort_lag, pop)) %>%
                      ungroup() %>%
                      filter(week < 12)),
         mod_data = map(data, filter, week>0),
         iptw = map(mod_data, ~iptw_pipeline_long(df_obs = .,
                                                  df_interv = intervene(.),
                                                  den_formula = 'open ~ ns(week, df=3) + log(case_gr)*open_lag',
                                                  num_formula = 'open ~ ns(week, df=3) + open_lag',
                                                  family = binomial,
                                                  yvar = 'mort',
                                                  ylagvar = 'mort_lag',
                                                  idvar='state',
                                                  timevar = 'week',
                                                  pt_link_fun = qlogis,
                                                  binomial_n = .$pop)),
         ice = map(mod_data, ~ice_pipeline_long(df_obs = .,
                                                df_interv = intervene(.),
                                                inside_formula_t = 'mort_bin ~ ns(week, df=3) + log(case_gr) + week_by_open ',
                                                inside_formula_tmin1 = 'mort_bin_lag ~ ns(week, df=3) + lag(case_gr) + week_by_open ',
                                                outside_formula = 'preds ~  ns(week, df={pos(Tt-n-5)+1}) + log(case_gr_lag{n}) + week_by_open_lag{n}',
                                                Tt = 11,
                                                t_col = .$week,
                                                n_nested = 7,
                                                inside_family = binomial,
                                                pt_link_fun = qlogis,
                                                binomial_n = .$pop)),
         obs = map(mod_data, obs_pipeline))

#calculate standard errors for the relevant quantities
boot_ses = boot_agg %>%
  pivot_longer(iptw:obs, names_to='method') %>%
  mutate(gEY0 = map(data, ~with(.[.$week==0, ], qlogis(sum(mort)/sum(pop)))),
         rates = map2(value, gEY0, make_rates),
         lives = map2_dbl(value, data, make_lives_saved)) %>%
  unnest(rates) %>%
  group_by(t, method) %>%
  summarise(rate_mean = mean(estimate),
            rate_se = sd(estimate),
            lives_mean = mean(lives),
            lives_se = sd(lives))


## ----lives_estimates-------------------------------------------------------------------------------------------------------

lives_estimates_se = tibble(method = c('iptw','ice'),
                            lives = sapply(list(iptw_estimate, ice_estimate),
                                           make_lives_saved, sah)) %>%
  left_join(boot_ses %>% ungroup() %>% filter(method != 'obs') %>% select(method, se=lives_se) %>% unique()) %>%
  mutate(lwr = lives - 1.96*se, upr = lives + 1.96*se)

  
lives_estimates_se

## ----fig.width=6, fig.height=4---------------------------------------------------------------------------------------------
df_figure = list(ice=ice_estimate,
                 iptw=iptw_estimate,
                 obs = obs_pipeline(sah)) %>%
  map(make_rates, gEY0) %>%
  bind_rows(.id='method')  %>%
  left_join(boot_ses %>% select(method, t, se=rate_se)) %>%
  mutate(lwr = estimate - 1.96*se, upr = estimate+1.96*se)

df_figure %>%
  mutate(across(estimate:upr, ~.*1e5),
         Estimator = case_when(method == 'ice' ~'Iterated conditional',
                               method == 'iptw' ~'Inverse probability weighted',
                               method == 'obs' ~ 'Observed rates')) %>%
  ggplot(aes(x=t, y=estimate, ymin=lwr, ymax=upr, color=Estimator, fill=Estimator)) +
  geom_line() +
  geom_ribbon(alpha=0.2) +
  theme(legend.position = c(1,1), legend.justification = c(1,1)) +
  labs(x='Weeks since April 4, 2020', y='U.S. all-cause mortality rate per 100k (95% CI)')



