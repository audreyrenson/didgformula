## code to prepare `stayathome2020` dataset goes here




library(tidyverse)
library(lubridate)
library(glue)

# load data ---------------------------------------------------------------
path = 'data-raw/stayathome2020_rawfiles/'
#from https://statepolicies.com/data/graphs/stay-at-home-order/
lockdown_dates_raw <- read_csv(glue("{path}COVID-19 US state policy database (CUSP) - Stay at Home.csv"))

#from https://data.cdc.gov/NCHS/Weekly-Provisional-Counts-of-Deaths-by-State-and-S/muzy-jte6
mortality_weekly_raw <- read_csv(glue("{path}Weekly_Provisional_Counts_of_Deaths_by_State_and_Select_Causes__2020-2021.csv"))

#case data
cases_daily_by_county <- read_csv(glue("{path}time_series_covid19_confirmed_US.csv"))

#state denominators
state_pop <- read_csv(glue("{path}nhgis0042_ds244_20195_2019_state.csv"))


# Merge data --------------------------------------------------------------

lockdown_dates_formatted <- lockdown_dates_raw %>%
  mutate(start = as.Date(`Stay at home/shelter in place`, format = "%m/%d/%Y"),
         start_less_restrictive = as.Date(`Stay-at-home order issued but did not specifically restrict movement of the general public`,
                                          format = "%m/%d/%Y"),
         start = ifelse(is.na(start), start_less_restrictive, start), #choosing to include the less restrictive orders
         start = as.Date(start, origin = as.Date("1970-01-01")),
         stop = as.Date(`End stay at home/shelter in place`, format="%m/%d/%Y")) %>%
  select(state=State, start, stop)

mortality_weekly_plus <- mortality_weekly_raw %>%
  rename(state = `Jurisdiction of Occurrence`,
         week_end = `Week Ending Date`,
         mort = `All Cause`,
         mort_covid_mult = `COVID-19 (U071, Multiple Cause of Death)`,
         mort_covid_undr = `COVID-19 (U071, Underlying Cause of Death)`) %>%
  mutate(week_start = week_end - 6) %>%
  select(state, starts_with("week"), starts_with("mort"))
mortality_weekly_us <- mortality_weekly_plus %>%
  filter(state == "United States")
mortality_weekly <- mortality_weekly_plus %>%
  filter(! state %in% c("United States", "New York City", "Puerto Rico")) #don't have lockdown data on these jurisdictions (yet)




lockdown_weekly <- lockdown_dates_formatted %>%
  inner_join(mortality_weekly) %>%
  mutate(stayathome = 1*(start < week_end)*(stop > week_start)) %>%
  select(state, week_end, stayathome)


cases_weekly <- cases_daily_by_county %>%
  select(state = Province_State, county = Admin2, `1/22/20`:`7/1/20`) %>%
  pivot_longer(-state:-county, names_to = "date", values_to="weekly_cumulative_cases_by_county") %>%
  mutate(week_end = as.Date(date, format = "%m/%d/%y")) %>%
  filter(wday(week_end) == 7) %>%
  group_by(state, week_end) %>%
  summarise(weekly_cumulative_cases_by_state = sum(weekly_cumulative_cases_by_county)) %>%
  mutate(case_abs_growth_1wk = weekly_cumulative_cases_by_state - lag(weekly_cumulative_cases_by_state),
         case_abs_growth_4wk = weekly_cumulative_cases_by_state - lag(weekly_cumulative_cases_by_state, n=4)) %>%
  ungroup()


state_pop_formatted <- state_pop %>%
  select(state = STATE, pop = ALUBE001)

stayathome2020 <- mortality_weekly %>%
  inner_join(lockdown_weekly) %>%
  inner_join(cases_weekly) %>%
  inner_join(state_pop_formatted) %>%
  group_by(state) %>%
  mutate(case_abs_growth_1wk_per100k = lag(1e5 * case_abs_growth_1wk / pop), #these are lagged because we want the case growth in the week / 4 weeks previous
         case_abs_growth_4wk_per100k = lag(1e5 * case_abs_growth_4wk / pop)) %>%
  mutate(across(starts_with('case'), ~ifelse(.<0.001, 0.001, .))) %>% #this is only for rhode island, t=10, a=1, this avoids NAs from I(A=0) in IPW formula
  ungroup() %>%
  filter(!is.na(stayathome) & week_end >= as.Date("2020-04-05")) %>%
  mutate(week = week(week_end) - min( week(week_end))) %>%
  select(state, week, case_gr_1wk=case_abs_growth_1wk_per100k, case_gr_4wk=case_abs_growth_4wk_per100k, stayathome, mort, pop)

# # showing that the correlation is .91
# with(stayathome2020, cbind(case_gr_1wk, case_gr_4wk)) %>% log() %>% cor()

usethis::use_data(stayathome2020, overwrite = TRUE)
