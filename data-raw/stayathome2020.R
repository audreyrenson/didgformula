## code to prepare `stayathome2020` dataset goes here




library(tidyverse)
library(lubridate)
library(splines)
library(didgformula)
theme_set(theme_minimal())


# load data ---------------------------------------------------------------
path = 'data-raw/stayathome2020_rawfiles/'
#from https://statepolicies.com/data/graphs/stay-at-home-order/
lockdown_dates_raw <- read_csv(glue::glue("{path}COVID-19 US state policy database (CUSP) - Stay at Home.csv"))

#from https://data.cdc.gov/NCHS/Weekly-Provisional-Counts-of-Deaths-by-State-and-S/muzy-jte6
mortality_weekly_raw <- read_csv(glue::glue("{path}Weekly_Provisional_Counts_of_Deaths_by_State_and_Select_Causes__2020-2021.csv"))

#case data
cases_daily <- read_csv(glue::glue("{path}time_series_covid19_confirmed_US.csv"))

#state denominators
state_pop <- read_csv(glue::glue("{path}nhgis0042_ds244_20195_2019_state.csv"))


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


cases_weekly <- cases_daily %>%
  select(state = Province_State, county = Admin2, `1/22/20`:`7/1/20`) %>%
  pivot_longer(-state:-county, names_to = "date", values_to="daily_cumulative_cases_by_county") %>%
  mutate(date = as.Date(date, format = "%m/%d/%y"),
         week_end = date + 1 - wday(date + 1)) %>%
  group_by(state, county, week_end) %>%
  mutate(weekly_cumulative_cases_by_county = max(daily_cumulative_cases_by_county)) %>%
  select(-daily_cumulative_cases_by_county, -date) %>%
  unique() %>%
  group_by(state, week_end) %>%
  summarise(weekly_cumulative_cases_by_state = sum(weekly_cumulative_cases_by_county)) %>%
  rename(cases_cumulative = weekly_cumulative_cases_by_state) %>%
  mutate(cases_abs_growth = cases_cumulative - lag(cases_cumulative),
         cases_perc_growth = 100 * (cases_abs_growth / cases_cumulative)) %>%
  ungroup()


state_pop_formatted <- state_pop %>%
  select(state = STATE, pop = ALUBE001)

full_weekly <- mortality_weekly %>%
  inner_join(lockdown_weekly) %>%
  inner_join(cases_weekly) %>%
  inner_join(state_pop_formatted) %>%
  group_by(state) %>%
  mutate(mort_rate_per100k = 1e5 * mort / pop,
         case_abs_growth_per100k = 1e5 * cases_abs_growth / pop,
         case_growth_lastweek = lag(case_abs_growth_per100k),
         case_growth_twoweeksago = lag(case_abs_growth_per100k, n=2)) %>%
  filter(!is.na(stayathome) & week_start >= as.Date("2020-04-05")) %>%
  mutate(weeknum = week(week_end) - min( week(week_end))) %>%
  ungroup()

stayathome2020 =full_weekly %>%
  mutate(week      = weeknum,
         case_gr = case_growth_lastweek %>%  ifelse(.<0.001, 0.001, .)) %>%#this is only for rhode island, t=10, a=1, this avoids NAs from I(A=0) in IPW formula
  select(state, week, case_gr, stayathome, mort, pop)



usethis::use_data(stayathome2020, overwrite = TRUE)
