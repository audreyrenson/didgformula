#' U.S. state-level stay-at-home and shelter-in-place orders during the COVID-19 pandemic in Spring 2020
#'
#' A dataset containing weekly stay-at-home/shelter-in-place order status for the 43 U.S. states that
#' implemented these orders in Spring 2020, along with all-cause mortality counts, COVID-19 case count growth
#' from the previous week, and state population counts.
#'
#'
#' @format A data frame with 516 rows and 6 variables:
#' \describe{
#'   \item{state}{State name}
#'   \item{week}{Week number, where week 0 is the week of April 5-11, 2020}
#'   \item{case_gr}{Number of new COVID-19 cases reported per 100k population in the state in the previous week}
#'   \item{stayathome}{Binary indicator of whether the state was under stay-at-home/shelter-in-place order at any time during the week (1=yes, 0=no)}
#'   \item{pop}{State population counts}
#' }
#' @source COVID-19 case counts come from \url{https://github.com/CSSEGISandData/COVID-19}, policy status
#' from \url{https://hrs.isr.umich.edu/data-products/restricted-data/available-products/10799}, mortality counts from
#' \url{https://data.cdc.gov/NCHS/Weekly-Provisional-Counts-of-Deaths-by-State-and-S/muzy-jte6}, and population
#' counts from the American Community Survey, 2015-2019.
"stayathome2020"
