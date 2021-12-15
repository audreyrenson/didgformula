check_link_function <- function(link_fun) {
  if(is.null(link_fun)) {
    return (function(x) x)
  } else {
    stopifnot(is.function(link_fun))
    return( link_fun )
  }
}

make_time_dummies <- function(df, timevar) {
  new_timevars = paste0(timevar, min(df[[timevar]]):max(df[[timevar]]))

  result = df %>%
    dplyr::mutate(dummy = 1) %>%
    tidyr::pivot_wider(names_from = all_of(timevar), values_from = dummy, names_pre = 't') %>%
    dplyr::mutate(across(new_timevars, ~ifelse(is.na(.x), 0, .x)))

  result[[timevar]] = df[[timevar]] #keep the original

  attr(result, 'timevars') = new_timevars

  result
}

append_lags = function(data, n_lags, lag_vars, default=0) { #perhaps add this to the package??
  for(n in 1:n_lags)
    for(lag_var in lag_vars)
      #this uses mutate because you can't take advantage of dplyr's groups otherwise
      data[[glue('{lag_var}_lag{n}')]] = dplyr::mutate(data, v=dplyr::lag(.data[[lag_var]], n, default))$v
    return (data)
}
