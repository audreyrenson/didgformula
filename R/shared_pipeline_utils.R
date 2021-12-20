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

  time_columns = model.matrix(as.formula(glue::glue('~-1 +factor({timevar})')), data=df)
  colnames(time_columns) = new_timevars

  result = dplyr::bind_cols(df, tibble::as_tibble(time_columns))
  attr(result, 'timevars') = new_timevars

  return(result)
}

append_lags = function(data, n_lags, lag_vars, default=0) { #perhaps add this to the package??
  for(n in 1:n_lags)
    for(lag_var in lag_vars)
      #this uses mutate because you can't take advantage of dplyr's groups otherwise
      data[[glue('{lag_var}_lag{n}')]] = dplyr::mutate(data, v=dplyr::lag(.data[[lag_var]], n, default))$v
    return (data)
}

