% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{pt_deviation_histograms}
\alias{pt_deviation_histograms}
\title{View histograms of deviations from parallel trends in simulated data}
\usage{
pt_deviation_histograms(df_po, Tt, k = NULL, link_fun = NULL)
}
\arguments{
\item{df_po}{Dataset generated using generate_data(potential_outcomes=TRUE)}

\item{Tt}{int. max period.}

\item{k}{int. period index in conditioning event. If null, all periods are show.}

\item{link}{function. Corresponds to the scale of parallel trends. Typically either I (identity-scale) or qlogis (logit-scale)}
}
\description{
View histograms of deviations from parallel trends in simulated data
}
