% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metrics.R
\name{adrr}
\alias{adrr}
\title{Calculate Average Daily Risk Range (ADRR)}
\usage{
adrr(x, times, unit = "mg", method = "manuscript")
}
\arguments{
\item{x}{vector of glucose readings}

\item{times}{vector of corresponding times, in minutes}

\item{unit}{"mg" if the units are mg/dL or "mmol" if the units are mmol/L. Null value is "mg".}

\item{method}{"manuscript", "corrected", or "easy". Null value is "manuscript".}
}
\value{
The numeric ADRR value for a given dataset of glucose measurements and times.
}
\description{
Calculate Average Daily Risk Range (ADRR)
}
\examples{
adrr(x=c(rep(100, 10), rep(120, 10), 105, 85), times=seq(0, 1260, 60),
unit="mg", method='manuscript')
}
