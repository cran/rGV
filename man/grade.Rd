% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metrics.R
\name{grade}
\alias{grade}
\title{Calculate Glycemic Risk Assessment Diabetes Equation (GRADE)}
\usage{
grade(
  x,
  unit = "mg",
  method = "manuscript",
  c1 = ifelse(unit == "mg", 70.2, 3.9),
  c2 = ifelse(unit == "mg", 140.4, 7.8)
)
}
\arguments{
\item{x}{vector of glucose readings}

\item{unit}{"mg" if the units are mg/dL or "mmol" if the units are mmol/L. Null value is "mg".}

\item{method}{"manuscript" or "easy". Null value is "manuscript".}

\item{c1}{glucose value below which readings are considered hypoglycemic. Default is 70.2 mg/dL.}

\item{c2}{glucose value above which readings are considered hyperglycemic. Default is 140.4 mg/dL.}
}
\value{
A list containing the GRADE value and the percentage of the GRADE value due to euglycemia, hypoglycemia, and hyperglycemia for a given dataset of glucose measurements
}
\description{
Calculate Glycemic Risk Assessment Diabetes Equation (GRADE)
}
\examples{
grade(x=c(rep(100, 10), rep(120, 10), 105, 85), unit='mg', method='manuscript', c1=70.2, c2=140.4)
}
