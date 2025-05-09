% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metric.R
\name{init_meta_regr}
\alias{init_meta_regr}
\title{Initialize Regression Meta Information}
\usage{
init_meta_regr()
}
\description{
This function initializes the metadata for regression metrics, checking whether
the `mlr3` package is installed. If `mlr3` is available, it loads the regression
metrics from `mlr3` and formats them. If not, it returns `NULL`.
}
\keyword{internal}
