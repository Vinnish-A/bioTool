% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metric.R
\name{add_custom_regr}
\alias{add_custom_regr}
\title{Add Custom Regression Metric}
\usage{
add_custom_regr(fun_, key_, description_ = "")
}
\arguments{
\item{fun_}{A function representing the custom regression metric.}

\item{key_}{A string representing the key for the custom metric.}

\item{description_}{A description of the custom metric (optional).}
}
\value{
NULL (invisible).
}
\description{
This function allows you to add a custom regression metric to the environment,
which can be used with other functions that expect regression metrics.
}
