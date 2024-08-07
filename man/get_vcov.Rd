% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{get_vcov}
\alias{get_vcov}
\title{Get Variance-Covariance Matrix}
\usage{
get_vcov(model_fitted, multistatemodels)
}
\arguments{
\item{model_fitted}{A fitted multistate models object}

\item{multistatemodels}{Loaded multistatemodels Julia environment}
}
\value{
Variance-covariance matrix
}
\description{
Return the variance covariance matrix at the maximum likelihood estimate.
}
\examples{
get_vcov(model_fitted=model_fitted, multistatemodels=multistatemodels)
}
