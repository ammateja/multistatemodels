% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{multistatemodel}
\alias{multistatemodel}
\title{Construct a multistate model}
\usage{
multistatemodel(
  hazard,
  data,
  multistatemodels,
  SamplingWeights = NULL,
  CensoringPatterns = NULL,
  nPatterns = NULL,
  verbose = FALSE
)
}
\arguments{
\item{hazard}{A vector of cause specific hazards}

\item{data}{Data frame}

\item{multistatemodels}{Loaded multistatemodels Julia environment}

\item{SamplingWeights}{Sampling weights}

\item{CensoringPatterns}{Censoring patterns}

\item{nPatterns}{Number of user-defined censoring patterns}

\item{verbose}{logical; print messages, defaults to FALSE}
}
\value{
A Julia MultistateModel environment that can be used for simulation and inference
}
\description{
Construct a multistate model
}
\examples{
CensoringPatterns <- matrix(c(3, 1, 0, 1, 0, 0), nrow=1, ncol=6)
model <- multistatemodel(hazard = c(h12, h23, h24, h45), data = d[[1]], SamplingWeights = d[[2]],
CensoringPatterns = CensoringPatterns, nPatterns = 1, multistatemodels =multistatemodels)
}
