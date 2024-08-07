% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{compute_cumulative_hazard}
\alias{compute_cumulative_hazard}
\title{Compute cumulative hazard}
\usage{
compute_cumulative_hazard(
  tstart,
  tstop,
  model,
  hazard,
  multistatemodels,
  subj = 1
)
}
\arguments{
\item{tstart}{starting times}

\item{tstop}{stopping times}

\item{model}{multistate model}

\item{hazard}{string specifying the hazard, e.g., "h12" for the hazard for transitioning from state 1 to state 2.}

\item{multistatemodels}{Loaded multistatemodels Julia environment}

\item{subj}{subject id}
}
\value{
An integer with the cumulative hazard over (tsart, tstop)
}
\description{
Compute the cumulative hazard over (tstart,tstop).
}
\examples{
compute_cumulative_hazard(tstart=2, tstop=3, model=model, hazard="h12", multistatemodels=multistatemodels)
}
