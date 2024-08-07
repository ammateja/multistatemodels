% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{fit}
\alias{fit}
\title{Fit a multistate model}
\usage{
fit(
  model,
  multistatemodels,
  verbose = TRUE,
  compute_vcov = TRUE,
  constraints = NULL,
  optimize_surrogate = TRUE,
  surrogate_constraints = NULL,
  surrogate_parameters = NULL,
  maxiter = 200,
  tol = 0.001,
  alpha = 0.05,
  gamma = 0.05,
  kappa = 4/3,
  ess_target_initial = 100,
  MaxSamplingEffort = 20,
  npaths_additional = 10,
  return_ConvergenceRecords = TRUE,
  return_ProposedPaths = FALSE
)
}
\arguments{
\item{model}{A model created in the multistatemodel function}

\item{multistatemodels}{Loaded multistatemodels Julia environment}

\item{verbose}{logical; print messages, defaults to TRUE}

\item{compute_vcov}{logical; compute variance-covariance matrix, defaults to TRUE}

\item{constraints}{constraints on model parameters}

\item{optimize_surrogate}{logical; should the parameters Markov surrogate for
proposing paths be set to the MLE? defaults to TRUE}

\item{surrogate_constraints}{parameter constraints for the Markov surrogate}

\item{surrogate_parameters}{surrogate parameters}

\item{maxiter}{maximum number of MCEM iterations, defaults to 200}

\item{tol}{tolerance for the change in the MLL, i.e., upper bound of the
stopping rule to be ruled out, defaults to 1e-3}

\item{alpha}{standard normal quantile for asymptotic lower bound for ascent,
defaults to 0.05}

\item{gamma}{standard normal quantile for stopping the MCEM algorithm,
defaults to 0.05}

\item{kappa}{Inflation factor for target ESS per person, ESS_new = ESS_cur*kappa,
defaults to 4/3}

\item{ess_target_initial}{initial number of particles per participant for MCEM,
defaults to 100}

\item{MaxSamplingEffort}{factor of the ESS at which to break the loop for sampling additional paths,
defaults to 20}

\item{npaths_additional}{increment for number of additional paths when augmenting the pool of paths,
defaults to 10}

\item{return_ConvergenceRecords}{logical; save history throughout the run, defaults to TRUE}

\item{return_ProposedPaths}{logical; save latent paths and importance weights, defaults to FALSE}
}
\value{
A fitted multistate models environment
}
\description{
Fit a multistate model
}
\examples{
model_fitted <- fit(model, multistatemodels = multistatemodels, tol=1e-3, ess_target_initial=500)
}
