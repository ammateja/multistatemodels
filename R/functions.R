#' Check Julia setup
#'
#' @param ... Ignored
#'
#' @export
#'
#' @examples
#' multistatemodels_setup()
multistatemodels_setup <- function(...){
  if (JuliaConnectoR::juliaSetupOk()){
    JuliaConnectoR::juliaEval('
       import Pkg
       Pkg.add(url = "https://github.com/fintzij/MultistateModels.jl.git")
       Pkg.add("CSV")
       Pkg.add("DataFrames")
       Pkg.add("Distributions")
       Pkg.add("LinearAlgebra")
       Pkg.add("JLD2")
       Pkg.add("ArraysOfArrays")
       Pkg.add("Random")')
  }
  else {
    stop("Julia setup incorrect.
         Ensure Julia version >= 1.0 is properly installed.")
  }
}


#' Load Julia
#'
#' @description Load needed Julia packages. Run to use multistatemodels.
#'
#' @param ... Ignored
#'
#' @export
#'
#' @examples
#' multistatemodels <- multistatemodels_load()
multistatemodels_load <- function(...){
  tryCatch({
    JuliaConnectoR::juliaEval('
          using MultistateModels
          using CSV
          using DataFrames
          using Distributions
          using LinearAlgebra
          using JLD2
          using ArraysOfArrays
          using Random')
    multistatemodels <- JuliaConnectoR::juliaImport("MultistateModels")
    return(multistatemodels)
  },
  error=function(e) {
    stop("Run multistatemodels.multistatemodels_setup() once
            to complete installation and reload multistatemodels.")
  }
  )
}




check_data <- function(dat) {

  JuliaConnectoR::juliaLet('if any(names(dat)[1:6] .!== ["id", "tstart", "tstop", "statefrom", "stateto", "obstype"])
        error("The first 6 columns of the data should be id, tstart, tstop, statefrom, stateto, obstype.")
    end', dat=dat)

  JuliaConnectoR::juliaLet("dat.id = convert(Vector{Int64}, dat.id)
    dat.tstart = convert(Vector{Float64}, dat.tstart)
    dat.tstop = convert(Vector{Float64}, dat.tstop)
    dat.obstype = convert(Vector{Int64}, dat.obstype)
    dat.statefrom = convert(Vector{Union{Missing,Int64}}, dat.statefrom)
    dat.stateto = convert(Vector{Union{Missing, Int64}}, dat.stateto)", dat=dat)

  JuliaConnectoR::juliaLet('unique_id = unique(dat.id)
    nsubj = length(unique_id)
    if any(unique_id .!= 1:nsubj)
        error("The subject ids should be 1, 2, 3, ... .")
    end', dat=dat)

}

#' Collapse subjects and recompute sampling weights
#'
#' @description Collapse subjects to create an internal representation of a
#' dataset and optionally recompute a vector of sampling weights.
#'
#' @param dat A data frame
#' @param multistatemodels Loaded multistatemodels Julia environment
#' @param SamplingWeights Sampling weights, defaults to 1 for each id in the data
#'
#' @return A 2-element list. The first element is a data frame and the second
#' is a vector of sampling weights
#' @export
#'
#' @examples
#' multistatemodels <- multistatemodels_load()
#' df <- read.csv("~/derived_states_regen.csv")
#' d <- collapse_data(df, multistatemodels)
#' \dontshow{
#' JuliaConnectoR::stopJulia()
#' }
collapse_data <- function(dat, multistatemodels, SamplingWeights = rep(1, length(unique(dat$id)))) {

  df <- readr::format_csv(dat)

  df2 <- juliaLet("CSV.read(IOBuffer(df), DataFrame)", df=df)
  df2 <- JuliaConnectoR::juliaCall("DataFrame", df2)

  check_data(df2)

  d <- multistatemodels$collapse_data(df2, SamplingWeights = SamplingWeights)
  return(list(as.data.frame(d[[1]]), d[[2]]))

}


#' Specify hazard
#'
#' @description Specify a parametric or semi-parametric baseline cause-specific
#' hazard function.
#'
#'
#' @param formula Formula for log-hazard. Covariates have a multiplicative
#' effect on the baseline cause specific hazard.
#' Must be specified with "0 ~" on the left hand side.
#' @param statefrom Integer specifying the origin state
#' @param stateto Integer specifying the destination state
#' @param multistatemodels Loaded multistatemodels Julia environment
#' @param family string; One of "exp" or "wei" for exponential or Weibull cause-specific
#' baseline hazard functions, or "sp" for a semi-parametric spline basis for the
#' baseline hazard (defaults to M-splines).
#' @param degree For family="sp": degree of the spline polynomial basis
#' @param natural_spline For family="sp": logical; Restrict the second derivative to zero at the boundaries, defaults to TRUE.
#' @param extrapolation For family="sp": string; Either "linear" or "flat", defaults to "linear"
#' @param monotone For family="sp": 0, -1, or 1 for non-monotone, monotone decreasing, or monotone increasing.
#' @param boundaryknots For family="sp": Optional vector of boundary knots,
#' defaults to the range of possible sojourn times if not supplied.
#' @param knots For family="sp": Optional vector of knots. Defaults to the range of sojourns in
#' the data with no interior knots if not supplied.
#'
#' @return An environment with a Julia hazard
#' @export
#'
#' @examples
#' multistatemodels <- multistatemodels_load()
#' h12 <- Hazard(formula = 0 ~ 1+mab, family = "sp", statefrom = 1, stateto=2,
#' multistatemodels=multistatemodels, degree=1, knots=5/7)
#' h23 <- Hazard(formula = 0 ~ 1+mab, family = "exp", statefrom = 2, stateto=3,
#' multistatemodels=multistatemodels)
#' h24 <- Hazard(formula = 0 ~ 1+mab, family = "exp", statefrom = 2, stateto=4,
#' multistatemodels=multistatemodels)
#' h45 <- Hazard(formula = 0 ~ 1+mab, family = "exp", statefrom = 4, stateto=5,
#' multistatemodels=multistatemodels)
#' \dontshow{
#' JuliaConnectoR::stopJulia()
#' }
Hazard <- function(formula, statefrom, stateto, multistatemodels, family = c("exp", "wei", "sp"),
                   degree=3, natural_spline=T, extrapolation="linear", monotone=0, knots=NULL,
                   boundaryknots=NULL) {

  if (!family %in% c("exp", "wei", "sp")) {
    stop("family must be one of 'exp', 'wei', or 'sp'")
  }

  rhs <- paste(c(1, labels(terms(formula))), collapse = " + ")
  form <- JuliaConnectoR::juliaEval(paste0("@formula(", as.character(terms(formula)[[2]]), as.character(terms(formula)[[1]]), rhs, ")"))

  if (family == "sp") {
    knots <- JuliaConnectoR::juliaEval(paste0("vec([", knots, "])"))

    haz <- multistatemodels$Hazard(form, family, JuliaConnectoR::juliaCall("Int", statefrom), JuliaConnectoR::juliaCall("Int", stateto),
                                   degree = JuliaConnectoR::juliaCall("Int", degree), knots = knots, extrapolation = extrapolation, natural_spline = natural_spline,
                                   monotone = monotone, boundaryknots=boundaryknots)
  } else {
    haz <- multistatemodels$Hazard(form, family, JuliaConnectoR::juliaCall("Int", statefrom), JuliaConnectoR::juliaCall("Int", stateto))
  }

  return(haz)

}



#' Construct a multistate model
#'
#' @param hazard A vector of cause specific hazards
#' @param data Data frame
#' @param multistatemodels Loaded multistatemodels Julia environment
#' @param SamplingWeights Sampling weights
#' @param CensoringPatterns Censoring patterns
#' @param nPatterns Number of user-defined censoring patterns
#' @param verbose logical; print messages, defaults to FALSE
#'
#' @return A Julia MultistateModel environment that can be used for simulation and inference
#' @export
#'
#' @examples
#' CensoringPatterns <- matrix(c(3, 1, 0, 1, 0, 0), nrow=1, ncol=6)
#' model <- multistatemodel(hazard = c(h12, h23, h24, h45), data = d[[1]], SamplingWeights = d[[2]],
#' CensoringPatterns = CensoringPatterns, nPatterns = 1, multistatemodels =multistatemodels)

multistatemodel <- function(hazard, data, multistatemodels,
                            SamplingWeights=NULL, CensoringPatterns=NULL, nPatterns=NULL, verbose = FALSE) {

  if (!is.data.frame(data)) {
    stop("Data must be a data frame.")
  }

  if (!is.null(CensoringPatterns)) {

    CensoringPatterns <- matrix(as.integer(CensoringPatterns), nrow=nPatterns)

  }

  htuple <- JuliaConnectoR::juliaCall("Tuple", hazard)
  out <- "multistatemodel(htuple..., data=data,SamplingWeights=SamplingWeights, CensoringPatterns=CensoringPatterns, verbose=verbose)"
  mod <- JuliaConnectoR::juliaLet(out, htuple=htuple, data=JuliaConnectoR::juliaCall("DataFrame", data), SamplingWeights = SamplingWeights,
                                  CensoringPatterns = CensoringPatterns,
                                  verbose=verbose)
  #JuliaConnectoR::juliaGet(mod)
  return(mod)
}


#' Initialize parameters
#'
#' @description
#' Modify the parameter values in a MultistateProcess object,
#' calibrate to the MLE of a Markov surrogate.
#'
#' @param model A model created in the multistatemodel function
#' @param multistatemodels Loaded multistatemodels Julia environment
#' @param constraints Constraints on model parameters
#' @param surrogate_constraints Surrogate constraints
#' @param surrogate_parameters Surrogate parameters
#' @param crude logical; Defaults to FALSE
#'
#' @export
#'
#' @examples
#' initialize_parameters(model=model, multistatemodels = multistatemodels)
initialize_parameters <- function(model, multistatemodels, constraints = NULL, surrogate_constraints = NULL, surrogate_parameters=NULL, crude=FALSE) {

  multistatemodels$`initialize_parameters!`(model, constraints = constraints, surrogate_constraints = surrogate_constraints, surrogate_parameters=surrogate_parameters, crude=crude)

}


#' Fit a Markov surrogate model
#'
#' @param model A model created in the multistatemodel function
#' @param multistatemodels Loaded multistatemodels Julia environment
#' @param surrogate_parameters Surrogate parameters
#' @param surrogate_constraints Surrogate constraints on model parameters
#' @param crude_inits logical; Defaults to TRUE
#' @param verbose logical; print messages, defaults to TRUE
#'
#' @return A Julia environment with a fitted Markov surrogate model
#' @export
#'
#' @examples
#' model_fitted_surrogate <- fit_surrogate(model=model, multistatemodels=multistatemodels)
fit_surrogate <- function(model, multistatemodels, surrogate_parameters = NULL, surrogate_constraints = NULL, crude_inits = TRUE, verbose = TRUE) {

  multistatemodels$fit_surrogate(model, surrogate_parameters = surrogate_parameters, surrogate_constraints = surrogate_constraints, crude_inits = crude_inits, verbose = verbose)

}


#' Fit a multistate model
#'
#' @param model A model created in the multistatemodel function
#' @param multistatemodels Loaded multistatemodels Julia environment
#' @param verbose logical; print messages, defaults to TRUE
#' @param compute_vcov logical; compute variance-covariance matrix, defaults to TRUE
#' @param constraints constraints on model parameters
#' @param optimize_surrogate logical; should the parameters Markov surrogate for
#' proposing paths be set to the MLE? defaults to TRUE
#' @param surrogate_constraints parameter constraints for the Markov surrogate
#' @param surrogate_parameters surrogate parameters
#' @param maxiter maximum number of MCEM iterations, defaults to 200
#' @param tol tolerance for the change in the MLL, i.e., upper bound of the
#' stopping rule to be ruled out, defaults to 1e-3
#' @param alpha standard normal quantile for asymptotic lower bound for ascent,
#' defaults to 0.05
#' @param gamma standard normal quantile for stopping the MCEM algorithm,
#' defaults to 0.05
#' @param kappa Inflation factor for target ESS per person, ESS_new = ESS_cur*kappa,
#' defaults to 4/3
#' @param ess_target_initial initial number of particles per participant for MCEM,
#' defaults to 100
#' @param MaxSamplingEffort factor of the ESS at which to break the loop for sampling additional paths,
#' defaults to 20
#' @param npaths_additional increment for number of additional paths when augmenting the pool of paths,
#' defaults to 10
#' @param return_ConvergenceRecords logical; save history throughout the run, defaults to TRUE
#' @param return_ProposedPaths logical; save latent paths and importance weights, defaults to FALSE
#'
#' @return A fitted multistate models environment
#' @export
#'
#' @examples
#' model_fitted <- fit(model, multistatemodels = multistatemodels, tol=1e-3, ess_target_initial=500)
fit <- function(model, multistatemodels, verbose=TRUE, compute_vcov=TRUE, constraints = NULL,
                optimize_surrogate=TRUE, surrogate_constraints = NULL, surrogate_parameters = NULL, maxiter=200, tol=1e-3,
                alpha=0.05, gamma=0.05, kappa=4/3, ess_target_initial = 100, MaxSamplingEffort=20, npaths_additional = 10,
                return_ConvergenceRecords = TRUE, return_ProposedPaths = FALSE) {

  type <- juliaCall("typeof", model)[1]

  if (type %in% c("MultistateModels.MultistateMarkovModelCensored", "MultistateModels.MultistateMarkovModel", "MultistateModels.MultistateModel")) {

    multistatemodels$fit(model, verbose = verbose, compute_vcov = compute_vcov, constraints = constraints)

  } else {

    multistatemodels$fit(model, optimize_surrogate = optimize_surrogate, constraints = constraints, surrogate_constraints = surrogate_constraints, surrogate_parameters = surrogate_parameters,
                         maxiter=JuliaConnectoR::juliaCall("Int", maxiter), α=alpha, γ=gamma, κ=kappa, ess_target_initial = JuliaConnectoR::juliaCall("Int", ess_target_initial),
                         MaxSamplingEffort = JuliaConnectoR::juliaCall("Int", MaxSamplingEffort), npaths_additional = JuliaConnectoR::juliaCall("Int", npaths_additional),
                         return_ConvergenceRecords=return_ConvergenceRecords, return_ProposedPaths=return_ProposedPaths, compute_vcov=compute_vcov, verbose=verbose)

  }

}



#' AIC
#'
#' @description
#' Akaike's Information Criterion, defined as -2*log(L) + 2*k,
#' where L is the likelihood and k is the number of consumed degrees of freedom.
#'
#' @param model_fitted A fitted multistate models object
#' @param multistatemodels Loaded multistatemodels Julia environment
#' @param estimate_likelihood logical; whether to estimate the log-likelihood, defaults to TRUE.
#' @param min_ess minimum effective sample size per subject, defaults to 100.
#' @param loglik value of the log-likelihood to use, if provided.
#'
#' @return An integer with calculated AIC
#' @export
#'
#' @examples
#' aic(model_fitted=model_fitted, multistatemodels=multistatemodels)
aic <- function(model_fitted, multistatemodels, estimate_likelihood = TRUE, min_ess=100, loglik=NULL) {

  if (!JuliaConnectoR::juliaLet("isa(model_fitted, MultistateModels.MultistateModelFitted)", model_fitted = model_fitted)) {
    stop("Fitted model of type MultistateModels.MultistateModelFitted must be provided")
  }

  multistatemodels$aic(model_fitted, estimate_likelihood = estimate_likelihood, min_ess = JuliaConnectoR::juliaCall("Int",min_ess), loglik=loglik)

}

#' BIC
#'
#' @description
#' Bayesian Information Criterion, defined as -2*log(L) + k*log(n),
#' where L is the likelihood and k is the number of consumed degrees of freedom.
#'
#' @param model_fitted A fitted multistate models object
#' @param multistatemodels Loaded multistatemodels Julia environment
#' @param estimate_likelihood logical; whether to estimate the log-likelihood, defaults to TRUE.
#' @param min_ess minimum effective sample size per subject, defaults to 100.
#' @param loglik value of the log-likelihood to use, if provided.
#'
#' @return An integer with calculated BIC
#' @export
#'
#' @examples
#' bic(model_fitted=model_fitted, multistatemodels=multistatemodels)
bic <- function(model_fitted, multistatemodels, estimate_likelihood = TRUE, min_ess=100, loglik=NULL) {

  if (!JuliaConnectoR::juliaLet("isa(model_fitted, MultistateModels.MultistateModelFitted)", model_fitted = model_fitted)) {
    stop("Fitted model of type MultistateModels.MultistateModelFitted must be provided")
  }

  multistatemodels$bic(model_fitted, estimate_likelihood = estimate_likelihood, min_ess = JuliaConnectoR::juliaCall("Int",min_ess), loglik=loglik)

}

#' Get Convergence Records
#'
#' @description
#' Return the convergence records for the fit.
#'
#'
#' @param model_fitted A fitted multistate models object
#' @param multistatemodels Loaded multistatemodels Julia environment
#'
#' @return A Julia environment with status, candidate solution, algorithm,
#' convergence records, and work counters
#' @export
#'
#' @examples
#' get_ConvergenceRecords(model_fitted=model_fitted, multistatemodels=multistatemodels)
get_ConvergenceRecords <- function(model_fitted, multistatemodels) {

  if (!JuliaConnectoR::juliaLet("isa(model_fitted, MultistateModels.MultistateModelFitted)", model_fitted = model_fitted)) {
    stop("Fitted model of type MultistateModels.MultistateModelFitted must be provided")
  }

  multistatemodels$get_ConvergenceRecords(model_fitted)

}

#' Get Log-Likelihood
#'
#' @description
#'  Return the log-likelihood at the maximum likelihood estimates.
#'
#'
#' @param model_fitted A fitted multistate models object
#' @param multistatemodels Loaded multistatemodels Julia environment
#' @param ll string; one of "loglik" (default) for the observed data log-likelihood or
#' "subj_lml" for log marginal likelihood at the subject level
#'
#' @return An integer with the log-likelihood at the MLE
#' @export
#'
#' @examples
#' get_loglik(model_fitted=model_fitted, multistatemodels=multistatemodels)
get_loglik <- function(model_fitted, multistatemodels, ll="loglik") {

  if (!JuliaConnectoR::juliaLet("isa(model_fitted, MultistateModels.MultistateModelFitted)", model_fitted = model_fitted)) {
    stop("Fitted model of type MultistateModels.MultistateModelFitted must be provided")
  }

  if (!ll %in% c("loglik", "subj_lml")){
    stop("ll must be either 'loglik' or 'subj_lml'")
  }

  multistatemodels$get_loglik(model_fitted, ll=ll)

}

#' Get Parameters
#'
#' @description
#' Return the model parameters
#'
#'
#' @param model_fitted A fitted multistate models object
#' @param multistatemodels Loaded multistatemodels Julia environment
#'
#' @return A matrix of model parameters
#' @export
#'
#' @examples
#' get_parameters(model_fitted=model_fitted, multistatemodels=multistatemodels)
get_parameters <- function(model_fitted, multistatemodels) {

  if (!JuliaConnectoR::juliaLet("isa(model_fitted, MultistateModels.MultistateModelFitted)", model_fitted = model_fitted)) {
    stop("Fitted model of type MultistateModels.MultistateModelFitted must be provided")
  }

  pars <- multistatemodels$get_parameters(model_fitted)
  pars_list <- JuliaConnectoR::juliaGet(JuliaConnectoR::juliaCall("Tuple", pars))

  matrix(unlist(pars_list), nrow=length(pars_list), byrow=T)

}

#' Get Parameter Names
#'
#' @description
#'  Return the parameter names
#'
#'
#' @param model_fitted A fitted multistate models object
#' @param multistatemodels Loaded multistatemodels Julia environment
#'
#' @return A matrix with parameter names as strings
#' @export
#'
#' @examples
#' get_parnames(model_fitted=model_fitted, multistatemodels=multistatemodels)
get_parnames <- function(model_fitted, multistatemodels) {

  if (!JuliaConnectoR::juliaLet("isa(model_fitted, MultistateModels.MultistateModelFitted)", model_fitted = model_fitted)) {
    stop("Fitted model of type MultistateModels.MultistateModelFitted must be provided")
  }

  parnames <- multistatemodels$get_parnames(model_fitted)
  parnames_list <- JuliaConnectoR::juliaGet(JuliaConnectoR::juliaCall("Tuple", parnames))

  matrix(as.character(unlist(parnames_list)), nrow=length(parnames_list), byrow=T)

}

#' Get Variance-Covariance Matrix
#'
#' @description
#' Return the variance covariance matrix at the maximum likelihood estimate.
#'
#'
#' @param model_fitted A fitted multistate models object
#' @param multistatemodels Loaded multistatemodels Julia environment
#'
#' @return Variance-covariance matrix
#' @export
#'
#' @examples
#' get_vcov(model_fitted=model_fitted, multistatemodels=multistatemodels)
get_vcov <- function(model_fitted, multistatemodels) {

  if (!JuliaConnectoR::juliaLet("isa(model_fitted, MultistateModels.MultistateModelFitted)", model_fitted = model_fitted)) {
    stop("Fitted model of type MultistateModels.MultistateModelFitted must be provided")
  }

  multistatemodels$get_vcov(model_fitted)

}

#' Make parameter constraints
#'
#' @param cons string or vector of strings, name of constraints
#' @param lcons lcons
#' @param ucons ucons
#' @param multistatemodels Loaded multistatemodels Julia environment
#'
#' @return A Julia environment with a named Tuple of constraints
#' @export
#'
#' @examples
#' constraints_weibull = make_constraints(cons = "h12_shape", lcons = 0, ucons = 0, ultistatemodels=multistatemodels)
#' model_fitted=fit(model=model, multistatemodels=multistatemodels, constraints=constraints_weibull)

make_constraints <- function(cons, lcons, ucons, multistatemodels) {

  if((length(cons) != length(ucons)) | (length(cons) != length(lcons)) | (length(ucons) != length(lcons))) {
    stop("cons, lcons, and ucons must all be the same length")
  }

  if (length(cons) == 1) {
    cons <- JuliaConnectoR::juliaEval(paste0("vec([:(", cons, ";)])"))
    lcons <- JuliaConnectoR::juliaCall("float", JuliaConnectoR::juliaEval(paste0("vec([", lcons, "])")))
    ucons <- JuliaConnectoR::juliaCall("float", JuliaConnectoR::juliaEval(paste0("vec([", ucons, "])")))
    multistatemodels$make_constraints(cons = cons, lcons = lcons, ucons = ucons)
  } else {
    cons <- paste0(":(", cons, ")", collapse = ",")
    cons <- paste0("vec([", cons, "])")
    multistatemodels$make_constraints(cons=JuliaConnectoR::juliaEval(cons), lcons=lcons, ucons = ucons)
  }

}


#' Compute hazard
#'
#' @description
#' Compute the hazard at time t
#'
#'
#' @param t time or vector of times
#' @param model multistate model
#' @param hazard string specifying the hazard, e.g., "h12" for the hazard for transitioning from state 1 to state 2.
#' @param multistatemodels Loaded multistatemodels Julia environment
#' @param subj subject id
#'
#' @return An integer or vector of integers with the hazard at the specified time(s)
#' @export
#'
#' @examples
#' compute_hazard(t=2, model=model, hazard="h12", multistatemodels=multistatemodels)
compute_hazard <- function(t, model, hazard, multistatemodels, subj=1) {

  hazards <- NULL

  hazard <- paste0(":", hazard)

  for (i in 1:length(t)) {
    hazards <- append(hazards, multistatemodels$compute_hazard(t[i], model, JuliaConnectoR::juliaEval(hazard), JuliaConnectoR::juliaCall("Int", subj))[1])
  }

  return(hazards)

}

#' Compute cumulative hazard
#'
#' @description
#' Compute the cumulative hazard over (tstart,tstop).
#'
#'
#' @param tstart starting times
#' @param tstop stopping times
#' @param model multistate model
#' @param hazard string specifying the hazard, e.g., "h12" for the hazard for transitioning from state 1 to state 2.
#' @param multistatemodels Loaded multistatemodels Julia environment
#' @param subj subject id
#'
#' @return An integer with the cumulative hazard over (tsart, tstop)
#' @export
#'
#' @examples
#' compute_cumulative_hazard(tstart=2, tstop=3, model=model, hazard="h12", multistatemodels=multistatemodels)
compute_cumulative_hazard <- function(tstart, tstop, model, hazard, multistatemodels, subj = 1) {

  if (length(tstart) != length(tstop)) {
    stop("tstart and tstop must be the same length")
  }
  hazards <- NULL

  hazard <- paste0(":", hazard)

  for (i in 1:length(tstart)) {
    hazards <- append(hazards, multistatemodels$compute_cumulative_hazard(tstart[i], tstop[i], model, JuliaConnectoR::juliaEval(hazard), JuliaConnectoR::juliaCall("Int", subj))[1])
  }

  return(hazards)

}



#' Estimate log-likelihood
#'
#' @description
#' Estimate the log marginal likelihood for a fitted multistate model.
#' Require that the minimum effective sample size per subject is greater than min_ess.
#'
#' @param model multistate model
#' @param multistatemodels Loaded multistatemodels Julia environment
#' @param min_ess minimum effective sample size, defaults to 100.
#'
#' @return A Julia environment with a named Tuple
#' @export
#'
#' @examples
#' estimate_loglik(model=model, multistatemodels=multistatemodels)
estimate_loglik <- function(model, multistatemodels, min_ess = 100) {

  multistatemodels$estimate_loglik(model, min_ess = JuliaConnectoR::juliaCall("Int", min_ess))

}

#' Draw paths
#'
#' @description
#' Draw sample paths conditional on the data.
#' Require that the minimum effective sample size is greater than min_ess.
#'
#'
#' @param model multistate model
#' @param multistatemodels Loaded multistatemodels Julia environment
#' @param min_ess minimum effective sample size, defaults to 100.
#' @param paretosmooth logical; pareto smooth importance weights, defaults to TRUE unless min_ess < 25.
#' @param return_logliks logical; defaults to FALSE
#'
#' @return A Julia environment with a named Tuple of sample paths
#' @export
#'
#' @examples
#' draw_paths(model=model, multistatemodels = multistatemodels)
draw_paths <- function(model, multistatemodels, min_ess = 100, paretosmooth=TRUE, return_logliks = FALSE) {

  multistatemodels$draw_paths(model, min_ess = JuliaConnectoR::juliaCall("Int", min_ess), paretosmooth = paretosmooth, return_logliks = return_logliks)

}


#' Cumulative incidence
#'
#' @description
#' Compute the cumulative incidence for each possible transition as a function of time since state entry.
#' Assumes the subject starts their observation period at risk and saves cumulative incidence at the supplied vector of times, t.
#'
#' @param t time or a vector of times
#' @param model multistate model
#' @param multistatemodels Loaded multistatemodels Julia environment
#' @param subj subject id
#'
#' @return A matrix with the cumulative incidence where columns represent each
#' possible transition and rows represent t
#' @export
#'
#' @examples
#' cumulative_incidence(t=2, model=model, multistatemodels=multistatemodels)
cumulative_incidence <- function(t, model, multistatemodels, subj=1) {

  multistatemodels$cumulative_incidence(t, model, JuliaConnectoR::juliaCall("Int", subj))

}


#' Set parameters
#'
#' @description Set model parameters given a vector of values.
#' Copies newvalues to model.parameters.
#'
#' @param model multistate model
#' @param multistatemodels Loaded multistatemodels Julia environment
#' @param newvalues A list, each element of the list is new parameters for each hazard
#'
#' @export
#'
#' @examples
#' set_parameters(model, multistatemodels, newvalues = list(h12 = c(log(1.5), log(1)),
#' h23=c(log(2/3), log(1)), h24=c(log(1), log(1.25)), h45=c(log(1), log(1.25))))
set_parameters <- function(model, multistatemodels, newvalues) {

  values <- NULL

  for (i in 1:(length(newvalues)-1)) {

    x <- newvalues[[i]]
    x <- paste0(x, collapse = ",")
    x <- paste0("vec([", x, "]),")
    values <- paste0(values, x)
  }

  y <- newvalues[[length(newvalues)]]
  y <- paste0(y, collapse = ",")
  y <- paste0("vec([", y, "])")
  values <- paste0(values, y)



  multistatemodels$`set_parameters!`(model, JuliaConnectoR::juliaCall("Tuple", JuliaConnectoR::juliaEval(values)))

}


#' Simulate data
#'
#' @description
#' Simulate n data sets or collections of sample paths from a multistate model.
#' If data = TRUE (the default) discretely observed sample paths are returned, possibly subject to measurement error.
#' If paths = FALSE (the default), continuous-time sample paths are not returned.
#'
#' @param model multistate model
#' @param multistatemodels Loaded multistatemodels Julia environment
#' @param nsim number of sample paths to simulate
#' @param data logical; if TRUE then return discretely observed sample paths
#' @param paths logical; if FALSE then continuous-time sample paths not returned
#' @param delta_u minimum cumulative incidence increment per-jump, defaults to the larger of (1 / N subjects)^2 or sqrt(eps())
#' @param delta_t minimum time increment per-jump, defaults to sqrt(eps())
#'
#' @return A data frame with simulated data
#' @export
#'
#' @examples
#' simulate(model=model, multistatemodels=multistatemodels)
simulate <- function(model, multistatemodels, nsim=1, data=TRUE, paths=FALSE, delta_u = sqrt(.Machine$double.eps), delta_t = sqrt(.Machine$double.eps)) {

  sim_dat <- multistatemodels$simulate(model, nsim = JuliaConnectoR::juliaCall("Int", nsim), data=data, paths=paths, delta_u=delta_u, delta_t=delta_t)
  as.data.frame(sim_dat[[1]])

}


#' Summary of model output
#'
#' @param model_fitted fitted multistate models object
#' @param multistatemodels Loaded multistatemodels Julia environment
#' @param confidence_level confidence level of the confidence intervals, defaults to 0.95
#' @param estimate_likelihood logical; estimate likelihood, defaults to FALSE
#' @param min_ess minimum effective sample size, defaults to 100.
#'
#' @return A Julia environment with a named Tuple
#' Estimate, SE, upper and lower confidence limits for each hazard parameter
#' loglik, AIC, BIC, MCSE_loglik
#' @export
#'
#' @examples
#' model_summary(model_fitted=model_fitted, multistatemodels=multistatemodels, estimate_likelihood=TRUE)
model_summary <- function(model_fitted, multistatemodels, confidence_level = 0.95, estimate_likelihood=FALSE, min_ess = 100) {

  multistatemodels$summary(model_fitted, confidence_level = confidence_level, estimate_likelihood=estimate_likelihood, min_ess = JuliaConnectoR::juliaCall("Int", min_ess))

}





