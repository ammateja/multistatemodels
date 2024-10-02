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

#' Setup Multistate Models with Julia packages
#'
#' @description
#' This function should be run any time there are updates to the Julia  packages
#' used (specifically MultistateModels).
#' This allows the necessary Julia packages to be updated to use in R.
#'
#' @export
#'
#' @examples
#' multistatemodels_setup()
multistatemodels_setup <- function() {

  if (JuliaConnectoR::juliaSetupOk()){
    JuliaConnectoR::juliaEval('
       import Pkg
       Pkg.add(url = "https://github.com/fintzij/MultistateModels.jl.git")
       Pkg.add("CSV")
       Pkg.add("DataFrames")')
  }
  else {
    stop("Julia setup incorrect.
         Ensure Julia version >= 1.10 is properly installed.")
  }

}

#' Collapse subjects and recompute sampling weights
#'
#' @description Collapse subjects to create an internal representation of a
#' dataset and optionally recompute a vector of sampling weights.
#'
#' @param dat A data frame
#' @param SamplingWeights Sampling weights, defaults to 1 for each id in the data
#'
#' @return A 2-element list. The first element is a data frame and the second
#' is a vector of sampling weights
#' @export
#'
#' @examples
#' d_collpased <- collapse_data(dat = illness_death_dat)
#' \dontshow{
#' JuliaConnectoR::stopJulia()
#' }
collapse_data <- function(dat, SamplingWeights = rep(1, length(unique(dat$id)))) {

  df <- readr::format_csv(dat)

  df2 <- JuliaConnectoR::juliaLet("CSV.read(IOBuffer(df), DataFrame)", df=df)
  df2 <- JuliaConnectoR::juliaCall("DataFrame", df2)

  check_data(df2)

  d <- multistatemodels_env$collapse_data(df2, SamplingWeights = SamplingWeights)
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
#' h12 <- Hazard(formula = 0~1, statefrom = 1, stateto=2, family="wei")
#' h13 <- Hazard(formula = 0~1, statefrom = 1, stateto=3, family="wei")
#' h23 <- Hazard(formula = 0~1, statefrom = 2, stateto=3, family="wei")
#' \dontshow{
#' JuliaConnectoR::stopJulia()
#' }
Hazard <- function(formula, statefrom, stateto, family = c("exp", "wei", "sp"),
                   degree=3, natural_spline=TRUE, extrapolation="linear", monotone=0, knots=NULL,
                   boundaryknots=NULL) {

  if (!family %in% c("exp", "wei", "sp")) {
    stop("family must be one of 'exp', 'wei', or 'sp'")
  }

  rhs <- paste(c(1, labels(stats::terms(formula))), collapse = " + ")
  form <- JuliaConnectoR::juliaEval(paste0("@formula(", as.character(stats::terms(formula)[[2]]), as.character(stats::terms(formula)[[1]]), rhs, ")"))

  if (family == "sp") {

    if (length(knots) == 1) {
      knots <- JuliaConnectoR::juliaEval(paste0("vec([", knots, "])"))
    }

    haz <- multistatemodels_env$Hazard(form, family, JuliaConnectoR::juliaCall("Int", statefrom), JuliaConnectoR::juliaCall("Int", stateto),
                                   degree = JuliaConnectoR::juliaCall("Int", degree), knots = knots, extrapolation = extrapolation, natural_spline = natural_spline,
                                   monotone = monotone, boundaryknots=boundaryknots)
  } else {
    haz <- multistatemodels_env$Hazard(form, family, JuliaConnectoR::juliaCall("Int", statefrom), JuliaConnectoR::juliaCall("Int", stateto))
  }

  return(haz)
}



#' Construct a multistate model
#'
#' @param hazard A vector of cause specific hazards
#' @param data Data frame
#' @param SamplingWeights Sampling weights
#' @param CensoringPatterns Censoring patterns
#' @param verbose logical; print messages, defaults to FALSE
#'
#' @return A Julia MultistateModel environment that can be used for simulation and inference
#' @export
#'
#' @examples
#' h12 <- Hazard(formula = 0~1, statefrom = 1, stateto=2, family="wei")
#' h13 <- Hazard(formula = 0~1, statefrom = 1, stateto=3, family="wei")
#' h23 <- Hazard(formula = 0~1, statefrom = 2, stateto=3, family="wei")
#' model_fit <- multistatemodel(hazard = c(h12, h13, h23), data=illness_death_dat)
#' \dontshow{
#' JuliaConnectoR::stopJulia()
#' }
multistatemodel <- function(hazard, data,
                            SamplingWeights=NULL, CensoringPatterns=NULL, verbose = FALSE) {

  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
    #stop("Data must be a data frame.")
  }

  if (!is.null(CensoringPatterns)) {

    CensoringPatterns <- matrix(as.integer(CensoringPatterns), nrow=nrow(CensoringPatterns))

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
#' @param constraints Constraints on model parameters
#' @param surrogate_constraints Surrogate constraints
#' @param surrogate_parameters Surrogate parameters
#' @param crude logical; Defaults to FALSE
#'
#' @export
#'
#' @examples
#' h12 <- Hazard(formula = 0~1, statefrom = 1, stateto=2, family="wei")
#' h13 <- Hazard(formula = 0~1, statefrom = 1, stateto=3, family="wei")
#' h23 <- Hazard(formula = 0~1, statefrom = 2, stateto=3, family="wei")
#' model_fit <- multistatemodel(hazard = c(h12, h13, h23), data=illness_death_dat)
#' model_fit <- initialize_parameters(model = model_fit)
#' \dontshow{
#' JuliaConnectoR::stopJulia()
#' }
initialize_parameters <- function(model, constraints = NULL, surrogate_constraints = NULL, surrogate_parameters=NULL, crude=FALSE) {

  #mod <- rlang::duplicate(model)
  multistatemodels_env$initialize_parameters(model, constraints = constraints, surrogate_constraints = surrogate_constraints, surrogate_parameters=surrogate_parameters, crude=crude)
  #return(mod)

}


#' Fit a Markov surrogate model
#'
#' @param model A model created in the multistatemodel function
#' @param surrogate_parameters Surrogate parameters
#' @param surrogate_constraints Surrogate constraints on model parameters
#' @param crude_inits logical; Defaults to TRUE
#' @param verbose logical; print messages, defaults to TRUE
#'
#' @return A Julia environment with a fitted Markov surrogate model
#' @export
#'
#' @examples
#' h12 <- Hazard(formula = 0~1, statefrom = 1, stateto=2, family="wei")
#' h13 <- Hazard(formula = 0~1, statefrom = 1, stateto=3, family="wei")
#' h23 <- Hazard(formula = 0~1, statefrom = 2, stateto=3, family="wei")
#' model_fit <- multistatemodel(hazard = c(h12, h13, h23), data=illness_death_dat)
#' model_fitted_surr <- fit_surrogate(model = model_fit)
#' \dontshow{
#' JuliaConnectoR::stopJulia()
#' }
fit_surrogate <- function(model, surrogate_parameters = NULL, surrogate_constraints = NULL, crude_inits = TRUE, verbose = TRUE) {

  multistatemodels_env$fit_surrogate(model, surrogate_parameters = surrogate_parameters, surrogate_constraints = surrogate_constraints, crude_inits = crude_inits, verbose = verbose)

}


#' Fit a multistate model
#'
#' @param model A model created in the multistatemodel function
#' @param verbose logical; print messages, defaults to TRUE
#' @param compute_vcov logical; compute variance-covariance matrix, defaults to TRUE
#' @param constraints constraints on model parameters
#' @param optimize_surrogate logical; should the parameters Markov surrogate for
#' proposing paths be set to the MLE? defaults to TRUE
#' @param surrogate_constraints parameter constraints for the Markov surrogate
#' @param surrogate_parameters surrogate parameters
#' @param maxiter maximum number of MCEM iterations, defaults to 100
#' @param tol tolerance for the change in the MLL, i.e., upper bound of the
#' stopping rule to be ruled out, defaults to 1e-2
#' @param ascent_threshold standard normal quantile for asymptotic lower bound for ascent,
#' defaults to 0.1
#' @param stopping_threshold standard normal quantile for stopping the MCEM algorithm,
#' defaults to 0.1
#' @param ess_increase Inflation factor for target ESS per person, ESS_new = ESS_cur*ess_increase,
#' defaults to 2
#' @param ess_target_initial initial number of particles per participant for MCEM,
#' defaults to 50
#' @param  max_ess maximum ess after which the mcem is stopped for nonconvergence, defaults to 10000
#' @param max_sampling_effort factor of the ESS at which to break the loop for sampling additional paths,
#' defaults to 20
#' @param npaths_additional increment for number of additional paths when augmenting the pool of paths,
#' defaults to 10
#' @param return_convergence_records logical; save history throughout the run, defaults to TRUE
#' @param return_proposed_paths logical; save latent paths and importance weights, defaults to FALSE
#' @param vcov_threshold if true, the variance covariance matrix calculation only inverts singular values
#' of the fisher information matrix that are greater than 1/sqrt(log(n) * k) where k is the number of parameters
#' and n is the number of subjects in the dataset. otherwise, the absolute tolerance is set to the square root of eps().
#' Defaults to TRUE
#'
#' @return A fitted multistate models environment
#' @export
#'
#' @examples
#' h12 <- Hazard(formula = 0~1, statefrom = 1, stateto=2, family="wei")
#' h13 <- Hazard(formula = 0~1, statefrom = 1, stateto=3, family="wei")
#' h23 <- Hazard(formula = 0~1, statefrom = 2, stateto=3, family="wei")
#' model_fit <- multistatemodel(hazard = c(h12, h13, h23), data=illness_death_dat)
#' model_fitted <- fit(model = model_fit, verbose=TRUE, compute_vcov = TRUE, ess_target_initial = 50,
#' ascent_threshold = 0.2, stopping_threshold = 0.2, tol = 0.001)
#' \dontshow{
#' JuliaConnectoR::stopJulia()
#' }
fit <- function(model, verbose=TRUE, compute_vcov=TRUE, constraints = NULL,
                optimize_surrogate=TRUE, surrogate_constraints = NULL, surrogate_parameters = NULL, maxiter=100, tol=1e-2,
                ascent_threshold=0.1, stopping_threshold=0.1, ess_increase=2, ess_target_initial = 50, max_ess=10000, max_sampling_effort=20, npaths_additional = 10,
                return_convergence_records = TRUE, return_proposed_paths = FALSE, vcov_threshold = TRUE) {

  type <- JuliaConnectoR::juliaCall("typeof", model)[1]

  if (type %in% c("MultistateModels.MultistateMarkovModelCensored", "MultistateModels.MultistateMarkovModel", "MultistateModels.MultistateModel")) {

    multistatemodels_env$fit(model, verbose = verbose, compute_vcov = compute_vcov, constraints = constraints)

  } else {

    multistatemodels_env$fit(model, optimize_surrogate = optimize_surrogate, constraints = constraints, surrogate_constraints = surrogate_constraints, surrogate_parameters = surrogate_parameters,
                         maxiter=JuliaConnectoR::juliaCall("Int", maxiter), ascent_threshold=ascent_threshold, stopping_threshold=stopping_threshold, ess_increase=ess_increase, ess_target_initial = JuliaConnectoR::juliaCall("Int", ess_target_initial),
                         max_sampling_effort = JuliaConnectoR::juliaCall("Int", max_sampling_effort), npaths_additional = JuliaConnectoR::juliaCall("Int", npaths_additional), max_ess = JuliaConnectoR::juliaCall("Int", max_ess),
                         return_convergence_records=return_convergence_records, return_proposed_paths=return_proposed_paths, compute_vcov=compute_vcov, verbose=verbose, vcov_threshold = vcov_threshold)

  }

}



#' AIC
#'
#' @description
#' Akaike's Information Criterion, defined as -2*log(L) + 2*k,
#' where L is the likelihood and k is the number of consumed degrees of freedom.
#'
#' @param model_fitted A fitted multistate models object
#' @param estimate_likelihood logical; whether to estimate the log-likelihood, defaults to TRUE.
#' @param min_ess minimum effective sample size per subject, defaults to 100.
#' @param loglik value of the log-likelihood to use, if provided.
#'
#' @return An integer with calculated AIC
#' @export
#'
#' @examples
#' h12 <- Hazard(formula = 0~1, statefrom = 1, stateto=2, family="wei")
#' h13 <- Hazard(formula = 0~1, statefrom = 1, stateto=3, family="wei")
#' h23 <- Hazard(formula = 0~1, statefrom = 2, stateto=3, family="wei")
#' model_fit <- multistatemodel(hazard = c(h12, h13, h23), data=illness_death_dat)
#' model_fitted <- fit(model = model_fit, verbose=TRUE, compute_vcov = TRUE, ess_target_initial = 50,
#' ascent_threshold = 0.2, stopping_threshold = 0.2, tol = 0.001)
#' aic(model_fitted = model_fitted)
#' \dontshow{
#' JuliaConnectoR::stopJulia()
#' }
aic <- function(model_fitted, estimate_likelihood = TRUE, min_ess=100, loglik=NULL) {

  if (!JuliaConnectoR::juliaLet("isa(model_fitted, MultistateModels.MultistateModelFitted)", model_fitted = model_fitted)) {
    stop("Fitted model of type MultistateModels.MultistateModelFitted must be provided")
  }

  multistatemodels_env$aic(model_fitted, estimate_likelihood = estimate_likelihood, min_ess = JuliaConnectoR::juliaCall("Int",min_ess), loglik=loglik)

}

#' BIC
#'
#' @description
#' Bayesian Information Criterion, defined as -2*log(L) + k*log(n),
#' where L is the likelihood and k is the number of consumed degrees of freedom.
#'
#' @param model_fitted A fitted multistate models object
#' @param estimate_likelihood logical; whether to estimate the log-likelihood, defaults to TRUE.
#' @param min_ess minimum effective sample size per subject, defaults to 100.
#' @param loglik value of the log-likelihood to use, if provided.
#'
#' @return An integer with calculated BIC
#' @export
#'
#' @examples
#' h12 <- Hazard(formula = 0~1, statefrom = 1, stateto=2, family="wei")
#' h13 <- Hazard(formula = 0~1, statefrom = 1, stateto=3, family="wei")
#' h23 <- Hazard(formula = 0~1, statefrom = 2, stateto=3, family="wei")
#' model_fit <- multistatemodel(hazard = c(h12, h13, h23), data=illness_death_dat)
#' model_fitted <- fit(model = model_fit, verbose=TRUE, compute_vcov = TRUE, ess_target_initial = 50,
#' ascent_threshold = 0.2, stopping_threshold = 0.2, tol = 0.001)
#' bic(model_fitted = model_fitted)
#' \dontshow{
#' JuliaConnectoR::stopJulia()
#' }
bic <- function(model_fitted, estimate_likelihood = TRUE, min_ess=100, loglik=NULL) {

  if (!JuliaConnectoR::juliaLet("isa(model_fitted, MultistateModels.MultistateModelFitted)", model_fitted = model_fitted)) {
    stop("Fitted model of type MultistateModels.MultistateModelFitted must be provided")
  }

  multistatemodels_env$bic(model_fitted, estimate_likelihood = estimate_likelihood, min_ess = JuliaConnectoR::juliaCall("Int",min_ess), loglik=loglik)

}

#' Get Convergence Records
#'
#' @description
#' Return the convergence records for the fit.
#'
#'
#' @param model_fitted A fitted multistate models object
#'
#' @return A Julia environment with status, candidate solution, algorithm,
#' convergence records, and work counters
#' @export
#'
#' @examples
#' h12 <- Hazard(formula = 0~1, statefrom = 1, stateto=2, family="wei")
#' h13 <- Hazard(formula = 0~1, statefrom = 1, stateto=3, family="wei")
#' h23 <- Hazard(formula = 0~1, statefrom = 2, stateto=3, family="wei")
#' model_fit <- multistatemodel(hazard = c(h12, h13, h23), data=illness_death_dat)
#' model_fitted <- fit(model = model_fit, verbose=TRUE, compute_vcov = TRUE, ess_target_initial = 50,
#' ascent_threshold = 0.2, stopping_threshold = 0.2, tol = 0.001)
#' get_convergence_records(model_fitted=model_fitted)
#' \dontshow{
#' JuliaConnectoR::stopJulia()
#' }
get_convergence_records <- function(model_fitted) {

  if (!JuliaConnectoR::juliaLet("isa(model_fitted, MultistateModels.MultistateModelFitted)", model_fitted = model_fitted)) {
    stop("Fitted model of type MultistateModels.MultistateModelFitted must be provided")
  }

  multistatemodels_env$get_convergence_records(model_fitted)

}

#' Get Log-Likelihood
#'
#' @description
#'  Return the log-likelihood at the maximum likelihood estimates.
#'
#'
#' @param model_fitted A fitted multistate models object
#' @param ll string; one of "loglik" (default) for the observed data log-likelihood or
#' "subj_lml" for log marginal likelihood at the subject level
#'
#' @return An integer with the log-likelihood at the MLE
#' @export
#'
#' @examples
#' h12 <- Hazard(formula = 0~1, statefrom = 1, stateto=2, family="wei")
#' h13 <- Hazard(formula = 0~1, statefrom = 1, stateto=3, family="wei")
#' h23 <- Hazard(formula = 0~1, statefrom = 2, stateto=3, family="wei")
#' model_fit <- multistatemodel(hazard = c(h12, h13, h23), data=illness_death_dat)
#' model_fitted <- fit(model = model_fit, verbose=TRUE, compute_vcov = TRUE, ess_target_initial = 50,
#' ascent_threshold = 0.2, stopping_threshold = 0.2, tol = 0.001)
#' get_loglik(model_fitted=model_fitted)
#' \dontshow{
#' JuliaConnectoR::stopJulia()
#' }
get_loglik <- function(model_fitted, ll="loglik") {

  if (!JuliaConnectoR::juliaLet("isa(model_fitted, MultistateModels.MultistateModelFitted)", model_fitted = model_fitted)) {
    stop("Fitted model of type MultistateModels.MultistateModelFitted must be provided")
  }

  if (!ll %in% c("loglik", "subj_lml")){
    stop("ll must be either 'loglik' or 'subj_lml'")
  }

  multistatemodels_env$get_loglik(model_fitted, ll=ll)

}

#' Get Parameters
#'
#' @description
#' Return the model parameters
#'
#'
#' @param model_fitted A fitted multistate models object
#'
#' @return A matrix of model parameters
#' @export
#'
#' @examples
#' h12 <- Hazard(formula = 0~1, statefrom = 1, stateto=2, family="wei")
#' h13 <- Hazard(formula = 0~1, statefrom = 1, stateto=3, family="wei")
#' h23 <- Hazard(formula = 0~1, statefrom = 2, stateto=3, family="wei")
#' model_fit <- multistatemodel(hazard = c(h12, h13, h23), data=illness_death_dat)
#' model_fitted <- fit(model = model_fit, verbose=TRUE, compute_vcov = TRUE, ess_target_initial = 50,
#' ascent_threshold = 0.2, stopping_threshold = 0.2, tol = 0.001)
#' get_parameters(model_fitted=model_fitted)
#' \dontshow{
#' JuliaConnectoR::stopJulia()
#' }
get_parameters <- function(model_fitted) {

  if (!JuliaConnectoR::juliaLet("isa(model_fitted, MultistateModels.MultistateModelFitted)", model_fitted = model_fitted)) {
    stop("Fitted model of type MultistateModels.MultistateModelFitted must be provided")
  }

  pars <- multistatemodels_env$get_parameters(model_fitted)
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
#'
#' @return A matrix with parameter names as strings
#' @export
#'
#' @examples
#' h12 <- Hazard(formula = 0~1, statefrom = 1, stateto=2, family="wei")
#' h13 <- Hazard(formula = 0~1, statefrom = 1, stateto=3, family="wei")
#' h23 <- Hazard(formula = 0~1, statefrom = 2, stateto=3, family="wei")
#' model_fit <- multistatemodel(hazard = c(h12, h13, h23), data=illness_death_dat)
#' model_fitted <- fit(model = model_fit, verbose=TRUE, compute_vcov = TRUE, ess_target_initial = 50,
#' ascent_threshold = 0.2, stopping_threshold = 0.2, tol = 0.001)
#' get_parnames(model_fitted=model_fitted)
#' \dontshow{
#' JuliaConnectoR::stopJulia()
#' }
get_parnames <- function(model_fitted) {

  if (!JuliaConnectoR::juliaLet("isa(model_fitted, MultistateModels.MultistateModelFitted)", model_fitted = model_fitted)) {
    stop("Fitted model of type MultistateModels.MultistateModelFitted must be provided")
  }

  parnames <- multistatemodels_env$get_parnames(model_fitted)
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
#'
#' @return Variance-covariance matrix
#' @export
#'
#' @examples
#' h12 <- Hazard(formula = 0~1, statefrom = 1, stateto=2, family="wei")
#' h13 <- Hazard(formula = 0~1, statefrom = 1, stateto=3, family="wei")
#' h23 <- Hazard(formula = 0~1, statefrom = 2, stateto=3, family="wei")
#' model_fit <- multistatemodel(hazard = c(h12, h13, h23), data=illness_death_dat)
#' model_fitted <- fit(model = model_fit, verbose=TRUE, compute_vcov = TRUE, ess_target_initial = 50,
#' ascent_threshold = 0.2, stopping_threshold = 0.2, tol = 0.001)
#' get_vcov(model_fitted=model_fitted)
#' \dontshow{
#' JuliaConnectoR::stopJulia()
#' }
get_vcov <- function(model_fitted) {

  if (!JuliaConnectoR::juliaLet("isa(model_fitted, MultistateModels.MultistateModelFitted)", model_fitted = model_fitted)) {
    stop("Fitted model of type MultistateModels.MultistateModelFitted must be provided")
  }

  multistatemodels_env$get_vcov(model_fitted)

}

#' Make parameter constraints
#'
#' @param cons string or vector of strings, name of constraints
#' @param lcons lcons
#' @param ucons ucons
#'
#' @return A Julia environment with a named Tuple of constraints
#' @export
#'
#' @examples
#' constraints_weibull <- make_constraints(cons = "h12_shape", lcons = 0, ucons = 0)
#' \dontshow{
#' JuliaConnectoR::stopJulia()
#' }
make_constraints <- function(cons, lcons, ucons) {

  if((length(cons) != length(ucons)) | (length(cons) != length(lcons)) | (length(ucons) != length(lcons))) {
    stop("cons, lcons, and ucons must all be the same length")
  }

  if (length(cons) == 1) {
    cons <- JuliaConnectoR::juliaEval(paste0("vec([:(", cons, ";)])"))
    lcons <- JuliaConnectoR::juliaCall("float", JuliaConnectoR::juliaEval(paste0("vec([", lcons, "])")))
    ucons <- JuliaConnectoR::juliaCall("float", JuliaConnectoR::juliaEval(paste0("vec([", ucons, "])")))
    multistatemodels_env$make_constraints(cons = cons, lcons = lcons, ucons = ucons)
  } else {
    cons <- paste0(":(", cons, ")", collapse = ",")
    cons <- paste0("vec([", cons, "])")
    multistatemodels_env$make_constraints(cons=JuliaConnectoR::juliaEval(cons), lcons=lcons, ucons = ucons)
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
#' @param subj subject id
#'
#' @return An integer or vector of integers with the hazard at the specified time(s)
#' @export
#'
#' @examples
#' h12 <- Hazard(formula = 0~1, statefrom = 1, stateto=2, family="wei")
#' h13 <- Hazard(formula = 0~1, statefrom = 1, stateto=3, family="wei")
#' h23 <- Hazard(formula = 0~1, statefrom = 2, stateto=3, family="wei")
#' model_fit <- multistatemodel(hazard = c(h12, h13, h23), data=illness_death_dat)
#' compute_hazard(t=0.5, model=model_fit, hazard = "h12")
#' \dontshow{
#' JuliaConnectoR::stopJulia()
#' }
compute_hazard <- function(t, model, hazard, subj=1) {

  hazards <- NULL

  hazard <- paste0(":", hazard)

  for (i in 1:length(t)) {
    hazards <- append(hazards, multistatemodels_env$compute_hazard(t[i], model, JuliaConnectoR::juliaEval(hazard), JuliaConnectoR::juliaCall("Int", subj))[1])
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
#' @param subj subject id
#'
#' @return An integer with the cumulative hazard over (tsart, tstop)
#' @export
#'
#' @examples
#' h12 <- Hazard(formula = 0~1, statefrom = 1, stateto=2, family="wei")
#' h13 <- Hazard(formula = 0~1, statefrom = 1, stateto=3, family="wei")
#' h23 <- Hazard(formula = 0~1, statefrom = 2, stateto=3, family="wei")
#' model_fit <- multistatemodel(hazard = c(h12, h13, h23), data=illness_death_dat)
#' compute_cumulative_hazard(tstart = 0, tstop =0.5, model=model_fit, hazard = "h12")
#' \dontshow{
#' JuliaConnectoR::stopJulia()
#' }
compute_cumulative_hazard <- function(tstart, tstop, model, hazard, subj = 1) {

  if (length(tstart) != length(tstop)) {
    stop("tstart and tstop must be the same length")
  }
  hazards <- NULL

  hazard <- paste0(":", hazard)

  for (i in 1:length(tstart)) {
    hazards <- append(hazards, multistatemodels_env$compute_cumulative_hazard(tstart[i], tstop[i], model, JuliaConnectoR::juliaEval(hazard), JuliaConnectoR::juliaCall("Int", subj))[1])
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
#' @param min_ess minimum effective sample size, defaults to 100.
#'
#' @return A Julia environment with a named Tuple
#' @export
#'
#' @examples
#' h12 <- Hazard(formula = 0~1, statefrom = 1, stateto=2, family="wei")
#' h13 <- Hazard(formula = 0~1, statefrom = 1, stateto=3, family="wei")
#' h23 <- Hazard(formula = 0~1, statefrom = 2, stateto=3, family="wei")
#' model_fit <- multistatemodel(hazard = c(h12, h13, h23), data=illness_death_dat)
#' estimate_loglik(model = model_fit)
#' \dontshow{
#' JuliaConnectoR::stopJulia()
#' }
estimate_loglik <- function(model, min_ess = 100) {

  multistatemodels_env$estimate_loglik(model, min_ess = JuliaConnectoR::juliaCall("Int", min_ess))

}

#' Draw paths
#'
#' @description
#' Draw sample paths conditional on the data.
#' Require that the minimum effective sample size is greater than min_ess.
#'
#'
#' @param model multistate model
#' @param min_ess minimum effective sample size, defaults to 100.
#' @param paretosmooth logical; pareto smooth importance weights, defaults to TRUE unless min_ess < 25.
#' @param return_logliks logical; defaults to FALSE
#'
#' @return A Julia environment with a named Tuple of sample paths
#' @export
#'
#' @examples
#' h12 <- Hazard(formula = 0~1, statefrom = 1, stateto=2, family="wei")
#' h13 <- Hazard(formula = 0~1, statefrom = 1, stateto=3, family="wei")
#' h23 <- Hazard(formula = 0~1, statefrom = 2, stateto=3, family="wei")
#' model_fit <- multistatemodel(hazard = c(h12, h13, h23), data=illness_death_dat)
#' draw_paths(model = model_fit)
#' \dontshow{
#' JuliaConnectoR::stopJulia()
#' }
draw_paths <- function(model, min_ess = 100, paretosmooth=TRUE, return_logliks = FALSE) {

  multistatemodels_env$draw_paths(model, min_ess = JuliaConnectoR::juliaCall("Int", min_ess), paretosmooth = paretosmooth, return_logliks = return_logliks)

}


#' Cumulative incidence
#'
#' @description
#' Compute the cumulative incidence for each possible transition as a function of time since state entry.
#' Assumes the subject starts their observation period at risk and saves cumulative incidence at the supplied vector of times, t.
#'
#' @param t time or a vector of times
#' @param model multistate model
#' @param subj subject id
#'
#' @return A matrix with the cumulative incidence where columns represent each
#' possible transition and rows represent t
#' @export
#'
#' @examples
#' h12 <- Hazard(formula = 0~1, statefrom = 1, stateto=2, family="wei")
#' h13 <- Hazard(formula = 0~1, statefrom = 1, stateto=3, family="wei")
#' h23 <- Hazard(formula = 0~1, statefrom = 2, stateto=3, family="wei")
#' model_fit <- multistatemodel(hazard = c(h12, h13, h23), data=illness_death_dat)
#' cumulative_incidence(t=0.5, model=model_fit)
#' \dontshow{
#' JuliaConnectoR::stopJulia()
#' }
cumulative_incidence <- function(t, model, subj=1) {

  multistatemodels_env$cumulative_incidence(t, model, JuliaConnectoR::juliaCall("Int", subj))

}


#' Set parameters
#'
#' @description Set model parameters given a vector of values.
#' Copies newvalues to model.parameters.
#'
#' @param model multistate model
#' @param newvalues A list, each element of the list is new parameters for each hazard
#'
#' @export
#'
#' @examples
#' h12 <- Hazard(formula = 0~1, statefrom = 1, stateto=2, family="wei")
#' h13 <- Hazard(formula = 0~1, statefrom = 1, stateto=3, family="wei")
#' h23 <- Hazard(formula = 0~1, statefrom = 2, stateto=3, family="wei")
#' model_fit <- multistatemodel(hazard = c(h12, h13, h23), data=illness_death_dat)
#' model_fit <- set_parameters(model = model_fit, newvalues = list(h12 = c(log(1.25), log(1.5)),
#' h13 = c(log(1.25), log(1)), h23 = c(log(1.25), log(2))))
#' \dontshow{
#' JuliaConnectoR::stopJulia()
#' }
set_parameters <- function(model, newvalues) {

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

  #mod <- rlang::duplicate(model)

  multistatemodels_env$set_parameters(model, JuliaConnectoR::juliaCall("Tuple", JuliaConnectoR::juliaEval(values)))
  #return(mod)

}


#' Simulate data
#'
#' @description
#' Simulate n data sets or collections of sample paths from a multistate model.
#' If data = TRUE (the default) discretely observed sample paths are returned, possibly subject to measurement error.
#' If paths = FALSE (the default), continuous-time sample paths are not returned.
#'
#' @param model multistate model
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
#' h12 <- Hazard(formula = 0~1, statefrom = 1, stateto=2, family="wei")
#' h13 <- Hazard(formula = 0~1, statefrom = 1, stateto=3, family="wei")
#' h23 <- Hazard(formula = 0~1, statefrom = 2, stateto=3, family="wei")
#' model_fit <- multistatemodel(hazard = c(h12, h13, h23), data=illness_death_dat)
#' d <- simulate(model = model_fit)
#' \dontshow{
#' JuliaConnectoR::stopJulia()
#' }
simulate <- function(model, nsim=1, data=TRUE, paths=FALSE, delta_u = sqrt(.Machine$double.eps), delta_t = sqrt(.Machine$double.eps)) {

  sim_dat <- multistatemodels_env$simulate(model, nsim = JuliaConnectoR::juliaCall("Int", nsim), data=data, paths=paths, delta_u=delta_u, delta_t=delta_t)
  if (isFALSE(paths)) {
    as.data.frame(sim_dat[[1]])
  } else {
    JuliaConnectoR::juliaGet(sim_dat)
  }

}


#' Summary of model output
#'
#' @param model_fitted fitted multistate models object
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
#' h12 <- Hazard(formula = 0~1, statefrom = 1, stateto=2, family="wei")
#' h13 <- Hazard(formula = 0~1, statefrom = 1, stateto=3, family="wei")
#' h23 <- Hazard(formula = 0~1, statefrom = 2, stateto=3, family="wei")
#' model_fit <- multistatemodel(hazard = c(h12, h13, h23), data=illness_death_dat)
#' model_fitted <- fit(model = model_fit, verbose=TRUE, compute_vcov = TRUE, ess_target_initial = 50,
#' ascent_threshold = 0.2, stopping_threshold = 0.2, tol = 0.001)
#' model_summary(model_fitted = model_fitted, estimate_likelihood = TRUE)
#' \dontshow{
#' JuliaConnectoR::stopJulia()
#' }
model_summary <- function(model_fitted, confidence_level = 0.95, estimate_likelihood=FALSE, min_ess = 100) {

  multistatemodels_env$summary(model_fitted, confidence_level = confidence_level, estimate_likelihood=estimate_likelihood, min_ess = JuliaConnectoR::juliaCall("Int", min_ess))

}





