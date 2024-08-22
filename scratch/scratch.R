library(JuliaConnectoR)
library(readr)

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
    multistatemodels <- JuliaConnectoR::juliaImport("MultistateModels", all = FALSE)
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



collapse_data <- function(dat, multistatemodels, SamplingWeights = rep(1, length(unique(dat$id)))) {



  df <- readr::format_csv(dat)

  df2 <- juliaLet("CSV.read(IOBuffer(df), DataFrame)", df=df)
  df2 <- JuliaConnectoR::juliaCall("DataFrame", df2)

  check_data(df2)

  d <- multistatemodels$collapse_data(df2, SamplingWeights = SamplingWeights)
  return(list(as.data.frame(d[[1]]), d[[2]]))


}


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

initialize_parameters <- function(model, multistatemodels, constraints = NULL, surrogate_constraints = NULL, surrogate_parameters=NULL, crude=FALSE) {

  multistatemodels$`initialize_parameters!`(model, constraints = constraints, surrogate_constraints = surrogate_constraints, surrogate_parameters=surrogate_parameters, crude=crude)

}


fit_surrogate <- function(model, multistatemodels, surrogate_parameters = NULL, surrogate_constraints = NULL, crude_inits = T, verbose = T) {

  multistatemodels$fit_surrogate(model, surrogate_parameters = surrogate_parameters, surrogate_constraints = surrogate_constraints, crude_inits = crude_inits, verbose = verbose)

}


fit <- function(model, multistatemodels, verbose=T, compute_vcov=T, constraints = NULL,
                optimize_surrogate=T, surrogate_constraints = NULL, surrogate_parameters = NULL, maxiter=200, tol=1e-3,
                alpha=0.05, gamma=0.05, kappa=4/3, ess_target_initial = 100, MaxSamplingEffort=20, npaths_additional = 10,
                return_ConvergenceRecords = T, return_ProposedPaths = F) {

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



aic <- function(model_fitted, multistatemodels, estimate_likelihood = TRUE, min_ess=100, loglik=NULL) {

  if (!JuliaConnectoR::juliaLet("isa(model_fitted, MultistateModels.MultistateModelFitted)", model_fitted = model_fitted)) {
    stop("Fitted model of type MultistateModels.MultistateModelFitted must be provided")
  }

  multistatemodels$aic(model_fitted, estimate_likelihood = estimate_likelihood, min_ess = JuliaConnectoR::juliaCall("Int",min_ess), loglik=loglik)

}

bic <- function(model_fitted, multistatemodels, estimate_likelihood = T, min_ess=100, loglik=NULL) {

  if (!JuliaConnectoR::juliaLet("isa(model_fitted, MultistateModels.MultistateModelFitted)", model_fitted = model_fitted)) {
    stop("Fitted model of type MultistateModels.MultistateModelFitted must be provided")
  }

  multistatemodels$bic(model_fitted, estimate_likelihood = estimate_likelihood, min_ess = JuliaConnectoR::juliaCall("Int",min_ess), loglik=loglik)

}

get_ConvergenceRecords <- function(model_fitted, multistatemodels) {

  if (!JuliaConnectoR::juliaLet("isa(model_fitted, MultistateModels.MultistateModelFitted)", model_fitted = model_fitted)) {
    stop("Fitted model of type MultistateModels.MultistateModelFitted must be provided")
  }

  multistatemodels$get_ConvergenceRecords(model_fitted)

}

get_loglik <- function(model_fitted, multistatemodels, ll="loglik") {

  if (!JuliaConnectoR::juliaLet("isa(model_fitted, MultistateModels.MultistateModelFitted)", model_fitted = model_fitted)) {
    stop("Fitted model of type MultistateModels.MultistateModelFitted must be provided")
  }

  if (!ll %in% c("loglik", "subj_lml")){
    stop("ll must be either 'loglik' or 'subj_lml'")
  }

  multistatemodels$get_loglik(model_fitted, ll=ll)

}

get_parameters <- function(model_fitted, multistatemodels) {

  if (!JuliaConnectoR::juliaLet("isa(model_fitted, MultistateModels.MultistateModelFitted)", model_fitted = model_fitted)) {
    stop("Fitted model of type MultistateModels.MultistateModelFitted must be provided")
  }

  multistatemodels$get_parameters(model_fitted)

}

get_parnames <- function(model_fitted, multistatemodels) {

  if (!JuliaConnectoR::juliaLet("isa(model_fitted, MultistateModels.MultistateModelFitted)", model_fitted = model_fitted)) {
    stop("Fitted model of type MultistateModels.MultistateModelFitted must be provided")
  }

  multistatemodels$get_parnames(model_fitted)

}

get_vcov <- function(model_fitted, multistatemodels) {

  if (!JuliaConnectoR::juliaLet("isa(model_fitted, MultistateModels.MultistateModelFitted)", model_fitted = model_fitted)) {
    stop("Fitted model of type MultistateModels.MultistateModelFitted must be provided")
  }

  multistatemodels$get_vcov(model_fitted)

}

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


compute_hazard <- function(t, model, hazard, multistatemodels, subj=1) {

  hazards <- NULL

  hazard <- paste0(":", hazard)

  for (i in 1:length(t)) {
    hazards <- append(hazards, multistatemodels$compute_hazard(t[i], model, JuliaConnectoR::juliaEval(hazard), JuliaConnectoR::juliaCall("Int", subj))[1])
  }

  return(hazards)

}

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


estimate_loglik <- function(model, multistatemodels, min_ess = 100) {

  multistatemodels$estimate_loglik(model, min_ess = JuliaConnectoR::juliaCall("Int", min_ess))

}

draw_paths <- function(model, multistatemodels, min_ess = 100, paretosmooth=TRUE, return_logliks = FALSE) {

  multistatemodels$draw_paths(model, min_ess = JuliaConnectoR::juliaCall("Int", min_ess), paretosmooth = paretosmooth, return_logliks = return_logliks)

}


cumulative_incidence <- function(t, model, multistatemodels, subj=1) {

  multistatemodels$cumulative_incidence(t, model, JuliaConnectoR::juliaCall("Int", subj))

}


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


simulate <- function(model, multistatemodels, nsim=1, data=T, paths=F, delta_u = sqrt(.Machine$double.eps), delta_t = sqrt(.Machine$double.eps)) {

  multistatemodels$simulate(model, nsim = JuliaConnectoR::juliaCall("Int", nsim), data=data, paths=paths, delta_u=delta_u, delta_t=delta_t)

}


model_summary <- function(model_fitted, multistatemodels, confidence_level = 0.95, estimate_likelihood=F, min_ess = 100) {

  multistatemodels$summary(model_fitted, confidence_level = confidence_level, estimate_likelihood=estimate_likelihood, min_ess = JuliaConnectoR::juliaCall("Int", min_ess))

}



#multistatemodels_setup()
multistatemodels <- multistatemodels_load()
df <- read.csv("H:/Projects/Jon/R package and Julia/Julia/derived_states_regen.csv")
d <- collapse_data(df, multistatemodels)

h12 <- Hazard(formula = 0 ~ 1+mab, family = "exp", statefrom = 1, stateto=2, multistatemodels=multistatemodels)
#h12 <- Hazard(formula = 0 ~ 1+mab, family = "sp", statefrom = 1, stateto=2, multistatemodels=multistatemodels, degree=1, knots=5/7)
h23 <- Hazard(formula = 0 ~ 1+mab, family = "exp", statefrom = 2, stateto=3, multistatemodels=multistatemodels)
h24 <- Hazard(formula = 0 ~ 1+mab, family = "exp", statefrom = 2, stateto=4, multistatemodels=multistatemodels)
h45 <- Hazard(formula = 0 ~ 1+mab, family = "exp", statefrom = 4, stateto=5, multistatemodels=multistatemodels)
#h45 <- Hazard(formula = 0 ~ 1+mab, family = "exp", statefrom = 4, stateto=5, multistatemodels=multistatemodels, degree=1)

CensoringPatterns <- matrix(c(3, 1, 0, 1, 0, 0), nrow=1, ncol=6)
#CensoringPatterns <- matrix(c(3, 1, 1, 1, 0), nrow=1, ncol=5) #For cav data. For pseudo-observation before directly observed death
#First is the obstype. Then are the 4 states.
#Can go between any states before death

model <- multistatemodel(hazard = c(h12, h23, h24, h45), data = d[[1]], SamplingWeights = d[[2]],
                         CensoringPatterns = CensoringPatterns, nPatterns = 1, multistatemodels =multistatemodels)

initialize_parameters(model, multistatemodels = multistatemodels)

model_fitted <- fit(model, multistatemodels = multistatemodels, tol=1e-3, ess_target_initial=500)


m <- juliaGet(model)

collapse_data_v2 <- JuliaConnectoR::juliaFun("collapse_data")
dat <- readr::format_csv(df)

dat2 <- juliaLet("CSV.read(IOBuffer(df), DataFrame)", df=dat)
dat2 <- JuliaConnectoR::juliaCall("DataFrame", dat2)

d <- collapse_data_v2(dat2)
d1 <- as.data.frame(d[[1]])
d2 <- d[[2]]




functions <- names(multistatemodels)

for (i in 1:length(functions)) {

  assign(functions[i], juliaFun(functions[i]))

}


i <- 7
assign(gsub("!", "", functions[i]), juliaFun(functions[i]))


assign("model_summary", juliaFun("summary"))




hlist <- c(h12, h23, h24, h45)
htuple <- juliaCall("Tuple", hlist)

hazinfo <- juliaLet("MultistateModels.enumerate_hazards(htuple...)", htuple=htuple)
tmat <- juliaLet("MultistateModels.create_tmat(hazinfo)", hazinfo=hazinfo)
juliaLet("MultistateModels.check_data!(data, tmat, CensoringPatterns, verbose=true)", tmat=tmat, data=JuliaConnectoR::juliaCall("DataFrame", df),
         CensoringPatterns = matrix(as.integer(CensoringPatterns), nrow=1))


