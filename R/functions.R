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




check_missing <- function(df, ...) {
  numMissing <- sum(apply(df, 2, function(x) {sum(is.na(x))}))
  if (numMissing == 0) {
    FALSE
  } else {
    TRUE
  }
}

#' Collapse subjects and recompute sampling weights
#'
#' @description Collapse subjects to create an internal representation of a
#' dataset and optionally recompute a vector of sampling weights.
#'
#' @param df A data frame
#' @param multistatemodels Loaded multistatemodels Julia environment
#' @param SamplingWeights Sampling weights
#'
#' @return A 2-element list. The first element is a data frame and the second
#' is a vector of sampling weights
#' @export
#'
#' @examples
#' multistatemodels <- multistatemodels_load()
#' df <- read.csv("~/derived_states_regen.csv")
#' d <- collapse_data(df[ ,-12], multistatemodels)
#' \dontshow{
#' JuliaConnectoR::stopJulia()
#' }
collapse_data <- function(df, multistatemodels, SamplingWeights = rep(1, length(unique(df$id)))) {
  if (check_missing(df)) {
    stop("No missing data allowed in data frame")
  }

  d <- multistatemodels$collapse_data(JuliaConnectoR::juliaCall("DataFrame", df), SamplingWeights = SamplingWeights)
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
#' @param family One of "exp" or "wei" for exponential or Weibull cause-specific
#' baseline hazard functions, or "sp" for a semi-parametric spline basis for the
#' baseline hazard (defaults to M-splines).
#' @param degree For family="sp": degree of the spline polynomial basis
#' @param knots For family="sp": Vector of knots
#' @param ... ADD HERE
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
Hazard <- function(formula, statefrom, stateto, multistatemodels, family = c("exp", "wei", "sp"), degree=3, knots=NULL, ...) {
  if (family == "sp" & (is.null(degree) | is.null(knots))) {
    stop("Degree and knots must be provided for splines")
  }

  if (family != "sp" & (!is.null(knots))) {
    warning("Degree and knots not required for this model, will be ignored")
  }

  if (!family %in% c("exp", "wei", "sp")) {
    stop("family must be one of 'exp', 'wei', or 'sp'")
  }

  rhs <- paste(c(1, labels(stats::terms(formula))), collapse = " + ")
  form <- JuliaConnectoR::juliaEval(paste0("@formula(", as.character(stats::terms(formula)[[2]]), as.character(stats::terms(formula)[[1]]), rhs, ")"))

  if (family == "sp") {
    knots <- JuliaConnectoR::juliaEval(paste0("vec([", knots, "])"))
    haz <- multistatemodels$Hazard(form, family, JuliaConnectoR::juliaCall("Int", statefrom), JuliaConnectoR::juliaCall("Int", stateto),
                                   degree = JuliaConnectoR::juliaCall("Int", degree), knots = knots)
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
#' @param verbose Defaults to FALSE
#'
#' @return A Julia MultistateModel that can be used for simulation and inference
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


  htuple <- JuliaConnectoR::juliaCall("Tuple", hazard)
  out <- "multistatemodel(htuple..., data=data,SamplingWeights=SamplingWeights, CensoringPatterns=CensoringPatterns, verbose=verbose)"
  mod <- JuliaConnectoR::juliaLet(out, htuple=htuple, data=JuliaConnectoR::juliaCall("DataFrame", data), SamplingWeights = SamplingWeights,
                                  CensoringPatterns = matrix(as.integer(CensoringPatterns), nrow=nPatterns),
                                  verbose=verbose)
  #JuliaConnectoR::juliaGet(mod)
  return(mod)
}
