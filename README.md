
<!-- README.md is generated from README.Rmd. Please edit that file -->

# multistatemodels

<!-- badges: start -->
<!-- badges: end -->

multistatemodels provides a package for simulating and fitting
multistate models especially semi-Markov models, to panel data and
interval censored data. This is a wrapper to the Julia package
MultistateModels.jl. Add some more about why it is good to call Julia
from R…advantages, etc.

## Installation

You can install the development version of multistatemodels from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ammateja/multistatemodels")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(multistatemodels)
#> Starting Julia ...
#> 
#> Attaching package: 'multistatemodels'
#> The following object is masked from 'package:stats':
#> 
#>     simulate
```

``` r
library(JuliaConnectoR)
#> Warning: package 'JuliaConnectoR' was built under R version 4.3.3
```

``` r
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
```

``` r
library(binom)
#> Warning: package 'binom' was built under R version 4.3.3
```

``` r
library(survival)
#> Warning: package 'survival' was built under R version 4.3.3
```

``` r

#Function to make parameters
makepars <- function() {
  parameters <- list(h12 = c(log(1.25), log(1.5)), 
                     h13 = c(log(1.25), log(1)), 
                     h23 = c(log(1.25), log(2)))
  return(parameters)
}

#Function to make the assessment times
make_obstimes <- function(ntimes) {
  #observation times
  interval <- seq(0, 1, length.out = ntimes+1)
  times <- interval[-c(1, length(interval))] + (rbeta(length(interval[-c(1, length(interval))]), 1.5, 1.5) - 0.5)*diff(interval)[1]
  times <- c(0, times, 1)

  return(times)
}

#Function to set up the model
setup_model <- function(make_pars, data = NULL, nsubj = 300, family = "wei", ntimes = 6, spknots = NULL, sim_exp = FALSE) {
  
  if (family != "sp1" & family != "sp2") {
    
    h12 <- multistatemodels::Hazard(formula = 0~1, statefrom = 1, stateto=2, family=family)
    h13 <- multistatemodels::Hazard(formula = 0~1, statefrom = 1, stateto=3, family=family)
    h23 <- multistatemodels::Hazard(formula = 0~1, statefrom = 2, stateto=3, family=family)
    
  } else if (family == "sp1") {
      
    knots12 <- spknots[[1]]
    knots13 <- spknots[[2]]
    knots23 <- spknots[[3]]
    
    h12 <- multistatemodels::Hazard(formula = 0~1, statefrom = 1, stateto=2, family="sp", degree=1, knots = knots12[-c(1, length(knots12))], boundaryknots = knots12[c(1, length(knots12))], extrapolation = "flat")
    h13 <- multistatemodels::Hazard(formula = 0~1, statefrom = 1, stateto=3, family="sp", degree=1, knots = knots13[-c(1, length(knots13))], boundaryknots = knots13[c(1, length(knots13))], extrapolation = "flat")
    h23 <- multistatemodels::Hazard(formula = 0~1, statefrom = 2, stateto=3, family="sp", degree=1, knots = knots23[-c(1, length(knots23))], boundaryknots = knots23[c(1, length(knots23))], extrapolation = "flat")
    
  } else if (family == "sp2") {
    
    knots12 <- spknots[[1]]
    knots13 <- spknots[[2]]
    knots23 <- spknots[[3]]
    
    h12 <- multistatemodels::Hazard(formula = 0~1, statefrom = 1, stateto=2, family="sp", degree=3, knots = knots12[-c(1, length(knots12))], boundaryknots = knots12[c(1, length(knots12))], extrapolation = "flat", monotone=0)
    h13 <- multistatemodels::Hazard(formula = 0~1, statefrom = 1, stateto=3, family="sp", degree=3, knots = knots13[-c(1, length(knots13))], boundaryknots = knots13[c(1, length(knots13))], extrapolation = "flat", monotone=0)
    h23 <- multistatemodels::Hazard(formula = 0~1, statefrom = 2, stateto=3, family="sp", degree=3, knots = knots23[-c(1, length(knots23))], boundaryknots = knots23[c(1, length(knots23))], extrapolation = "flat")
    
  }
  
  if (is.null(data)) {
    
    for (i in 1:nsubj) {
      if(ntimes == 12) {
        visitdays <- make_obstimes(12)
      } else if (ntimes == 6) {
        visitdays <- make_obstimes(12)[c(1, 3, 5, 7, 9, 11, 13)]
      } else {
        visitdays <- make_obstimes(12)[c(1, 4, 7, 10, 13)]
      }
      id <- rep(i, ntimes)
      tstart <- visitdays[-length(visitdays)]
      tstop <- visitdays[-1]
      statefrom <- rep(1, ntimes)
      stateto <- rep(1, ntimes)
      obstype <- rep(2, ntimes)
      d <- data.frame(id=id, tstart=tstart, tstop=tstop, statefrom=statefrom, stateto=stateto, obstype=obstype)
      data <- rbind(data, d)
    }
    
  }
  
  model <- multistatemodels::multistatemodel(hazard = c(h12, h13, h23), data=data)
  
  if (make_pars & !sim_exp) {
    
    parameters <- makepars()
    model <- multistatemodels::set_parameters(model=model, newvalues=parameters)
    
  } else if (make_pars & sim_exp) {
    
    model <- multistatemodels::set_parameters(model=model, newvalues=list(h12 = log(1.5), 
                                                                          h13 = log(1), 
                                                                          h23 = log(2)))
    
  }
  
  return(model)
  
}


observe_subjdat <- function(path, model) {
  
  d <- as.data.frame(model$data)
  subjdat <- NULL
  
  for (i in 1:length(path)) {
    
    subj.dat.raw <- d[d$id == i, ]
    obstimes <- unique(sort(c(0, subj.dat.raw$tstop, path[[i]]$times[path[[i]]$states ==3])))
    
    if (sum(c(2,3) %in% path[[i]]$states) == 2) {
      
      new <- path[[i]]$times[length(path[[i]]$times)] - sqrt(.Machine$double.eps)
      obstimes <- sort(c(obstimes, new))
      
    }
    
    obsinds <- unlist(lapply(obstimes, function(x){max(which(path[[i]]$times <= x))}))
    obsstates <- path[[i]]$states[obsinds]
    
    data <- data.frame(id = rep(path[[i]]$subj, (length(obstimes)-1)), 
                      tstart = obstimes[-length(obstimes)], 
                      tstop = obstimes[-1], 
                      statefrom = obsstates[-length(obsstates)], 
                      stateto = obsstates[-1])
    data <- data[data$stateto != 3 | data$statefrom != 3, ]
    
    data$obstype <- ifelse(data$stateto == 3, 1, 2)
    data$obstype <- ifelse(data$statefrom == data$stateto, 1, data$obstype)
    
    subjdat <- rbind(subjdat, data)
    
  }
  
  subjdat <- subjdat %>% dplyr::arrange(id, tstart)
  
  return(subjdat)
  
}


get_estimates <- function(paths) {
  # progression-free survival
  # percent who progress
  # percent who die after progression
  # percent who die without progressing
  pfs <- mean(unlist(lapply(paths, function(x){all(x$states == 1)})), na.rm=T)
  prog <- mean(unlist(lapply(paths, function(x){2 %in% x$states})), na.rm=T)
  die_wprog <- mean(unlist(lapply(paths, function(x){sum(c(2,3) %in% x$states) == 2})), na.rm=T)
  die_noprog <- mean(unlist(lapply(paths, function(x){3 %in% x$states & !(2 %in% x$states)})), na.rm=T)
  
  # restricted mean progression-free survival time
  rmpfst <- mean(unlist(lapply(paths, function(x){x$times[2]})), na.rm=T)
  time2prog_all <- mean(unlist(lapply(paths, function(x){ifelse(x$states[2] == 2, x$times[2], 1)})), na.rm=T)
  
  # time to disease among the progressors
  proginds <- unlist(lapply(paths, function(x){2 %in% x$states & length(x$states) != 2}))
  time2prog_gprog <- sum(unlist(lapply(paths[proginds], function(x){x$times[2]})), na.rm=T)/sum(proginds, na.rm=T)
  illnessdur <- sum(unlist(lapply(paths[proginds], function(x){x$times[3] - x$times[2]})), na.rm=T)/sum(proginds, na.rm=T)
  
  ests <- data.frame(pfs=pfs, 
                     prog = prog, 
                     die_wprog = die_wprog, 
                     die_noprog = die_noprog, 
                     rmpfst = rmpfst, 
                     time2prog_gprog = time2prog_gprog, 
                     time2prog_all = time2prog_all, 
                     illnessdur = illnessdur)
  
  return(ests)
  
}


asymptotic_bootstrap <- function(model, pars, vcov, sims_per_subj, nboot) {
  
  npars <- length(pars)
  pardraws <- rep(0, npars)
  
  U <- svd(vcov)$u
  D <- svd(vcov)$d
  
  D[D < 0] <- 0
  
  S <- U%*%diag(sqrt(D))
  
  ests <- matrix(0, nrow=8, ncol=nboot)
  
  for (k in 1:nboot) {
    
    pardraws[1:npars] <- as.matrix(pars) + S%*%rnorm(npars)
    
    elem_ptr <- JuliaConnectoR::juliaGet(model$parameters)$elem_ptr
    newvalues <- vector("list", length(elem_ptr)-1)
    for (i in 1:(length(elem_ptr)-1)) {
      newvalues[[i]] <- pardraws[elem_ptr[i]:(elem_ptr[i+1]-1)]
    }
    
    model <- multistatemodels::set_parameters(model=model, newvalues = newvalues)
    paths_sim <- multistatemodels::simulate(model=model, nsim = sims_per_subj, paths=TRUE, data=FALSE)
    
    ests[ ,k] <- as.numeric(get_estimates(paths_sim))
    
  }
  
  
  return(t(apply(ests, 1, quantile, na.rm=T, probs = c(0.025, 0.975))))
  
}


# wrapper for one simulation
# nsim is the number of simulated paths per subject
# ndraws is the number of draws from the asymptotic normal distribution of the MLEs
# family:
#   1: "exp"
#   2: "wei"
#   3: "sp1" - degree 1 with knot at midpoint and range
#   4: "sp2" - degree 3 with knots at 0.05, 1/3 and 2/3, and 0.95 quantiles
work_function <- function(simnum, seed, family, ntimes, sims_per_subj, nboot, sim_exp=FALSE) {
  
  set.seed(seed)
  
  model_sim <- setup_model(make_pars=TRUE, data = NULL, family = "wei", ntimes = ntimes, nsubj = 250, sim_exp=FALSE)
  paths <- multistatemodels::simulate(model=model_sim, nsim=1, paths=TRUE, data=FALSE)
  dat <- observe_subjdat(paths, model_sim)
  
  if (family < 3) {
    spknots <- NULL
  } else if (family == 3) {
    
    q12 <- c(0, quantile(JuliaConnectoR::juliaLet("MultistateModels.extract_sojourns(1, 2, MultistateModels.extract_paths(dat))", dat=JuliaConnectoR::juliaCall("DataFrame", dat)), c(0.5, 1)))
    q13 <- c(0, quantile(JuliaConnectoR::juliaLet("MultistateModels.extract_sojourns(1, 3, MultistateModels.extract_paths(dat))", dat=JuliaConnectoR::juliaCall("DataFrame", dat)), c(0.5, 1)))
    q23 <- c(0, quantile(JuliaConnectoR::juliaLet("MultistateModels.extract_sojourns(2, 3, MultistateModels.extract_paths(dat))", dat=JuliaConnectoR::juliaCall("DataFrame", dat)), c(0.5, 1)))
    
    spknots <- list(q12, q13, q23)
    
  } else if (family == 4) {
    
    q12 <- c(0, quantile(JuliaConnectoR::juliaLet("MultistateModels.extract_sojourns(1, 2, MultistateModels.extract_paths(dat))", dat=JuliaConnectoR::juliaCall("DataFrame", dat)), c(1/3, 2/3, 1)))
    q13 <- c(0, quantile(JuliaConnectoR::juliaLet("MultistateModels.extract_sojourns(1, 3, MultistateModels.extract_paths(dat))", dat=JuliaConnectoR::juliaCall("DataFrame", dat)), c(1/3, 2/3, 1)))
    q23 <- c(0, quantile(JuliaConnectoR::juliaLet("MultistateModels.extract_sojourns(2, 3, MultistateModels.extract_paths(dat))", dat=JuliaConnectoR::juliaCall("DataFrame", dat)), c(1/3, 2/3, 1)))
    
    spknots <- list(q12, q13, q23)
    
  }
  
  model_fit <- setup_model(make_pars = FALSE, data = dat, family = c("exp", "wei", "sp1", "sp2")[family], ntimes=ntimes, spknots = spknots, sim_exp = sim_exp)
  model_fit <- multistatemodels::initialize_parameters(model = model_fit)
  model_fitted <- multistatemodels::fit(model = model_fit, verbose=TRUE, compute_vcov = TRUE, ess_target_initial = 50, ascent_threshold = 0.2, stopping_threshold = 0.2, tol = 0.001)
  
  model_sim2 <- setup_model(make_pars = FALSE, data = as.data.frame(model_sim$data), family = c("exp", "wei", "sp1", "sp2")[family], spknots = spknots, ntimes=ntimes, sim_exp=sim_exp)
  model_sim2 <- multistatemodels::set_parameters(model = model_sim2, newvalues = list(JuliaConnectoR::juliaGet(model_fitted$parameters[1])$data, 
                                                                                      JuliaConnectoR::juliaGet(model_fitted$parameters[2])$data, 
                                                                                      JuliaConnectoR::juliaGet(model_fitted$parameters[3])$data))
  paths_sim <- multistatemodels::simulate(model = model_sim2, nsim = sims_per_subj, paths = TRUE, data = FALSE)
  
  ests <- get_estimates(paths_sim)
  
  asymp_cis <- asymptotic_bootstrap(model_sim2, JuliaConnectoR::juliaGet(model_fitted$parameters)$data, model_fitted$vcov, sims_per_subj, nboot)
  
  results <- data.frame(simnum = rep(simnum, 8), 
                        family = rep(family, 8), 
                        ntimes = rep(ntimes, 8), 
                        var = names(ests), 
                        ests = as.numeric(ests), 
                        lower = asymp_cis[, 1], 
                        upper = asymp_cis[ ,2])
  
  return(results)
  
}




# sim_res_exp1 <- work_function(1, 1, 1, 6, 20, 100, TRUE)
# sim_res_exp2 <- work_function(1, 1, 1, 12, 20, 100, TRUE)
# sim_res_exp3 <- work_function(1, 1, 1, 4, 20, 100, TRUE)
# 
# sim_res_wei1 <- work_function(1, 1, 2, 6, 20, 100)
# sim_res_wei2 <- work_function(1, 1, 2, 12, 20, 100)
# sim_res_wei3 <- work_function(1, 1, 2, 4, 20, 100)
# 
# sim_res_sp11 <- work_function(1, 1, 3, 6, 20, 100)
# sim_res_sp12 <- work_function(1, 1, 3, 12, 20, 100)
# sim_res_sp13 <- work_function(1, 1, 3, 4, 20, 100)
# 
# sim_res_sp21 <- work_function(1, 1, 4, 6, 20, 100)
# sim_res_sp22 <- work_function(1, 1, 4, 12, 20, 100)
# sim_res_sp23 <- work_function(1, 1, 4, 4, 20, 100)
```
