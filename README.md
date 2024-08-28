
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
set.seed(1)

#Function to make parameters
makepars <- function() {
  parameters <- list(h12 = c(log(1.25), log(1)), 
                     h13 = c(log(0.8), log(1)), 
                     h23 = c(log(1), log(1.25)))
  return(parameters)
}

#Function to make the assessment times
make_obstimes <- function() {
  #observation times
  interval <- seq(0.1, 0.9, 0.1)
  times <- interval + (rbeta(length(interval), 1.5, 1.5) - 0.5)*0.1
  times <- c(0, times, 1)

  return(times)
}

#Function to set up the model
setup_model <- function(make_pars, data=NULL, nsubj=300, family="wei", ntimes=10, spknots=NULL) {
  
  if (family != "sp1" & family != "sp2") {
    
    h12 <- multistatemodels::Hazard(formula = 0~1, statefrom = 1, stateto=2, family=family)
    h13 <- multistatemodels::Hazard(formula = 0~1, statefrom = 1, stateto=3, family=family)
    h23 <- multistatemodels::Hazard(formula = 0~1, statefrom = 2, stateto=3, family=family)
    
  } else if (family == "sp1") {
      
    knots12 <- spknots[1]
    knots13 <- spknots[2]
    knots23 <- spknots[3]
    
    h12 <- multistatemodels::Hazard(formula = 0~1, statefrom = 1, stateto=2, family="sp", degree=1, knots = knots12[-c(1, length(knots12))], boundaryknots = knots12[c(1, length(knots12))], extrapolation = "flat")
    h13 <- multistatemodels::Hazard(formula = 0~1, statefrom = 1, stateto=3, family="sp", degree=1, knots = knots13[-c(1, length(knots13))], boundaryknots = knots13[c(1, length(knots13))], extrapolation = "flat")
    h23 <- multistatemodels::Hazard(formula = 0~1, statefrom = 2, stateto=3, family="sp", degree=1, knots = knots23[-c(1, length(knots23))], boundaryknots = knots23[c(1, length(knots23))], extrapolation = "flat")
    
  } else if (family == "sp2") {
    
    knots12 <- spknots[1]
    knots13 <- spknots[2]
    knots23 <- spknots[3]
    
    h12 <- multistatemodels::Hazard(formula = 0~1, statefrom = 1, stateto=2, family="sp", degree=3, knots = knots12[-c(1, length(knots12))], boundaryknots = knots12[c(1, length(knots12))], extrapolation = "flat", monotone=0)
    h13 <- multistatemodels::Hazard(formula = 0~1, statefrom = 1, stateto=3, family="sp", degree=3, knots = knots13[-c(1, length(knots13))], boundaryknots = knots13[c(1, length(knots13))], extrapolation = "flat", monotone=0)
    h23 <- multistatemodels::Hazard(formula = 0~1, statefrom = 2, stateto=3, family="sp", degree=1, knots = knots23[-c(1, length(knots23))], boundaryknots = knots23[c(1, length(knots23))], extrapolation = "flat")
    
  }
  
  if (is.null(data)) {
    
    for (i in 1:nsubj) {
      visitdays <- make_obstimes()
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
  
  if (make_pars) {
    
    parameters <- makepars()
    multistatemodels::set_parameters(model=model, parameters=parameters)
    
  }
  
  return(model)
  
}
```
