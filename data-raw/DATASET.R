library(multistatemodels)
library(dplyr)

set.seed(1)
nsubj <- 250

make_obstimes <- function(ntimes=12) {
  #observation times
  interval <- seq(0, 1, length.out = ntimes+1)
  times <- interval[-c(1, length(interval))] + (rbeta(length(interval[-c(1, length(interval))]), 1.5, 1.5) - 0.5)*diff(interval)[1]
  times <- c(0, times, 1)

  return(times)
}

h12 <- multistatemodels::Hazard(formula = 0~1, statefrom = 1, stateto=2, family="wei")
h13 <- multistatemodels::Hazard(formula = 0~1, statefrom = 1, stateto=3, family="wei")
h23 <- multistatemodels::Hazard(formula = 0~1, statefrom = 2, stateto=3, family="wei")

dat <- NULL
for (i in 1:nsubj) {
  visitdays <- make_obstimes()
  id <- rep(i, ntimes)
  tstart <- visitdays[-length(visitdays)]
  tstop <- visitdays[-1]
  statefrom <- rep(1, ntimes)
  stateto <- rep(1, ntimes)
  obstype <- rep(2, ntimes)
  d <- data.frame(id=id, tstart=tstart, tstop=tstop, statefrom=statefrom, stateto=stateto, obstype=obstype)
  dat <- rbind(dat, d)
}

model <- multistatemodels::multistatemodel(hazard = c(h12, h13, h23), data=dat)
parameters <- makepars()
model_sim <- multistatemodels::set_parameters(model=model, newvalues=list(h12 = c(log(1.25), log(1.5)),
                                                                          h13 = c(log(1.25), log(1)),
                                                                          h23 = c(log(1.25), log(2))))

paths <- multistatemodels::simulate(model=model_sim, nsim=1, paths=TRUE, data=FALSE)


d <- as.data.frame(model_sim$data)
illness_death_dat <- NULL

for (i in 1:length(paths)) {

  subj.dat.raw <- d[d$id == i, ]
  obstimes <- unique(sort(c(0, subj.dat.raw$tstop, paths[[i]]$times[paths[[i]]$states ==3])))

  if (sum(c(2,3) %in% paths[[i]]$states) == 2) {

    new <- paths[[i]]$times[length(paths[[i]]$times)] - sqrt(.Machine$double.eps)
    obstimes <- sort(c(obstimes, new))

  }

  obsinds <- unlist(lapply(obstimes, function(x){max(which(paths[[i]]$times <= x))}))
  obsstates <- paths[[i]]$states[obsinds]

  df <- data.frame(id = rep(paths[[i]]$subj, (length(obstimes)-1)),
                   tstart = obstimes[-length(obstimes)],
                   tstop = obstimes[-1],
                   statefrom = obsstates[-length(obsstates)],
                   stateto = obsstates[-1])
  df <- df[df$stateto != 3 | df$statefrom != 3, ]

  df$obstype <- ifelse(df$stateto == 3, 1, 2)
  df$obstype <- ifelse(df$statefrom == df$stateto, 1, df$obstype)

  illness_death_dat <- rbind(illness_death_dat, df)

}

illness_death_dat <- illness_death_dat %>% dplyr::arrange(id, tstart)

usethis::use_data(illness_death_dat, overwrite = TRUE)
