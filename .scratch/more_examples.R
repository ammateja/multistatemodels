summarize_crude <- function(paths, dat, model) {
  
  nsubj <- length(paths)
  events <- as.numeric(get_estimates(paths))
  cis <- binom::binom.confint(events[1:4]*nsubj, nsubj, method = "wilson")[ ,c("lower", "upper")]
  
  times <- unlist(lapply(paths, function(x) x$times[2]))
  statuses <- unlist(lapply(paths, function(x) ifelse(x$states[2] == 1, 0, 1)))
  
  sm1 <- survival:::survmean(survfit(Surv(times, statuses) ~ 1), rmean=max(times))[[1]][c("rmean", "se(rmean)")]
  rmst <- data.frame(est = sm1[1], lower = sm1[1] - 1.96*sm1[2], upper = sm1[1] + 1.96*sm1[2])
  
  times2prog <- unlist(lapply(paths, function(x) ifelse(x$states[2] == 2, x$times[2], 1)))
  prog_status <- unlist(lapply(paths, function(x) ifelse(x$states[2] == 2, 1, 0)))
  
  proginds <- unlist(lapply(paths, function(x){2 %in% x$states & length(x$states) != 2}))
  progtimes <- unlist(lapply(paths[proginds], function(x){x$times[2]}))
  illnessdurs <- unlist(lapply(paths[proginds], function(x){x$times[3] - x$times[2]}))
  
  sm2 <- survival:::survmean(survfit(Surv(progtimes, rep(1, length(progtimes))) ~ 1), rmean=max(times))[[1]][c("rmean", "se(rmean)")]
  time2prog_gprog <- data.frame(est = sm2[1], lower = sm2[1] - 1.96*sm2[2], upper = sm2[1] + 1.96*sm2[2])
  
  sm3 <- survival:::survmean(survfit(Surv(times2prog, prog_status) ~ 1), rmean=max(times))[[1]][c("rmean", "se(rmean)")]
  time2prog_all <- data.frame(est = sm3[1], lower = sm3[1] - 1.96*sm3[2], upper = sm3[1] + 1.96*sm3[2])
  
  sm4 <- survival:::survmean(survfit(Surv(illnessdurs, rep(1, length(illnessdurs))) ~ 1), rmean=max(times))[[1]][c("rmean", "se(rmean)")]
  illnessdur <- data.frame(est = sm4[1], lower = sm4[1] - 1.96*sm4[2], upper = sm4[1] + 1.96*sm4[2])
  
  ests_crude <- data.frame(ests = c(events[1:4], rmst[ ,1], time2prog_gprog[ ,1], time2prog_all[ ,1], illnessdur[ ,1]), 
                           lower = c(cis$lower, rmst[ ,2], time2prog_gprog[ ,2], time2prog_all[ ,2], illnessdur[ ,2]), 
                           upper = c(cis$upper, rmst[ ,3], time2prog_gprog[ ,3], time2prog_all[ ,3], illnessdur[ ,3]))
  
  return(ests_crude)
  
}


crude_ests <- function(seed, ntimes) {
  
  set.seed(seed)
  
  model_sim <- setup_model(make_pars=TRUE, data = NULL, family = "wei", ntimes = ntimes, nsubj = 250)
  
  paths <- multistatemodels::simulate(model_sim, nsim=1, paths=TRUE, data=FALSE)
  dat <- observe_subjdat(paths, model_sim)
  
  paths <- JuliaConnectoR::juliaGet(JuliaConnectoR::juliaLet("MultistateModels.extract_paths(dat)", dat=JuliaConnectoR::juliaCall("DataFrame", dat)))
  
  ests_crude <- summarize_crude(paths, dat, model_sim)
  
  results <- data.frame(simnum = rep(seed, 8), 
                        ntimes = rep(ntimes, 8), 
                        family = rep("crude", 8), 
                        var = c("pfs", "prog", "die_wprog", "die_noprog", "rmst", "time2prog_gprog", "time2prog_all", "illnessdur"), 
                        ests = ests_crude[ ,1], 
                        lower = ests_crude[, 2], 
                        upper = ests_crude[ ,3])
  
  return(results)
  
}