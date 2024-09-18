library(multistatemodels)

df <- read.csv("H:/Projects/Jon/R package and Julia/Julia/derived_states_regen.csv")
d <- collapse_data(df)

#h12 <- Hazard(formula = 0 ~ 1+mab, family = "exp", statefrom = 1, stateto=2)
h12 <- Hazard(formula = 0 ~ 1+mab, family = "sp", statefrom = 1, stateto=2, degree=1, knots=5/7)
h23 <- Hazard(formula = 0 ~ 1+mab, family = "exp", statefrom = 2, stateto=3)
h24 <- Hazard(formula = 0 ~ 1+mab, family = "exp", statefrom = 2, stateto=4)
h45 <- Hazard(formula = 0 ~ 1+mab, family = "exp", statefrom = 4, stateto=5)

h12_list <- JuliaConnectoR::juliaGet(h12)
h23_list <- JuliaConnectoR::juliaGet(h23)

CensoringPatterns <- matrix(c(3, 1, 0, 1, 0, 0), nrow=1, ncol=6)

model <- multistatemodel(hazard = c(h12, h23, h24, h45), data = d[[1]], SamplingWeights = d[[2]],
                         CensoringPatterns = CensoringPatterns)

mod <- initialize_parameters(model)

model_fitted <- fit(mod, tol=1e-3, ess_target_initial=500)
