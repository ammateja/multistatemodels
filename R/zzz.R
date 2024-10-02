#When loading package, load Julia packages
#Create Julia environment with user-facing functions from multistate models Julia package
#This environment is accessible in the package environment, so it can be called in the package functions

#If loading Julia packages does not work, install Julia packages first
.onLoad <- function(libname, pkgname) {
    JuliaConnectoR::juliaEval('
          using MultistateModels
          using CSV
          using DataFrames')
    multistatemodels_env <- JuliaConnectoR::juliaImport("MultistateModels", all = F)
    assign("multistatemodels_env", multistatemodels_env, envir = topenv())
}

#Stop Julia on unload
.onUnload <- function(libpath){
  #dev.off()
  JuliaConnectoR::stopJulia()
}
