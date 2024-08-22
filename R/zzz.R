#When loading package, load Julia packages
#Create Julia environment with user-facing functions from multistate models Julia package
#This environment is accessible in the package environment, so it can be called in the package functions

#If loading Julia packages does not work, install Julia packages first
.onLoad <- function(libname, pkgname) {
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
    multistatemodels_env <- JuliaConnectoR::juliaImport("MultistateModels", all = F)
    assign("multistatemodels_env", multistatemodels_env, envir = topenv())
  },
  error=function(e) {
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

      JuliaConnectoR::juliaEval('
          using MultistateModels
          using CSV
          using DataFrames
          using Distributions
          using LinearAlgebra
          using JLD2
          using ArraysOfArrays
          using Random')

      multistatemodels_env <- JuliaConnectoR::juliaImport("MultistateModels", all = F)
      assign("multistatemodels_env", multistatemodels_env, envir = topenv())
    }
    else {
      stop("Julia setup incorrect.
         Ensure Julia version >= 1.10 is properly installed.")
    }
  }
  )

}

#Stop Julia on unload
.onUnLoad <- function(libpath){
  JuliaConnectoR::stopJulia()
}
