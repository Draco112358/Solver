using Pkg

ENV["JULIA_APP_BUILD"] = "true"
Pkg.activate(".")
Pkg.instantiate()
ENV["JULIA_APP_BUILD"] = "false"

include("./electricFieldsSimulation/test_risSimulation.jl")