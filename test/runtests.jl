using Pkg

Pkg.activate(".")
Pkg.instantiate()

include("./electricFieldsSimulation/test_risSimulation.jl")