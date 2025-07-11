using Pkg
ENV["JULIA_APP_BUILD"] = "true"
Pkg.activate(".")
Pkg.instantiate()
ENV["JULIA_APP_BUILD"] = "false"

using Solver
Solver.julia_main()
