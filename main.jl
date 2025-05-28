using Pkg

## Questo if sarebbe per fare un check se il modulo è già attivo, in modo da chiudere l'istanza precedente prima di attivare la nuova.
## Ma è da rivedere perché non è chiaro dove il modulo attivato viene agganciato alla variabile instance.
# if @isdefined(instance) && instance !== nothing
#     Solver.stop(instance)
#     global instance = nothing
# end

Pkg.activate(".")
Pkg.instantiate()

#include("src/solver_start.jl")
include("src/solver_start3.jl")
#include("src/restore_message.jl")