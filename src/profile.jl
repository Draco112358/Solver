using PProf, Profile, JSON, ProfileView, ProfileCanvas
# include("lib/saveFiles.jl")
include("lib/solve.jl")

const stopComputation = []
const commentsEnabled = []

# per usare il profiler, entrare nella shell julia con julia --project=. e lanciare le istruzioni presenti nel file

function force_compile_solver()
    println("------ Precompiling routes...wait for solver to be ready ---------")
    data = open(JSON.parse, "crossbar.json")
    doSolving(data["mesherOutput"], data["solverInput"], data["solverAlgoParams"], data["solverType"], "init"; commentsEnabled=false)
    println("SOLVER READY")
end

#@profview_allocs force_compile_solver()

# Profile.Allocs.clear()

# Profile.Allocs.@profile sample_rate=0.0001 begin
#     force_compile_solver()
# end

# prof = Profile.Allocs.fetch();
# PProf.Allocs.pprof(prof; from_c=false)

Profile.@profile force_compile_solver()
prof = Profile.fetch();
PProf.pprof(prof; from_c=false)





