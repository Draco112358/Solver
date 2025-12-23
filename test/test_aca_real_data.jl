using Test
using LinearAlgebra
using JSON3, SparseArrays, StaticArrays
include("../src/lib/risACA/ACA_P.jl")
include("../src/lib/risACA/ACA_Lp.jl")
include("../src/lib/risACA/calcola_P_ACA_cpp.jl")
include("../src/lib/risACA/prepare_multilevel.jl")
include("../src/lib/risACA/last_level_fill.jl")
include("../src/lib/risACA/process_surface_ACA.jl")
include("../src/lib/risACA/calcola_P_Sym.jl")
include("../src/lib/risACA/calcola_P_noSym.jl")
include("../src/lib/sharedRis/Song_P_improved_Ivana_strategy.jl")
include("../src/lib/risACA/calcola_Lp_ACA_cpp.jl")
include("../src/lib/risACA/calcola_Lp_Sym.jl")
include("../src/lib/risACA/calcola_Lp_noSym.jl")
include("../src/lib/sharedRis/Song_improved_Ivana_strategy.jl")
include("../src/lib/sharedRis/mean_length_rev.jl")
include("../src/lib/sharedRis/Compute_Lp_Self.jl")
include("../src/lib/risACA/gmres_ACA.jl")
include("../src/lib/risACA/ACA_Lp_delays.jl")
include("../src/lib/risACA/ACA_P_delays.jl")
include("../src/lib/risACA/prepare_sparse_mats_prec_ACA.jl")
include("../src/lib/solverRis/build_Yle_S.jl")
include("../src/lib/risACA/iter_solver_S_type_ACA_ports_sym.jl")

include("../src/lib/sharedRis/calcola_Lp.jl")
include("../src/lib/sharedRis/calcola_P.jl")

function compress_matrix_times_vector(mat_compress, comp, V)
    num_blocks = length(comp)
    res = zeros(ComplexF64, length(V))

    for c1 in 1:num_blocks
        idx = comp[c1]
        if !isempty(idx)
            res[idx] .+= mat_compress[:D][c1] * V[idx]
        end
    end

    for c1 in 1:num_blocks-1
        for c2 in c1+1:num_blocks
            idx1, idx2 = comp[c1], comp[c2]
            if !isempty(idx1) && !isempty(idx2)
                # Use isassigned to avoid UndefRefError with #undef entries
                if isassigned(mat_compress[:U], c1, c2) && isassigned(mat_compress[:V], c1, c2)
                    U_blk = mat_compress[:U][c1, c2]
                    V_blk = mat_compress[:V][c1, c2]

                    if !isnothing(U_blk) && !isnothing(V_blk)
                        res[idx1] .+= U_blk * (V_blk * V[idx2])
                        # Symmetric part (transpose)
                        res[idx2] .+= transpose((transpose(V[idx1]) * U_blk) * V_blk)
                    end
                end
            end
        end
    end
    return res
end


function is_stop_requested(sim_id::String)
    return false
end



# Include necessary files if not already included by Solver
# assuming Solver exports necessary functions or we access them relatively
# For this test script, we'll try to use the functions directly as included in Solver

using Serialization
include("../src/lib/risACA/calcola_Lp_ACA_cpp.jl")

println("--- Starting ACA Test with Real Debug Data ---")

# Load data
debug_file = joinpath(@__DIR__, "aca_debug_data_v2.jls")
if !isfile(debug_file)
    error("Debug file not found: $debug_file")
end

println("Loading data from $debug_file...")
data = deserialize(debug_file)

incidence_selection = data["incidence_selection"]
volumi = data["volumi"]
superfici = data["superfici"]
escalings = data["escalings"]
freq = data["freq"]
ports = data["ports"]
lumped_elements = data["lumped_elements"]
GMRES_settings = data["GMRES_settings"]
use_Zs_in = data["use_Zs_in"]
QS_Rcc_FW = data["QS_Rcc_FW"]
ports_to_excite = [1, 2]
mat_imp_simm = []
id = "id"
chan = nothing
commentsEnabled = true
ACA_thres = 1e-4

println("ACA_thres: ", ACA_thres)


println("Data loaded.")
println("Superfici count: ", size(superfici["estremi_celle"], 1))
println("Volumi count: ", size(volumi[:coordinate], 1))

superfici["estremi_celle"] = hcat(superfici["estremi_celle"]...)
superfici["centri"] = hcat(superfici["centri"]...)


# for (key, value) in superfici
#     println("Key: $key, Value: $(size(value))")
# end
# println("Calling calcola_P...")
P_data = calcola_P_ACA_cpp(superfici, ACA_thres, escalings, "test_params_P")
#P_data_noACA = calcola_P(superfici, escalings, 1, "test_params_P")

# indR = [3, 4, 5]
# indC = [6, 7]
# P_noSym = calcola_P_noSym(indR, indC, superfici, escalings, "id")
# println(P_noSym)
# println(P_data_noACA[:P][indR, indC])
# P_Sym = calcola_P_Sym(indC, superfici, escalings, "id")
# println(P_Sym)
# println(P_data_noACA[:P][indC, indC])
# for (key, value) in superfici
#     println("Key: $key, Value: $(size(value))")
# end


# println("P_data size: ", Base.summarysize(P_data_aca) / (1024^2), " MB")
# println("Calling calcola_Lp...")
Lp_data = calcola_Lp_ACA_cpp(volumi, incidence_selection, ACA_thres, escalings, id)
#Lp_noACA = calcola_Lp(volumi, incidence_selection, escalings, 1, "test_params_Lp")
out = iter_solver_S_type_ACA_ports_sym(freq, escalings, incidence_selection, P_data, Lp_data, ports, lumped_elements, GMRES_settings, volumi, superfici, use_Zs_in, QS_Rcc_FW, ACA_thres, ports_to_excite, mat_imp_simm, id, chan, commentsEnabled)
println(out[:Z])
# indR = [3, 4, 5]
# indC = [6, 7]
# Lp_noSym = calcola_Lp_noSym(indR, indC, volumi, incidence_selection, escalings, "id", "x")
# println(Lp_noSym)
# println(Lp_noACA[:Lp_x][indR, indC])
# Lp_Sym = calcola_Lp_Sym(indC, volumi, incidence_selection, escalings, "id", "x")
# println(Lp_Sym)
# println(Lp_noACA[:Lp_x][indC, indC])

# Lp_noSym = calcola_Lp_noSym(indR, indC, volumi, incidence_selection, escalings, "id", "y")
# println(Lp_noSym)
# println(Lp_noACA[:Lp_y][indR, indC])
# Lp_Sym = calcola_Lp_Sym(indC, volumi, incidence_selection, escalings, "id", "y")
# println(Lp_Sym)
# println(Lp_noACA[:Lp_y][indC, indC])

# Lp_noSym = calcola_Lp_noSym(indR, indC, volumi, incidence_selection, escalings, "id", "z")
# println(Lp_noSym)
# println(Lp_noACA[:Lp_z][indR, indC])
# Lp_Sym = calcola_Lp_Sym(indC, volumi, incidence_selection, escalings, "id", "z")
# println(Lp_Sym)
# println(Lp_noACA[:Lp_z][indC, indC])
# mx = Lp_data[:mx]
# my = Lp_data[:my]
# mz = Lp_data[:mz]
# ns = P_data_aca[:ns]
# m = mx + my + mz
# random = rand(ComplexF64, m)

# randomP = rand(ComplexF64, ns)
# r1P = compress_matrix_times_vector(P_data_aca, P_data_aca[:P_comp], randomP)
# r2P = P_data_noACA[:P] * randomP
# error_norm_P = norm(r1P - r2P) / norm(r2P)
# println("Relative error (P): ", error_norm_P)
# if error_norm_P < 1e-3
#     println("Test passed (ACA product matches dense product)!")
# else
#     println("Test failed (significant difference between ACA and dense)!")
# end


# randomLpx = random[1:mx]
# r1Lpx = compress_matrix_times_vector(Lp_data[:Lpx], Lp_data[:Lpx_comp], randomLpx)
# r2Lpx = Lp_noACA[:Lp_x] * randomLpx
# error_norm_Lpx = norm(r1Lpx - r2Lpx) / norm(r2Lpx)
# println("Relative error (X component): ", error_norm_Lpx)
# if error_norm_Lpx < 1e-3
#     println("Test passed (ACA product matches dense product)!")
# else
#     println("Test failed (significant difference between ACA and dense)!")
# end

# randomLpy = random[mx+1:mx+my]
# r1Lpy = compress_matrix_times_vector(Lp_data[:Lpy], Lp_data[:Lpy_comp], randomLpy)
# r2Lpy = Lp_noACA[:Lp_y] * randomLpy
# error_norm_Lpy = norm(r1Lpy - r2Lpy) / norm(r2Lpy)
# println("Relative error (Y component): ", error_norm_Lpy)
# if error_norm_Lpy < 1e-3
#     println("Test passed (ACA product matches dense product)!")
# else
#     println("Test failed (significant difference between ACA and dense)!")
# end

# randomLpz = random[mx+my+1:end]
# r1Lpz = compress_matrix_times_vector(Lp_data[:Lpz], Lp_data[:Lpz_comp], randomLpz)
# r2Lpz = Lp_noACA[:Lp_z] * randomLpz
# error_norm_Lpz = norm(r1Lpz - r2Lpz) / norm(r2Lpz)
# println("Relative error (Z component): ", error_norm_Lpz)
# if error_norm_Lpz < 1e-3
#     println("Test passed (ACA product matches dense product)!")
# else
#     println("Test failed (significant difference between ACA and dense)!")
# end




