using JLD2, Test, MAT, LinearAlgebra, SparseArrays, StaticArrays, JSON, JSON3, AWS, AWSS3, Base.Threads, FFTW, Printf
include("../../src/lib/solverRisCampi/do_solving_electric_fields.jl")
include("../../src/lib/utility.jl") # Contiene download_json_gz, get_solverInput_from_s3, ecc.
include("../../src/lib/format_input_output_solver_functions.jl")

#SHAREDRIS INCLUDE START

include("../../src/lib/sharedRis/Compute_Lp_Self.jl")
include("../../src/lib/sharedRis/Song_P_improved_Ivana_strategy.jl")
include("../../src/lib/sharedRis/Song_improved_Ivana_strategy.jl")
include("../../src/lib/sharedRis/distfcm.jl")
include("../../src/lib/sharedRis/find_nodes_ports_or_le.jl")
include("../../src/lib/sharedRis/calcola_P.jl")
include("../../src/lib/sharedRis/calcola_Lp.jl")
#SHAREDRIS INCLUDE END

# SOLVER RSI CAMPI INCLUDE START
include("../../src/lib/solverRisCampi/get_signal.jl")
include("../../src/lib/solverRisCampi/build_trapezoidal_pulse.jl")
include("../../src/lib/solverRisCampi/genera_segnale_esponenziale.jl")
include("../../src/lib/solverRisCampi/genera_segnale_Gaussiano_modulato.jl")
include("../../src/lib/solverRisCampi/genera_segnale_sinusoidale.jl")
include("../../src/lib/solverRisCampi/crea_freqs.jl")
include("../../src/lib/solverRisCampi/compute_fields_components.jl")
include("../../src/lib/solverRisCampi/computeVs.jl")
include("../../src/lib/solverRisCampi/fft_UAq.jl")
include("../../src/lib/solverRisCampi/complex_matrix_to_float_array_matrix.jl")
include("../../src/lib/solverRisCampi/genera_punti_circonferenza.jl")
include("../../src/lib/solverRisCampi/get_punti_oss_3D.jl")
include("../../src/lib/solverRisCampi/build_Yle.jl")
include("../../src/lib/solverRisCampi/compute_Matrix_vector.jl")
include("../../src/lib/solverRisCampi/gmres_custom.jl")
include("../../src/lib/solverRisCampi/compute_Ec_Gauss.jl")
include("../../src/lib/solverRisCampi/compute_Ar_Gauss.jl")
include("../../src/lib/solverRisCampi/compute_lambda_numeric.jl")
include("../../src/lib/solverRisCampi/compute_E_field_Gauss.jl")
include("../../src/lib/solverRisCampi/compute_H_field_Gauss.jl")
include("../../src/lib/solverRisCampi/iter_solver_E_Gaussian_Is_type.jl")
include("../../src/lib/solverRisCampi/do_solving_electric_fields.jl")
# SOLVER RSI CAMPI INCLUDE END

function is_stop_requested(sim_id::String)
    lock(stop_computation_lock) do
        return haskey(stopComputation, sim_id) && stopComputation[sim_id][]
    end
end

const solver_overall_status = Ref("ready") # ready, busy, error
const active_simulations = Dict{String, Dict{String, Any}}() # ID simulazione -> {status, progress, start_time, etc.}
const simulations_lock = ReentrantLock() # Lock per proteggere `active_simulations`
const stopComputation = Dict{String, Ref{Bool}}() # ID simulazione -> Ref{Bool} per il flag di stop
const stop_computation_lock = ReentrantLock() # Aggiungi un lock per proteggere stopComputation
const commentsEnabled = []
const norm_treshold = 1e-10
test_suites = [
        ("Test Tx ris 8x8", () -> @testset "Test Tx ris 8x8" begin
        @load "./test/electricFieldsSimulation/test_input_Tx8x8.jld2"
        out = doSolvingElectricFields(incidence_selection, volumi, superfici, nodi_coord, escalings, solverInput, solverAlgoParams, solverType, theta, phi, e_theta, e_phi, baricentro, r_circ, times, signal_type_E, ind_freq_interest, id, aws_config, bucket_name; chan=nothing, commentsEnabled=false)
        outMAT = matread("./test/electricFieldsSimulation/Tx_8x8.mat")
        normExJulia = norm(out["Ex"])
        normExMAT = norm(outMAT["out"]["Ex"])
        normEyJulia = norm(out["Ey"])
        normEyMAT = norm(outMAT["out"]["Ey"])
        normEzJulia = norm(out["Ez"])
        normEzMAT = norm(outMAT["out"]["Ez"])
        @test (normExMAT-normExJulia)/normExMAT < norm_treshold
        @test (normEyMAT-normEyJulia)/normEyMAT < norm_treshold 
        @test (normEzMAT-normEzJulia)/normEzMAT < norm_treshold
        end),
        
        ("Test ris 6x6", () -> @testset "Test ris 6x6" begin
            @load "./test/electricFieldsSimulation/test_input_ris6x6.jld2"
            out = doSolvingElectricFields(incidence_selection, volumi, superfici, nodi_coord, escalings, solverInput, solverAlgoParams, solverType, theta, phi, e_theta, e_phi, baricentro, r_circ, times, signal_type_E, ind_freq_interest, id, aws_config, bucket_name; chan=nothing, commentsEnabled=false)
            outMAT = matread("./test/electricFieldsSimulation/ris_6x6.mat")
            normExJulia = norm(out["Ex"])
            normExMAT = norm(outMAT["out"]["Ex"])
            normEyJulia = norm(out["Ey"])
            normEyMAT = norm(outMAT["out"]["Ey"])
            normEzJulia = norm(out["Ez"])
            normEzMAT = norm(outMAT["out"]["Ez"])
            @test (normExMAT-normExJulia)/normExMAT < norm_treshold 
            @test (normEyMAT-normEyJulia)/normEyMAT < norm_treshold 
            @test (normEzMAT-normEzJulia)/normEzMAT < norm_treshold
        end),
        
        ("Test ris 8x8", () -> @testset "Test ris 8x8" begin
            @load "./test/electricFieldsSimulation/test_input_ris8x8.jld2"
            out = doSolvingElectricFields(incidence_selection, volumi, superfici, nodi_coord, escalings, solverInput, solverAlgoParams, solverType, theta, phi, e_theta, e_phi, baricentro, r_circ, times, signal_type_E, ind_freq_interest, id, aws_config, bucket_name; chan=nothing, commentsEnabled=false)
            outMAT = matread("./test/electricFieldsSimulation/ris_8x8.mat")
            normExJulia = norm(out["Ex"])
            normExMAT = norm(outMAT["out"]["Ex"])
            normEyJulia = norm(out["Ey"])
            normEyMAT = norm(outMAT["out"]["Ey"])
            normEzJulia = norm(out["Ez"])
            normEzMAT = norm(outMAT["out"]["Ez"])
            @test (normExMAT-normExJulia)/normExMAT < norm_treshold 
            @test (normEyMAT-normEyJulia)/normEyMAT < norm_treshold 
            @test (normEzMAT-normEzJulia)/normEzMAT < norm_treshold
        end)
    ]

function send_rabbitmq_feedback(data::Dict, routing_key::String)
    println("Feedback non inviato perchè siamo in fase di test")
end

function risSimulationTest(test_suites)
    for (i, (suite_name, test_func)) in enumerate(test_suites)
        println("Eseguendo $suite_name ($(i)/$(length(test_suites)))...")
        
        results = test_func()
        
        if results.anynonpass
            println("❌ $suite_name fallito.")
            println("   Passati: $(results.n_passed)")
            println("   Falliti: $(results.n_failed)")
            println("   Errori: $(results.n_error)")
            println("⏹️  Interrompo l'esecuzione delle suite successive.")
            return false
        end
        
        println("✅ $suite_name completato ($(results.n_passed) test passati)")
    end
    
    println("🎉 Tutte le $(length(test_suites)) test suite di risSimulation completate con successo!")
    return true
end

# Esegui i test
success = risSimulationTest(test_suites)