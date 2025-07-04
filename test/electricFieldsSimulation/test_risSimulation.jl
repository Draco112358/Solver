using JLD2, Test, MAT, LinearAlgebra
include("../../src/lib/solve3.jl")

function is_stop_requested(sim_id::String)
    lock(stop_computation_lock) do
        return haskey(stopComputation, sim_id) && stopComputation[sim_id][]
    end
end

const solver_overall_status = Ref("ready") # ready, busy, error
const active_simulations = Dict{String, Dict{String, Any}}() # ID simulazione -> {status, progress, start_time, etc.}
const simulations_lock = ReentrantLock() # Lock per proteggere `active_simulations`
# const stopComputation = []
const stopComputation = Dict{String, Ref{Bool}}() # ID simulazione -> Ref{Bool} per il flag di stop
const stop_computation_lock = ReentrantLock() # Aggiungi un lock per proteggere stopComputation
const commentsEnabled = []

function send_rabbitmq_feedback(data::Dict, routing_key::String)
    println("Feedback non inviato perchè siamo in fase di test")
end


@load "./test/electricFieldsSimulation/test_input_electricFieldsSimulation.jld2"
out = doSolvingElectricFields(incidence_selection, volumi, superfici, nodi_coord, escalings, solverInput, solverAlgoParams, solverType, theta, phi, e_theta, e_phi, baricentro, r_circ, times, signal_type_E, ind_freq_interest, id, aws_config, bucket_name)
outMAT = matread("./test/electricFieldsSimulation/ris_8x8_new.mat")
normExJulia = norm(out["Ex"])
normExMAT = norm(outMAT["out"]["Ex"])
normEyJulia = norm(out["Ey"])
normEyMAT = norm(outMAT["out"]["Ey"])
normEzJulia = norm(out["Ez"])
normEzMAT = norm(outMAT["out"]["Ez"])
@testset "Test Electric Fields Simulation" begin
    @test (normExMAT-normExJulia)/normExMAT < 1e-10
    @test (normEyMAT-normEyJulia)/normEyMAT < 1e-10
    @test (normEzMAT-normEzJulia)/normEzMAT < 1e-10
end