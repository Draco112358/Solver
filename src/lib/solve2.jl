include("FFT_solver_QS_S_type.jl")
include("create_volume_centers.jl")
include("create_Grids_externals.jl")
include("compute_FFT_mutual_coupling_mats.jl")
include("mesher_FFT.jl")
include("From_3D_to_1D.jl")
include("utility.jl")
include("iter_solver_QS_S_type2.jl")
include("find_nodes_ports_or_le.jl")
include("calcola_P.jl")
include("calcola_Lp2.jl")
include("genera_segnale_esponenziale.jl")
include("genera_segnale_Gaussiano_modulato.jl")
include("genera_segnale_sinusoidale.jl")
include("crea_freqs.jl")
include("build_trapezoidal_pulse.jl")
include("compute_fields_components.jl")
include("computeVs.jl")
include("fft_UAq.jl")
include("genera_punti_circonferenza.jl")
include("get_punti_oss_3D.jl")
include("iter_solver_E_Gaussian_Is_type2.jl")

using MKL
using JSON, AWSS3, BSON
using MLUtils: unsqueeze
function dump_json_data(matrix_Z, matrix_S, matrix_Y, num_ports, id; partial=false, freqIndex=nothing)
    z = [[[[0.1, 0.0]]]]
    pop!(z)
    s = similar(z)
    y = similar(z)
    matrix_Z = convert(Array{ComplexF64,3}, matrix_Z)
    matrix_S = convert(Array{ComplexF64,3}, matrix_S)
    matrix_Y = convert(Array{ComplexF64,3}, matrix_Y)
    for i in range(1, num_ports)
        for j in range(1, num_ports)
            elements_Z = map(v -> reinterpret(Float64, [v]), matrix_Z[i, j, :])
            elements_S = map(v -> reinterpret(Float64, [v]), matrix_S[i, j, :])
            elements_Y = map(v -> reinterpret(Float64, [v]), matrix_Y[i, j, :])
            push!(z, [elements_Z])
            push!(s, [elements_S])
            push!(y, [elements_Y])
        end
    end
    solver_matrices_dict = Dict(
        "matrices" => Dict(
            "matrix_Z" => JSON.json(z),
            "matrix_S" => JSON.json(s),
            "matrix_Y" => JSON.json(y),
        ),
        "isStopped" => false,
        "id" => id,
        "partial" => partial,
        "freqIndex" => freqIndex
    )
    return solver_matrices_dict
end


function material(dict_element)
    return Dict(
        "name" => dict_element["name"],
        "color" => dict_element["color"],
        "permeability" => dict_element["permeability"],
        "tangent_delta_permeability" => haskey(dict_element, "tangent_delta_permeability") ? dict_element["tangent_delta_permeability"] : 0,
        "custom_permeability" => haskey(dict_element, "custom_permeability") ? dict_element["custom_permeability"] : 0,
        "permittivity" => dict_element["permittivity"],
        "tangent_delta_permittivity" => haskey(dict_element, "tangent_delta_permittivity") ? dict_element["tangent_delta_permittivity"] : 0,
        "custom_permittivity" => haskey(dict_element, "custom_permittivity") ? dict_element["custom_permittivity"] : 0,
        "conductivity" => dict_element["conductivity"],
        "tangent_delta_conductivity" => haskey(dict_element, "tangent_delta_conductivity") ? dict_element["tangent_delta_conductivity"] : 0,
        "custom_conductivity" => haskey(dict_element, "custom_conductivity") ? dict_element["custom_conductivity"] : 0,
        "sigmar" => dict_element["conductivity"],
        "tan_D" => haskey(dict_element, "tangent_delta_conductivity") ? dict_element["tangent_delta_conductivity"] : 0,
        "eps_re" => dict_element["permittivity"],
        "mur" => dict_element["permeability"],
        "epsr" => 1,
        "Rx" => nothing,
        "Ry" => nothing,
        "Rz" => nothing,
        "Cx" => nothing,
        "Cy" => nothing,
        "Cz" => nothing
    )
end


function port_def(port_start, port_end, port_voxels, port_nodes, surf_s_port_nodes, surf_e_port_nodes)
    return Dict(
        "port_start" => port_start,
        "port_end" => port_end,
        "port_voxels" => port_voxels,
        "port_nodes" => port_nodes,
        "surf_s_port_nodes" => surf_s_port_nodes,
        "surf_e_port_nodes" => surf_e_port_nodes
    )
end

function le_def(value, type, le_start, le_end, le_voxels, le_nodes, surf_s_le_nodes, surf_e_le_nodes, R_value, L_value, C_value)
    return Dict(
        "value" => value,
        "type" => type,
        "R_value" => R_value,
        "L_value" => L_value,
        "C_value" => C_value,
        "le_start" => le_start,
        "le_end" => le_end,
        "le_voxels" => le_voxels,
        "le_nodes" => le_nodes,
        "surf_s_le_nodes" => surf_s_le_nodes,
        "surf_e_le_nodes" => surf_e_le_nodes
    )
end


function read_ports(port_objects, escal)
    #@assert inputData isa Dict
    input_positions = []
    output_positions = []
    N_PORTS = length(port_objects)

    for port_object in port_objects
        @assert length(port_object["inputElement"]) == 3
        ipos = zeros((1, 3))
        ipos[1, 1] = port_object["inputElement"][1] * escal
        ipos[1, 2] = port_object["inputElement"][2] * escal
        ipos[1, 3] = port_object["inputElement"][3] * escal
        push!(input_positions, ipos)
        @assert length(port_object["outputElement"]) == 3
        opos = zeros((1, 3))
        opos[1, 1] = port_object["outputElement"][1] * escal
        opos[1, 2] = port_object["outputElement"][2] * escal
        opos[1, 3] = port_object["outputElement"][3] * escal
        push!(output_positions, opos)
    end
    @assert length(input_positions) == N_PORTS && length(output_positions) == N_PORTS
    inp_pos = []
    for i in input_positions
        push!(inp_pos, unsqueeze([i], dims=2))
    end
    out_pos = []
    for i in output_positions
        push!(out_pos, unsqueeze([i], dims=2))
    end
    ports_out = port_def(inp_pos, out_pos, zeros(Int64, (N_PORTS, 2)), zeros(Int64, (N_PORTS, 2)), Array{Any}(undef, 0), Array{Any}(undef, 0))

    return ports_out
end


function read_lumped_elements(lumped_elements_objects, escal)

    #@assert inputData isa Dict
    input_positions = []
    output_positions = []
    values = []
    types = []
    R_values = []
    L_values = []
    C_values = []
    N_LUMPED_ELEMENTS = length(lumped_elements_objects)
    if N_LUMPED_ELEMENTS == 0
        lumped_elements_out = le_def(zeros(0), zeros(Int64, 0), zeros((0, 3)), zeros((0, 3)), zeros(Int64, (0, 2)), zeros(Int64, (0, 2)), Array{Any}(undef, 0), Array{Any}(undef, 0), R_values, L_values, C_values)
        @assert length(input_positions) == N_LUMPED_ELEMENTS && length(output_positions) == N_LUMPED_ELEMENTS && length(values) == N_LUMPED_ELEMENTS && length(types) == N_LUMPED_ELEMENTS
    else
        for lumped_element_object in lumped_elements_objects
            @assert length(lumped_element_object["inputElement"]) == 3
            ipos = zeros((1, 3))
            ipos[1, 1] = lumped_element_object["inputElement"][1] * escal
            ipos[1, 2] = lumped_element_object["inputElement"][2] * escal
            ipos[1, 3] = lumped_element_object["inputElement"][3] * escal
            push!(input_positions, ipos)
            @assert length(lumped_element_object["outputElement"]) == 3
            opos = zeros((1, 3))
            opos[1, 1] = lumped_element_object["outputElement"][1] * escal
            opos[1, 2] = lumped_element_object["outputElement"][2] * escal
            opos[1, 3] = lumped_element_object["outputElement"][3] * escal
            push!(output_positions, opos)
            lvalue = zeros(1)
            lvalue[1] = lumped_element_object["value"]
            append!(values, lvalue)

            ltype = zeros(Int64, 1)
            ltype[1] = lumped_element_object["type"]
            push!(types, ltype)

            r_value = zeros(Float64, 1)
            r_value[1] = haskey(lumped_element_object["rlcParams"], "resistance") ? lumped_element_object["rlcParams"]["resistance"] : 0.0
            push!(R_values, r_value)

            l_value = zeros(Float64, 1)
            l_value[1] = haskey(lumped_element_object["rlcParams"], "inductance") ? lumped_element_object["rlcParams"]["inductance"] : 0.0
            push!(L_values, l_value)

            c_value = zeros(Float64, 1)
            c_value[1] = haskey(lumped_element_object["rlcParams"], "capacitance") ? lumped_element_object["rlcParams"]["capacitance"] : 0.0
            push!(C_values, c_value)
        end

        @assert length(input_positions) == N_LUMPED_ELEMENTS && length(output_positions) == N_LUMPED_ELEMENTS && length(values) == N_LUMPED_ELEMENTS && length(types) == N_LUMPED_ELEMENTS && length(R_values) == N_LUMPED_ELEMENTS && length(L_values) == N_LUMPED_ELEMENTS && length(C_values) == N_LUMPED_ELEMENTS

        value = [i[1] for i in values]
        type = [i[1] for i in types]
        R_value = [i[1] for i in R_values]
        L_value = [i[1] for i in L_values]
        C_value = [i[1] for i in C_values]
        in_pos = [unsqueeze([i], dims=2) for i in input_positions]
        out_pos = [unsqueeze([i], dims=2) for i in output_positions]

        lumped_elements_out = le_def(value, type, in_pos, out_pos, zeros(Int64, (N_LUMPED_ELEMENTS, 2)), (Int64, (N_LUMPED_ELEMENTS, 2)), Array{Any}(undef, 0), Array{Any}(undef, 0), R_value, L_value, C_value)
    end
    return lumped_elements_out
end


function getEscalFrom(unit)
    escal = 1.0
    if (unit == "m")
        escal = 1.0
    end
    if (unit == "dm")
        escal = 1e-1
    end
    if (unit == "cm")
        escal = 1e-2
    end
    if (unit == "mm")
        escal = 1e-3
    end
    if (unit == "microm")
        escal = 1e-6
    end
    if (unit == "nanom")
        escal = 1e-9
    end
    return escal
end

function create_volumes_mapping_v2(grids)
    num_grids = length(grids)
    #println(size(grids))
    #Nx, Ny, Nz = size(grids[1])
    Nx = size(grids[1], 1)
    Ny = size(grids[1][1], 1)
    Nz = size(grids[1][1][1], 1)
    mapping = zeros(Nx * Ny * Nz)
    num_ele = 0
    for cont3 = 1:Nz
        for cont2 = 1:Ny
            for cont = 1:Nx
                for k = 1:num_grids
                    if grids[k][cont][cont2][cont3] != 0
                        num_ele += 1
                        mapping[from_3D_to_1D(cont, cont2, cont3, Nx, Ny)] = num_ele
                        break
                    end
                end
            end
        end
    end
    return mapping, num_ele
end

function isMaterialConductor(materialName::String, materials)::Bool
    material = materials[findfirst(m -> m["name"] == materialName, materials)]
    return material["conductivity"] > 0.0
end

mutable struct ElectricFieldsSolvingState
    id::String
    incidence_selection::Any
    volumi::Any
    superfici::Any
    nodi_coord::Any
    escalings::Any
    solverInput::Dict
    solverAlgoParams::Dict
    solverType::Int64
    theta::Float64
    phi::Float64
    e_theta::Float64
    e_phi::Float64
    baricentro::Vector{Float64}
    r_circ::Float64
    times::Vector{Float64}
    signal_type_E::String
    ind_freq_interest::Vector{Int}
    currentState::String
    currentStateComputeFields::String
    progress::Float64
    intermediateResults::Dict{String, Any} # Per memorizzare dati tra gli step
end

function state_to_dict(state::ElectricFieldsSolvingState)
    return Dict(
        "id" => state.id,
        "incidence_selection" => state.incidence_selection,
        "volumi" => state.volumi,
        "superfici" => state.superfici,
        "nodi_coord" => state.nodi_coord,
        "escalings" => state.escalings,
        "solverInput" => state.solverInput,
        "solverAlgoParams" => state.solverAlgoParams,
        "solverType" => state.solverType,
        "theta" => state.theta,
        "phi" => state.phi,
        "e_theta" => state.e_theta,
        "e_phi" => state.e_phi,
        "baricentro" => state.baricentro,
        "r_circ" => state.r_circ,
        "times" => state.times,
        "signal_type_E" => state.signal_type_E,
        "ind_freq_interest" => state.ind_freq_interest,
        "currentState" => state.currentState,
        "currentStateComputeFields" => state.currentStateComputeFields,
        "progress" => state.progress,
        "intermediateResults" => state.intermediateResults,
    )
end

function serialize_with_symbols_as_strings_recursive(data::Any)
    if data isa Dict
        serialized = Dict{String, Any}()
        for (key, value) in data
            new_key = key isa Symbol ? "_sym_" * String(key) : String(key)
            serialized[new_key] = serialize_with_symbols_as_strings_recursive(value)
        end
        return serialized
    else
        return data
    end
end

function upload_serialized_data(aws_config, bucket_name, file_name, data_to_save)
    io = IOBuffer()
    # Serialize the variable into the IOBuffer.
    Serialization.serialize(io, data_to_save)
    
    # Get the bytes from the IOBuffer.
    data_bytes = take!(io)
    s3_put(aws_config, bucket_name, file_name, data_bytes)
end

function doSolvingElectricFields2(state::ElectricFieldsSolvingState, aws_config, bucket_name; chan=nothing, commentsEnabled=true, stopComputation)
    task_id = state.id

    try
        if state.currentState == "initial"
            println("[$task_id] Generazione segnale...")
            time_delay_vs = 0.0
            f0 = 0.0
            dev_stand = 0.0
            tr = 0.0
            power = 0.0
            if state.signal_type_E == "exponential"
                time_delay_vs = 3e-9
                tw = 50 * 0.1 / 3e8
                power = 4.0
                vs, tr = genera_segnale_esponenziale(tw, power, state.times, time_delay_vs)
            elseif state.signal_type_E == "gaussian_modulated"
                time_delay_vs = 3e-9
                f0 = 1e9
                dev_stand = 10 / (4 * pi) * 1e-9
                vs = genera_segnale_Gaussiano_modulato(f0, dev_stand, state.times, time_delay_vs)
            elseif state.signal_type_E == "sinusoidal"
                time_delay_vs = 3e-9
                f0 = 1e8
                vs = genera_segnale_sinusoidale(f0, state.times, time_delay_vs)
            end
            state.intermediateResults["vs"] = vs
            state.currentState = "signal_generated"
            save_task_state(state, aws_config, bucket_name, "electric_fields")
        end

        if state.currentState == "signal_generated"
            println("[$task_id] Lettura porte...")
            inputDict = state.solverInput
            unit = state.solverInput["unit"]
            escal = getEscalFrom(unit)
            ports_scatter_value = haskey(inputDict, "ports_scattering_value") ? inputDict["ports_scattering_value"] : 50.0
            state.intermediateResults["escal"] = escal
            state.intermediateResults["ports_scatter_value"] = ports_scatter_value
            ports, lumped_elements = find_nodes_ports_or_le(inputDict["ports"], inputDict["lumped_elements"], state.nodi_coord, escal)
            state.intermediateResults["ports"] = ports
            state.intermediateResults["lumped_elements"] = lumped_elements
            state.currentState = "ports_read"
            save_task_state(state, aws_config, bucket_name, "electric_fields")
        end

        if state.currentState == "ports_read" || state.currentState == "compute fields"
            println("[$task_id] Calcolo P e Lp ed esecuzione GMRES...")
            QS_Rcc_FW = state.solverType
            P_data = state.currentState == "ports_read" ? calcola_P(state.superfici, state.escalings, QS_Rcc_FW) : []
            if !isnothing(chan)
                publish_data(Dict("computingP" => true, "id" => task_id), "solver_feedback", chan)
            end
            if is_stopped_computation(task_id, stopComputation)
                return false
            end
            Lp_data = state.currentState == "ports_read" ? calcola_Lp2(state.volumi, state.incidence_selection, state.escalings, QS_Rcc_FW) : []
            if !isnothing(chan)
                publish_data(Dict("computingLp" => true, "id" => task_id), "solver_feedback", chan)
            end
            if is_stopped_computation(task_id, stopComputation)
                return false
            end
            ind_freq_interest_adj = state.ind_freq_interest .+ 1
            freq_vs_is = crea_freqs(state.times)
            freq = freq_vs_is[ind_freq_interest_adj]
            n_freq = length(freq)
            GMRES_settings = Dict("Inner_Iter" => state.solverAlgoParams["innerIteration"], "Outer_Iter" => state.solverAlgoParams["outerIteration"], "tol" => state.solverAlgoParams["convergenceThreshold"] * ones((n_freq)))
            ind_low_freq = findall(x -> x < 1e5, freq)
            GMRES_settings["tol"][ind_low_freq] .= 1e-8
            use_Zs_in = true

            is_matrix = zeros(size(state.intermediateResults["ports"][:port_start], 1), length(state.times))
            for (index, signal_type) in enumerate(state.intermediateResults["ports"][:signals_port])
                if signal_type != "no_signal"
                    vs = getSignalbasedOn(signal_type, state.times)
                    is_matrix[index, :] = vs
                end
            end

            row_indices, col_indices, nz_values = findnz(state.incidence_selection[:A])
            A = sparse(row_indices, col_indices, nz_values)
            E, K, H, E_theta_v, E_phi_v = compute_fields_components(state.theta, state.phi, state.e_theta, state.e_phi)
            Vs = computeVs(state.times, 0.0, state.signal_type_E, state.volumi, state.nodi_coord, E, K, A, 0.0, 0.0, 0.0, 0.0, ind_freq_interest_adj) # time_delay_vs, tr, power, dev_stand, f0 sono gestiti prima

            Is = zeros(ComplexF64, size(state.intermediateResults["ports"][:port_nodes], 1), n_freq)
            for k in 1:size(state.intermediateResults["ports"][:port_nodes], 1)
                Trasformata = fft_UAq(state.times, is_matrix[k, :])
                Is[k, :] = Trasformata[2, ind_freq_interest_adj]
            end

            r_circ_scaled = state.r_circ * state.intermediateResults["escal"]
            baricentro_scaled = state.baricentro .* state.intermediateResults["escal"]
            punti_xy = genera_punti_circonferenza(r_circ_scaled, 100, baricentro_scaled, 1)
            punti_zx = genera_punti_circonferenza(r_circ_scaled, 100, baricentro_scaled, 2)
            punti_yz = genera_punti_circonferenza(r_circ_scaled, 100, baricentro_scaled, 3)
            centri_oss = [punti_xy; punti_zx; punti_yz]
            centri_oss_3D, distanze_3D, theta_vals, x_grid, y_grid, z_grid = get_punti_oss_3D(r_circ_scaled, 50, baricentro_scaled)
            state.intermediateResults["centri_oss"] = centri_oss
            state.intermediateResults["centri_oss_3D"] = centri_oss_3D
            state.intermediateResults["distanze_3D"] = distanze_3D
            state.intermediateResults["theta_vals"] = theta_vals
            state.intermediateResults["x_grid"] = x_grid
            state.intermediateResults["y_grid"] = y_grid
            state.intermediateResults["z_grid"] = z_grid
            out = iter_solver_E_Gaussian_Is_type2(
                freq, state.escalings, state.incidence_selection, P_data, Lp_data,
                state.intermediateResults["ports"], state.intermediateResults["lumped_elements"], GMRES_settings, state.volumi, state.superfici, use_Zs_in, state.solverType, state.intermediateResults["ports_scatter_value"], Vs, Is, centri_oss, centri_oss_3D, task_id, chan, commentsEnabled,
                state, aws_config, bucket_name
            )

            if out == false
                s3_key = "solver_state/electric_fields_$(state.id).json.gz"
                s3_delete(aws_config, bucket_name, s3_key)
                return false
            end
            state.currentState = "gmres_solved"
        end


        if state.currentState == "gmres_solved"
            println("[$task_id] Salvataggio risultati...")
            out = state.intermediateResults["out"]
            resultsToStoreOnS3 = Dict(
                "Vp" => JSON.json(complex_matrix_to_float_array_matrix(out["Vp"])),
                "Ex" => JSON.json(complex_matrix_to_float_array_matrix(out["Ex"])),
                "Ey" => JSON.json(complex_matrix_to_float_array_matrix(out["Ey"])),
                "Ez" => JSON.json(complex_matrix_to_float_array_matrix(out["Ez"])),
                "Ex_3D" => JSON.json(complex_matrix_to_float_array_matrix(out["Ex_3D"])),
                "Ey_3D" => JSON.json(complex_matrix_to_float_array_matrix(out["Ey_3D"])),
                "Ez_3D" => JSON.json(complex_matrix_to_float_array_matrix(out["Ez_3D"])),
                "Hx_3D" => JSON.json(complex_matrix_to_float_array_matrix(out["Hx_3D"])),
                "Hy_3D" => JSON.json(complex_matrix_to_float_array_matrix(out["Hy_3D"])),
                "Hz_3D" => JSON.json(complex_matrix_to_float_array_matrix(out["Hz_3D"])),
                "centri_oss_3D" => JSON.json(transpose(state.intermediateResults["centri_oss_3D"])),
                "distanze_3D" => JSON.json(state.intermediateResults["distanze_3D"]),
                "theta_vals" => JSON.json(state.intermediateResults["theta_vals"]),
                "x_grid" => JSON.json(transpose(state.intermediateResults["x_grid"])),
                "y_grid" => JSON.json(transpose(state.intermediateResults["y_grid"])),
                "z_grid" => JSON.json(transpose(state.intermediateResults["z_grid"])),
                "baricentro" => JSON.json(state.baricentro),
                "f" => JSON.json(out["f"])
            )
            filename = task_id * "_results.json.gz"
            saveOnS3GZippedResults(task_id, resultsToStoreOnS3, aws_config, bucket_name)
            if !isnothing(chan)
                publish_data(Dict("computation_completed" => true, "path" => filename, "id" => task_id), "solver_feedback", chan)
            end
            state.currentState = "finished"
            s3_key = "solver_state/electric_fields_$(state.id).json.gz"
            s3_delete(aws_config, bucket_name, s3_key)
        end
    catch e
        println("[$task_id] Errore durante l'elaborazione")
        #state.currentState = "failed"
        #save_task_state(state, aws_config, bucket_name, "electric_fields")
        if e isa OutOfMemoryError
            publish_data(Dict("error" => "out of memory", "id" => task_id, "isStopped" => false, "partial" => false), "solver_feedback", chan)
        elseif e isa AMQPClient.AMQPClientException
            connection(; virtualhost=VIRTUALHOST, host=HOST) do conn
                AMQPClient.channel(conn, AMQPClient.UNUSED_CHANNEL, true) do chan
                    publish_data(Dict("error" => "broker connection closed", "id" => task_id, "isStopped" => false, "partial" => false), "solver_feedback", chan)
                    stopSolver(VIRTUALHOST, HOST)
                end
            end
        else
            error_msg = sprint(showerror, e)
            st = sprint((io, v) -> show(io, "text/plain", v), stacktrace(catch_backtrace()))
            @warn "Trouble doing things:\n$(error_msg)\n$(st)"
            publish_data(Dict("error" => "error", "id" => task_id, "isStopped" => false, "partial" => false), "solver_feedback", chan)
        end
    end
end

function save_task_state(state::ElectricFieldsSolvingState, aws_config, bucket_name, task_type="electric_fields")
    s3_key = "solver_state/$(task_type)_$(state.id).json.gz"
    try
        upload_state_json_gz(aws_config, bucket_name, s3_key, serialize_with_symbols_as_strings_recursive(state_to_dict(state)))
        println("[$(state.id)] Stato del task $task_type salvato su S3: s3://$bucket_name/$s3_key")
    catch e
        println("Errore nel salvataggio dello stato per $(state.id) ($task_type): $e")
    end
end

function getSignalbasedOn(signal_type, times)
    if signal_type == "exponential"
        f0=0
        dev_stand=0;

        time_delay_vs=3e-9;
        tw = 50*0.1/3e8;
        power=4.0;
        vs, tr = genera_segnale_esponenziale(tw, power, times, time_delay_vs)
        return vs
    elseif signal_type == "gaussian_modulated"
        tr=0;
        power=0;
        
        time_delay_vs=3e-9;
        f0=1e9;
        dev_stand=10/(4*pi)*1e-9;
        vs = genera_segnale_Gaussiano_modulato(f0, dev_stand, times, time_delay_vs)
        return vs
    elseif signal_type == "sinusoidal"
        tr=0;
        dev_stand=0;
        power=0;
        time_delay_vs=3e-9;
        f0=1e8;
        vs = genera_segnale_sinusoidale(f0, times, time_delay_vs)
        return vs
    elseif signal_type == "trapezoidal_pulse"
        Amplitude=1e-9;
        initial_delay_time=2e-9;
        high_level_time=10e-9;
        raise_time=10e-9;
        falling_time=10e-9;
        vs = build_trapezoidal_pulse(initial_delay_time, Amplitude, high_level_time, raise_time, falling_time, times)
        return vs
    end
end

function complex_matrix_to_float_array_matrix(complex_matrix::Matrix{ComplexF64})
    rows, cols = size(complex_matrix)
    float_array_matrix = Array{Array{Float64, 1}, 2}(undef, rows, cols)
    for i in 1:rows
        for j in 1:cols
            float_array_matrix[i, j] = [real(complex_matrix[i, j]), imag(complex_matrix[i, j])]
        end
    end
    return float_array_matrix
end