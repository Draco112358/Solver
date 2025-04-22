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
include("iter_solver_E_Gaussian_Is_type.jl")

using MKL
using JSON, AWSS3
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

function doSolving(mesherOutput, solverInput, solverAlgoParams, solverType, id, aws_config, bucket_name; chan=nothing, commentsEnabled=true)
    try
        mesherDict = mesherOutput
        inputDict = solverInput
        unit = solverInput["unit"]
        escal = getEscalFrom(unit)
        ports_scatter_value = haskey(solverInput, "ports_scattering_value") ? solverInput["ports_scattering_value"] : 50.0

        sx, sy, sz = mesherDict["cell_size"]["cell_size_x"] * escal, mesherDict["cell_size"]["cell_size_y"] * escal, mesherDict["cell_size"]["cell_size_z"] * escal

        origin = (mesherDict["origin"]["origin_x"], mesherDict["origin"]["origin_y"], mesherDict["origin"]["origin_z"])

        if is_stopped_computation(id, chan)
            return false
        end

        MATERIALS = [material(el) for el in inputDict["materials"]]

        dominant_list = []
        grids = []
        conductors_index = 1
        for (mat, value) in mesherDict["mesher_matrices"]
            push!(grids, unsqueeze(copy(value), dims=2))
            if (isMaterialConductor(mat, MATERIALS))
                push!(dominant_list, conductors_index)
            end
            conductors_index += 1
        end
        # testarray = [copy(value) for (index, value) in mesherDict["mesher_matrices"]]

        # grids = [unsqueeze(values, dims=2) for values in testarray]

        if is_stopped_computation(id, chan)
            return false
        end

        frequencies = inputDict["frequencies"]
        freq = Array{Float64}(undef, 1, length(frequencies))
        for i in range(1, length(frequencies))
            freq[1, i] = frequencies[i]
        end
        #freq = convert(Array{Float64}, freq)

        n_freq = length(freq)
        println("reading ports")

        PORTS = read_ports(inputDict["ports"], escal)

        L_ELEMENTS = read_lumped_elements(inputDict["lumped_elements"], escal)

        MATERIALS = [material(el) for el in inputDict["materials"]]
        println("reading ports completed")
        #SIGNALS = [el for el in inputDict["signals"]]

        if is_stopped_computation(id, chan)
            return false
        end

        # # START SETTINGS--------------------------------------------
        # ind_low_freq= filter(i -> !iszero(freq[i]), findall(f -> f<1e5, frequencies))
        # tol[ind_low_freq] .= 1e-7
        GMRES_settings = Dict("Inner_Iter" => solverAlgoParams["innerIteration"], "Outer_Iter" => solverAlgoParams["outerIteration"], "tol" => solverAlgoParams["convergenceThreshold"] * ones((n_freq)))
        QS_Rcc_FW = solverType # 1 QS, 2 Rcc, 3 Taylor
        use_escalings = 1
        println("create_volumes_mapping_v2")
        mapping_vols, num_centri = create_volumes_mapping_v2(grids)

        if is_stopped_computation(id, chan)
            return false
        end

        println("create_volume_centers")
        centri_vox, id_mat = create_volume_centers(grids, mapping_vols, num_centri, sx, sy, sz, origin)

        if is_stopped_computation(id, chan)
            return false
        end
        println("create_Grids_externals")
        externals_grids = create_Grids_externals(grids)

        if is_stopped_computation(id, chan)
            return false
        end

        println("MesherFFT")
        escalings, incidence_selection, circulant_centers, diagonals, expansions, ports, lumped_elements, li_mats, Zs_info = mesher_FFT(use_escalings, MATERIALS, sx, sy, sz, grids, centri_vox, externals_grids, mapping_vols, PORTS, L_ELEMENTS, origin, commentsEnabled, dominant_list)
        # if length(stopComputation) > 0
        #     pop!(stopComputation)
        #     return Dict("id" => id, "isStopped" => true)
        # end

        if is_stopped_computation(id, chan)
            return false
        end

        println("P and Lp")
        FFTCP, FFTCLp = compute_FFT_mutual_coupling_mats(circulant_centers, escalings, Int64(mesherDict["n_cells"]["n_cells_x"]), Int64(mesherDict["n_cells"]["n_cells_y"]), Int64(mesherDict["n_cells"]["n_cells_z"]), QS_Rcc_FW, id, chan)

        if is_stopped_computation(id, chan)
            return false
        end

        #@profile FFT_solver_QS_S_type(freq,escalings,incidence_selection,FFTCP,FFTCLp,diagonals,ports,lumped_elements,expansions,GMRES_settings,Zs_info,QS_Rcc_FW);
        # if length(stopComputation) > 0
        #     pop!(stopComputation)
        #     return Dict("id" => id, "isStopped" => true)
        # end
        println("gmres")
        out = @time FFT_solver_QS_S_type(freq, escalings, incidence_selection, FFTCP, FFTCLp, diagonals, ports, ports_scatter_value, lumped_elements, expansions, GMRES_settings, Zs_info, QS_Rcc_FW, id, chan, commentsEnabled)
        println("data publish")
        if (out == false)
            return false
        end

        if is_stopped_computation(id, chan)
            return false
        end
        if (commentsEnabled == true)
            #publish_data(dump_json_data(out["Z"], out["S"], out["Y"], length(inputDict["ports"]), id), "solver_results", chan)
            filename = id * "_results.json.gz"
            resultsToStoreOnS3 = dump_json_data(out["Z"], out["S"], out["Y"], length(inputDict["ports"]), id)
            #s3_put(aws_config, bucket_name, filename, JSON.json(resultsToStoreOnS3))
            dataToReturn = Dict(
                "portIndex" => 0,
                "partial" => false,
                "results" => Dict(
                    "matrixZ" => JSON.parse(resultsToStoreOnS3["matrices"]["matrix_Z"])[1],
                    "matrixS" => JSON.parse(resultsToStoreOnS3["matrices"]["matrix_S"])[1],
                    "matrixY" => JSON.parse(resultsToStoreOnS3["matrices"]["matrix_Y"])[1],
                )
            )
            publish_data(dataToReturn, "solver_results", chan)
            saveOnS3GZippedResults(id, resultsToStoreOnS3, aws_config, bucket_name)
            if !isnothing(chan)
                publish_data(Dict("computation_completed" => true, "path" => filename, "id" => id), "solver_feedback", chan)
            end
        end
    catch e
        if e isa OutOfMemoryError
            publish_data(Dict("error" => "out of memory", "id" => id, isStopped => false, partial: false), "solver_feedback", chan)
        else
            error_msg = sprint(showerror, e)
            st = sprint((io,v) -> show(io, "text/plain", v), stacktrace(catch_backtrace()))
            @warn "Trouble doing things:\n$(error_msg)\n$(st)"
        end
    end

end

function doSolvingRis(incidence_selection, volumi, superfici, nodi_coord, escalings, solverInput, solverAlgoParams, solverType, id, aws_config, bucket_name; chan=nothing, commentsEnabled=true)
    try
        inputDict = solverInput
        unit = solverInput["unit"]
        escal = getEscalFrom(unit)
        #nodi_coord = hcat(nodi_coord...)
        ports_scatter_value = haskey(solverInput, "ports_scattering_value") ? solverInput["ports_scattering_value"] : 50.0
        println("scatter ", ports_scatter_value)
        # if is_stopped_computation(id, chan)
        #     return false
        # end

        freq = inputDict["frequencies"]

        n_freq = length(freq)
        println("reading ports")
        println(size(nodi_coord))
        ports, lumped_elements = find_nodes_ports_or_le(inputDict["ports"], inputDict["lumped_elements"], nodi_coord, escal)
        println("port_nodes ", ports[:port_nodes])
        println("lumped_nodes ", lumped_elements[:le_nodes])
        println("reading ports completed")
        #SIGNALS = [el for el in inputDict["signals"]]

        # if is_stopped_computation(id, chan)
        #     return falses
        # end

        # # START SETTINGS--------------------------------------------
        # ind_low_freq= filter(i -> !iszero(freq[i]), findall(f -> f<1e5, frequencies))
        # tol[ind_low_freq] .= 1e-7
        GMRES_settings = Dict("Inner_Iter" => solverAlgoParams["innerIteration"], "Outer_Iter" => solverAlgoParams["outerIteration"], "tol" => solverAlgoParams["convergenceThreshold"] * ones((n_freq)))
        ind_low_freq = findall(x -> x < 1e5, freq)
        GMRES_settings["tol"][ind_low_freq] .= 1e-8
        QS_Rcc_FW = solverType # 1 QS, 2 Rcc, 3 Taylor
        use_Zs_in = true

        println("P and Lp")
        P_data = calcola_P(superfici, escalings, QS_Rcc_FW)
        if !isnothing(chan)
            publish_data(Dict("computingP" => true, "id" => id), "solver_feedback", chan)
        end
        Lp_data = calcola_Lp2(volumi, incidence_selection, escalings, QS_Rcc_FW)
        if !isnothing(chan)
            publish_data(Dict("computingLp" => true, "id" => id), "solver_feedback", chan)
        end
        # if is_stopped_computation(id, chan)
        #     return false
        # end

        println("gmres")
        out = iter_solver_QS_S_type2(
            freq, escalings, incidence_selection, P_data, Lp_data,
                ports, lumped_elements, GMRES_settings, volumi, use_Zs_in, QS_Rcc_FW, ports_scatter_value, id, chan, commentsEnabled
        )
        println("data publish")
        if (out == false)
            return false
        end

        # if is_stopped_computation(id, chan)
        #     return false
        # end
        if (commentsEnabled == true)
            resultsToStoreOnS3 = dump_json_data(out[:Z], out[:S], out[:Y], length(inputDict["ports"]), id)
            dataToReturn = Dict(
                "portIndex" => 0,
                "partial" => false,
                "results" => Dict(
                    "matrixZ" => JSON.parse(resultsToStoreOnS3["matrices"]["matrix_Z"])[1],
                    "matrixS" => JSON.parse(resultsToStoreOnS3["matrices"]["matrix_S"])[1],
                    "matrixY" => JSON.parse(resultsToStoreOnS3["matrices"]["matrix_Y"])[1],
                )
            )
            publish_data(dataToReturn, "solver_results", chan)
            filename = id * "_results.json.gz"
            saveOnS3GZippedResults(id, resultsToStoreOnS3, aws_config, bucket_name)
            #s3_put(aws_config, bucket_name, filename, JSON.json(resultsToStoreOnS3))
            if !isnothing(chan)
                publish_data(Dict("computation_completed" => true, "path" => filename, "id" => id), "solver_feedback", chan)
            end
        end
    catch e
        if e isa OutOfMemoryError
            publish_data(Dict("error" => "out of memory", "id" => id, isStopped => false, partial: false), "solver_feedback", chan)
        else
            error_msg = sprint(showerror, e)
            st = sprint((io,v) -> show(io, "text/plain", v), stacktrace(catch_backtrace()))
            @warn "Trouble doing things:\n$(error_msg)\n$(st)"
        end
    end

end

function doSolvingElectricFields(incidence_selection, volumi, superfici, nodi_coord, escalings, solverInput, solverAlgoParams, solverType, theta, phi, e_theta, e_phi, baricentro, r_circ, times, signal_type_E, ind_freq_interest, id, aws_config, bucket_name; chan=nothing, commentsEnabled=true)
    try
        N_circ = 100
        N_circ_3D = 50
        time_delay_vs = 0.0
        f0=0.0
        dev_stand=0.0
        tr=0.0
        power=0.0
        if signal_type_E == "exponential"
            f0=0
            dev_stand=0;
    
            time_delay_vs=3e-9;
            tw = 50*0.1/3e8;
            power=4.0;
            vs, tr = genera_segnale_esponenziale(tw, power, times, time_delay_vs)
        elseif signal_type_E == "gaussian_modulated"
            tr=0;
            power=0;
            
            time_delay_vs=3e-9;
            f0=1e9;
            dev_stand=10/(4*pi)*1e-9;
            vs = genera_segnale_Gaussiano_modulato(f0, dev_stand, times, time_delay_vs)
        elseif signal_type_E == "sinusoidal"
            tr=0;
            dev_stand=0;
            power=0;
            time_delay_vs=3e-9;
            f0=1e8;
            vs = genera_segnale_sinusoidale(f0, times, time_delay_vs)
        end
        ind_freq_interest = ind_freq_interest .+ 1
        freq_vs_is=crea_freqs(times);
        freq = freq_vs_is[ind_freq_interest];
        n_freq = length(freq)
        inputDict = solverInput
        unit = solverInput["unit"]
        escal = getEscalFrom(unit)
        ports_scatter_value = haskey(solverInput, "ports_scattering_value") ? solverInput["ports_scattering_value"] : 50.0
        println("scatter ", ports_scatter_value)

        println("reading ports")
        println(size(nodi_coord))
        ports, lumped_elements = find_nodes_ports_or_le(inputDict["ports"], inputDict["lumped_elements"], nodi_coord, escal)
        println("port_nodes ", ports[:port_nodes])
        println("lumped_nodes ", lumped_elements[:le_nodes])
        println("reading ports completed")
        
        GMRES_settings = Dict("Inner_Iter" => solverAlgoParams["innerIteration"], "Outer_Iter" => solverAlgoParams["outerIteration"], "tol" => solverAlgoParams["convergenceThreshold"] * ones((n_freq)))
        ind_low_freq = findall(x -> x < 1e5, freq)
        GMRES_settings["tol"][ind_low_freq] .= 1e-8
        QS_Rcc_FW = solverType # 1 QS, 2 Rcc, 3 Taylor
        use_Zs_in = true

        is_matrix = zeros(size(ports[:port_start], 1), length(times))
        for (index, signal_type) in enumerate(ports[:signals_port])
            if signal_type != "no_signal"
                vs = getSignalbasedOn(signal_type, times)
                is_matrix[index, :] = vs
            end
        end
        #is=getSignalbasedOn()

        println("P and Lp")
        P_data = calcola_P(superfici, escalings, QS_Rcc_FW)
        if !isnothing(chan)
            publish_data(Dict("computingP" => true, "id" => id), "solver_feedback", chan)
        end
        Lp_data = calcola_Lp2(volumi, incidence_selection, escalings, QS_Rcc_FW)
        if !isnothing(chan)
            publish_data(Dict("computingLp" => true, "id" => id), "solver_feedback", chan)
        end
        # if is_stopped_computation(id, chan)
        #     return false
        # end
        row_indices, col_indices, nz_values = findnz(incidence_selection[:A])
        A = sparse(row_indices, col_indices, nz_values)
        E,K,H,E_theta_v,E_phi_v = compute_fields_components(theta, phi, e_theta, e_phi)
        Vs=computeVs(times,time_delay_vs,signal_type_E,volumi,nodi_coord,E,K,A,tr,power,dev_stand,f0,ind_freq_interest);

        Is = zeros(ComplexF64, size(ports[:port_nodes], 1), n_freq)
        for k in 1:size(ports[:port_nodes], 1)
            Trasformata=fft_UAq(times, is_matrix[k, :])
            Is[k,:]=Trasformata[2, ind_freq_interest]
        end

        punti_xy=genera_punti_circonferenza(r_circ,N_circ,baricentro,1);
        punti_zx=genera_punti_circonferenza(r_circ,N_circ,baricentro,2);
        punti_yz=genera_punti_circonferenza(r_circ,N_circ,baricentro,3);

        centri_oss=[punti_xy;punti_zx;punti_yz];
        centri_oss_3D, distanze_3D, theta_vals, x_grid, y_grid, z_grid = get_punti_oss_3D(r_circ, N_circ_3D, baricentro);

        println("gmres")
        out = iter_solver_E_Gaussian_Is_type(
            freq, escalings, incidence_selection, P_data, Lp_data,
                ports, lumped_elements, GMRES_settings, volumi, superfici, use_Zs_in, QS_Rcc_FW, ports_scatter_value, Vs, Is, centri_oss, centri_oss_3D, id, chan, commentsEnabled
        )
        println("data publish")
        if (out == false)
            return false
        end

        println(out)

        # if is_stopped_computation(id, chan)
        #     return false
        # end
        if (commentsEnabled == true)
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
                "f" => JSON.json(out["f"])
            )
            # publish_data(dataToReturn, "solver_results", chan)
            filename = id * "_results.json.gz"
            saveOnS3GZippedResults(id, resultsToStoreOnS3, aws_config, bucket_name)
            #s3_put(aws_config, bucket_name, filename, JSON.json(resultsToStoreOnS3))
            if !isnothing(chan)
                publish_data(Dict("computation_completed" => true, "path" => filename, "id" => id), "solver_feedback", chan)
            end
        end
    catch e
        if e isa OutOfMemoryError
            publish_data(Dict("error" => "out of memory", "id" => id, isStopped => false, partial: false), "solver_feedback", chan)
        else
            error_msg = sprint(showerror, e)
            st = sprint((io,v) -> show(io, "text/plain", v), stacktrace(catch_backtrace()))
            @warn "Trouble doing things:\n$(error_msg)\n$(st)"
        end
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

# DotEnv.load!()

# aws_access_key_id = ENV["AWS_ACCESS_KEY_ID"]
# aws_secret_access_key = ENV["AWS_SECRET_ACCESS_KEY"]
# aws_region = ENV["AWS_DEFAULT_REGION"]
# aws_bucket_name = ENV["AWS_BUCKET_NAME"]
# creds = AWSCredentials(aws_access_key_id, aws_secret_access_key)
# aws = global_aws_config(; region=aws_region, creds=creds)
# mesherOutput = download_json_gz(aws, aws_bucket_name, "417782681790578896_mesh.json.gz")
# surface = download_json_gz(aws, aws_bucket_name, "417782681790578896_surface.json.gz")
# data = Dict{String, Any}("solverAlgoParams" => Dict{String, Any}("innerIteration" => 100, "convergenceThreshold" => 0.0001, "outerIteration" => 1), "mesherFileId" => "417778305446445264_mesh.json.gz", "id" => "417778305446445264", "storage" => "local", "solverType" => 2, "solverInput" => Dict{String, Any}("lumped_elements" => Any[Dict{String, Any}("isSelected" => false, "name" => "lumped-1", "outputElement" => Any[1.5, 0, 1.05], "inputElement" => Any[1.5, 0, 0.05], "value" => 0, "category" => "lumped", "type" => 1, "rlcParams" => Dict{String, Any}("capacitance" => 0, "inductance" => 0, "resistance" => 50)), Dict{String, Any}("isSelected" => true, "name" => "lumped-2", "outputElement" => Any[1.5, 5, 1.05], "inputElement" => Any[1.5, 5, 0.05], "value" => 0, "category" => "lumped", "type" => 1, "rlcParams" => Dict{String, Any}("capacitance" => 0, "inductance" => 0, "resistance" => 50))], "materials" => Any[Dict{String, Any}("name" => "antennaMaterial", "permeability" => 1, "coll" => Dict{String, Any}("name" => "Materials"), "id" => "408755738118193360", "color" => "#f8b054", "conductivity" => 58000000, "permittivity" => 1, "ts" => Dict{String, Any}("isoString" => "2024-09-11T18:18:19.150Z")), Dict{String, Any}("name" => "antennaDielMaterial", "permeability" => 1, "coll" => Dict{String, Any}("name" => "Materials"), "id" => "408755797380563147", "color" => "#bcbcbc", "conductivity" => 0, "permittivity" => 5, "ts" => Dict{String, Any}("isoString" => "2024-09-11T18:19:15.650Z"))], "ports_scattering_value" => 50, "unit" => "mm", "ports" => Any[Dict{String, Any}("isSelected" => false, "name" => "port-1", "outputElement" => Any[1.5, 0, 1.05], "inputElement" => Any[1.5, 0, 0.05], "category" => "port"), Dict{String, Any}("isSelected" => false, "name" => "port-2", "outputElement" => Any[1.5, 5, 1.05], "inputElement" => Any[1.5, 5, 0.05], "category" => "port")], "frequencies" => Any[100, 316.2277660168379, 1000, 3162.2776601683795, 10000, 31622.776601683792, 100000, 316227.7660168379, 1000000, 3.1622776601683795e6, 10000000, 3.162277660168379e7, 100000000, 3.1622776601683795e8, 1000000000]))
# doSolvingRis(mesherOutput["incidence_selection"], mesherOutput["volumi"], surface, mesherOutput["nodi_coord"], mesherOutput["escalings"], data["solverInput"], data["solverAlgoParams"], 2, "417782681790578896", aws, aws_bucket_name)