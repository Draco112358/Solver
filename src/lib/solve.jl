include("FFT_solver_QS_S_type.jl")
include("create_volume_centers.jl")
include("create_Grids_externals.jl")
include("compute_FFT_mutual_coupling_mats.jl")
include("mesher_FFT.jl")
include("From_3D_to_1D.jl")
include("utility.jl")

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
            filename = id * "_results.json"
            resultsToStoreOnS3 = dump_json_data(out["Z"], out["S"], out["Y"], length(inputDict["ports"]), id)
            s3_put(aws_config, bucket_name, filename, JSON.json(resultsToStoreOnS3))
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
