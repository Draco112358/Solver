function doSolvingFFT(mesherOutput, solverInput, solverAlgoParams, solverType, id, aws_config, bucket_name; chan=nothing, commentsEnabled=true)
    try
        mesherDict = mesherOutput
        inputDict = solverInput
        unit = solverInput["unit"]
        escal = getEscalFrom(unit)
        ports_scatter_value = haskey(solverInput, "ports_scattering_value") ? solverInput["ports_scattering_value"] : 50.0

        sx, sy, sz = mesherDict["cell_size"]["cell_size_x"] * escal, mesherDict["cell_size"]["cell_size_y"] * escal, mesherDict["cell_size"]["cell_size_z"] * escal

        origin = (mesherDict["origin"]["origin_x"], mesherDict["origin"]["origin_y"], mesherDict["origin"]["origin_z"])

        # if is_stopped_computation(id, chan)
        #     return false
        # end

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

        # if is_stopped_computation(id, chan)
        #     return false
        # end

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
        println("create_volumes_mapping")
        mapping_vols, num_centri = create_volumes_mapping(grids)

        # if is_stopped_computation(id, chan)
        #     return false
        # end
        if is_stop_requested(id)
            println("Simulazione $(id) interrotta per richiesta stop.")
            return nothing # O un altro valore che indica interruzione
        end

        println("create_volume_centers")
        centri_vox, id_mat = create_volume_centers(grids, mapping_vols, num_centri, sx, sy, sz, origin)

        # if is_stopped_computation(id, chan)
        #     return false
        # end
        if is_stop_requested(id)
            println("Simulazione $(id) interrotta per richiesta stop.")
            return nothing # O un altro valore che indica interruzione
        end
        println("create_Grids_externals")
        externals_grids = create_Grids_externals(grids)

        # if is_stopped_computation(id, chan)
        #     return false
        # end
        if is_stop_requested(id)
            println("Simulazione $(id) interrotta per richiesta stop.")
            return nothing # O un altro valore che indica interruzione
        end

        println("MesherFFT")
        escalings, incidence_selection, circulant_centers, diagonals, expansions, ports, lumped_elements, li_mats, Zs_info = mesher_FFT(use_escalings, MATERIALS, sx, sy, sz, grids, centri_vox, externals_grids, mapping_vols, PORTS, L_ELEMENTS, origin, commentsEnabled, dominant_list, id)
        if isnothing(escalings)
            return nothing
        end
        # if length(stopComputation) > 0
        #     pop!(stopComputation)
        #     return Dict("id" => id, "isStopped" => true)
        # end

        # if is_stopped_computation(id, chan)
        #     return false
        # end
        if is_stop_requested(id)
            println("Simulazione $(id) interrotta per richiesta stop.")
            return nothing # O un altro valore che indica interruzione
        end

        println("P and Lp")
        FFTCP, FFTCLp = compute_FFT_mutual_coupling_mats(circulant_centers, escalings, Int64(mesherDict["n_cells"]["n_cells_x"]), Int64(mesherDict["n_cells"]["n_cells_y"]), Int64(mesherDict["n_cells"]["n_cells_z"]), QS_Rcc_FW, id, chan)
        if isnothing(FFTCP)
            return nothing
        end
        if isnothing(FFTCLp)
            return nothing
        end
        # if is_stopped_computation(id, chan)
        #     return false
        # end

        #@profile FFT_solver_QS_S_type(freq,escalings,incidence_selection,FFTCP,FFTCLp,diagonals,ports,lumped_elements,expansions,GMRES_settings,Zs_info,QS_Rcc_FW);
        # if length(stopComputation) > 0
        #     pop!(stopComputation)
        #     return Dict("id" => id, "isStopped" => true)
        # end
        println("gmres")
        out = @time FFT_solver_QS_S_type(freq, escalings, incidence_selection, FFTCP, FFTCLp, diagonals, ports, ports_scatter_value, lumped_elements, expansions, GMRES_settings, Zs_info, QS_Rcc_FW, id, chan, commentsEnabled)
        println("data publish")
        if isnothing(out)
            return nothing
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
            send_rabbitmq_feedback(dataToReturn, "solver_results")
            saveOnS3GZippedResults(id, resultsToStoreOnS3, aws_config, bucket_name)
            send_rabbitmq_feedback(Dict("computation_completed" => true, "path" => filename, "id" => id), "solver_feedback")
        end
    catch e
        if e isa OutOfMemoryError
            send_rabbitmq_feedback(Dict("error" => "out of memory", "id" => id, isStopped => false, partial: false), "solver_feedback")
        else
            error_msg = sprint(showerror, e)
            st = sprint((io,v) -> show(io, "text/plain", v), stacktrace(catch_backtrace()))
            @warn "Trouble doing things:\n$(error_msg)\n$(st)"
            send_rabbitmq_feedback(Dict("error" => "Internal Server Error", "id" => id, isStopped => false, partial: false), "solver_feedback")
        end
    finally
        # Pulizia del flag di stop indipendentemente da come la simulazione finisce
        lock(stop_computation_lock) do
            if haskey(stopComputation, id)
                delete!(stopComputation, id)
                println("Flag di stop per simulazione $(id) rimosso.")
            end
        end
    end

end
