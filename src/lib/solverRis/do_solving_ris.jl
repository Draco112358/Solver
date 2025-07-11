function doSolvingRis(incidence_selection, volumi, superfici, nodi_coord, escalings, solverInput, solverAlgoParams, solverType, id, aws_config, bucket_name; chan=nothing, commentsEnabled=true)
    try
        inputDict = solverInput
        unit = solverInput["unit"]
        escal = getEscalFrom(unit)
        #nodi_coord = hcat(nodi_coord...)
        ports_scatter_value = haskey(solverInput, "ports_scattering_value") ? solverInput["ports_scattering_value"] : 50.0
        # if is_stopped_computation(id, chan)
        #     return false
        # end

        frequencies = inputDict["frequencies"]
        freq = Array{Float64}(undef, 1, length(frequencies))
        for i in range(1, length(frequencies))
            freq[1, i] = frequencies[i]
        end
        n_freq = length(freq)
        println("reading ports")
        ports, lumped_elements = find_nodes_ports_or_le(inputDict["ports"], inputDict["lumped_elements"], nodi_coord, escal)
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
        P_data = calcola_P(superfici, escalings, QS_Rcc_FW, id)
        if isnothing(P_data)
            return nothing
        end
        send_rabbitmq_feedback(Dict("computingP" => true, "id" => id), "solver_feedback")
        # if !isnothing(chan)
        #     publish_data(Dict("computingP" => true, "id" => id), "solver_feedback", chan)
        # end
        # if is_stopped_computation(id, chan)
        #     return false
        # end
        Lp_data = calcola_Lp(volumi, incidence_selection, escalings, QS_Rcc_FW, id)
        if isnothing(Lp_data)
            return nothing
        end
        send_rabbitmq_feedback(Dict("computingLp" => true, "id" => id), "solver_feedback")
        # if !isnothing(chan)
        #     publish_data(Dict("computingLp" => true, "id" => id), "solver_feedback", chan)
        # end
        # if is_stopped_computation(id, chan)
        #     return false
        # end

        println("gmres")
        out = iter_solver_QS_S_type(
            freq, escalings, incidence_selection, P_data, Lp_data,
                ports, lumped_elements, GMRES_settings, volumi, use_Zs_in, QS_Rcc_FW, ports_scatter_value, id, chan, commentsEnabled
        )
        if (isnothing(out))
            return nothing
        end
        println("data publish")

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
            #publish_data(dataToReturn, "solver_results", chan)
            filename = id * "_results.json.gz"
            saveOnS3GZippedResults(id, resultsToStoreOnS3, aws_config, bucket_name)
            send_rabbitmq_feedback(dataToReturn, "solver_results")
            send_rabbitmq_feedback(Dict("computation_completed" => true, "path" => filename, "id" => id), "solver_feedback")
            #s3_put(aws_config, bucket_name, filename, JSON.json(resultsToStoreOnS3))
            # if !isnothing(chan)
            #     publish_data(Dict("computation_completed" => true, "path" => filename, "id" => id), "solver_feedback", chan)
            # end
        end
    catch e
        if e isa OutOfMemoryError
            send_rabbitmq_feedback(Dict("error" => "out of memory", "id" => id, isStopped => false, partial: false), "solver_feedback")
            #publish_data(Dict("error" => "out of memory", "id" => id, isStopped => false, partial: false), "solver_feedback", chan)
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