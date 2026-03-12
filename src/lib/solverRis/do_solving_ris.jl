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
        freq = Vector{Float64}(undef, length(frequencies))
        for i in range(1, length(frequencies))
            freq[i] = frequencies[i]
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
        # ----------------------------------------------------------------------
        # CHECKPOINT LOAD P Matrix
        # ----------------------------------------------------------------------
        checkpoint_P = load_checkpoint(id, "P_matrix")
        
        P_data = nothing
        Lp_data = nothing
        
        superfici["estremi_celle"] = hcat(superfici["estremi_celle"]...)
        superfici["centri"] = hcat(superfici["centri"]...)
        
        # --- P Matrix Checkpoint Logic ---
        if !isnothing(checkpoint_P) && haskey(checkpoint_P, :P) && haskey(checkpoint_P, :R_cc) && !haskey(checkpoint_P, :P_start_block)
            println("P_data fully recovered from checkpoint.")
            P_data = Dict{Symbol, Any}(:P => checkpoint_P[:P], :R_cc => checkpoint_P[:R_cc])
            send_rabbitmq_feedback(Dict("computingP" => true, "id" => id), "solver_feedback")
        else
            println("P and Lp")
            
            # Use partial checkpointing if block exists
            start_block_P = (!isnothing(checkpoint_P) && haskey(checkpoint_P, :P_start_block)) ? checkpoint_P[:P_start_block] : 1
            
            function checkpoint_cb_P(P_mat, R_cc_mat, next_block::Int)
                state_P = Dict{Symbol, Any}(:P => P_mat, :R_cc => R_cc_mat)
                if next_block != typemax(Int)
                    state_P[:P_start_block] = next_block
                end
                save_checkpoint(id, "P_matrix", state_P)
            end

            P_data = calcola_P(superfici, escalings, QS_Rcc_FW, id, start_block_P, checkpoint_cb_P, checkpoint_P)
            if isnothing(P_data)
                return nothing
            end
            send_rabbitmq_feedback(Dict("computingP" => true, "id" => id), "solver_feedback")
            
            # Save state immediately after P computes cleanly
            save_checkpoint(id, "P_matrix", Dict{Symbol, Any}(:P => P_data[:P], :R_cc => P_data[:R_cc]))
        end
        
        # --- Lp Matrix Checkpoint Logic ---
        checkpoint_Lp = load_checkpoint(id, "Lp_matrix")
        
        if !isnothing(checkpoint_Lp) && (!haskey(checkpoint_Lp, :Lp_start_block)) && haskey(checkpoint_Lp, :Lp_x)
            println("Lp_data fully recovered from checkpoint.")
            Lp_data = Dict{Symbol, Any}(:Lp_x => checkpoint_Lp[:Lp_x], :Lp_y => checkpoint_Lp[:Lp_y], :Lp_z => checkpoint_Lp[:Lp_z], 
                           :Rx => checkpoint_Lp[:Rx], :Ry => checkpoint_Lp[:Ry], :Rz => checkpoint_Lp[:Rz])
            send_rabbitmq_feedback(Dict("computingLp" => true, "id" => id), "solver_feedback")
        else
            start_block_Lp = (!isnothing(checkpoint_Lp) && haskey(checkpoint_Lp, :Lp_start_block)) ? checkpoint_Lp[:Lp_start_block] : 1
            
            function checkpoint_cb_Lp(Lp_x_mat, Lp_y_mat, Lp_z_mat, Rx_mat, Ry_mat, Rz_mat, next_block::Int, curr_component::Int)
                state_Lp = Dict{Symbol, Any}(
                    :Lp_x => Lp_x_mat, :Lp_y => Lp_y_mat, :Lp_z => Lp_z_mat,
                    :Rx => Rx_mat, :Ry => Ry_mat, :Rz => Rz_mat
                )
                if next_block != typemax(Int)
                    state_Lp[:Lp_start_block] = next_block
                    state_Lp[:curr_component] = curr_component
                end
                save_checkpoint(id, "Lp_matrix", state_Lp)
            end
            
            Lp_data = calcola_Lp(volumi, incidence_selection, escalings, QS_Rcc_FW, id, start_block_Lp, checkpoint_cb_Lp, checkpoint_Lp)
            if isnothing(Lp_data)
                return nothing
            end
            send_rabbitmq_feedback(Dict("computingLp" => true, "id" => id), "solver_feedback")
            
            # Save state immediately after Lp computes cleanly
            save_checkpoint(id, "Lp_matrix", Dict{Symbol, Any}(
                :Lp_x => Lp_data[:Lp_x], :Lp_y => Lp_data[:Lp_y], :Lp_z => Lp_data[:Lp_z],
                :Rx => Lp_data[:Rx], :Ry => Lp_data[:Ry], :Rz => Lp_data[:Rz]
            ))
        end

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
        return out
    catch e
        if e isa OutOfMemoryError
            send_rabbitmq_feedback(Dict("error" => "out of memory", "id" => id, "isStopped" => false, "partial" => false), "solver_feedback")
            #publish_data(Dict("error" => "out of memory", "id" => id, isStopped => false, partial: false), "solver_feedback", chan)
        else
            error_msg = sprint(showerror, e)
            st = sprint((io, v) -> show(io, "text/plain", v), stacktrace(catch_backtrace()))
            @warn "Trouble doing things:\n$(error_msg)\n$(st)"
            send_rabbitmq_feedback(Dict("error" => "Internal Server Error", "id" => id, "isStopped" => false, "partial" => false), "solver_feedback")
        end
    finally
        # Nessuna pulizia necessaria: il meccanismo di stop è ora basato su processi
    end

end