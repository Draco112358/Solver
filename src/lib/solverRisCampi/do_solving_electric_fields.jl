include("genera_segnale_esponenziale.jl")
include("genera_segnale_Gaussiano_modulato.jl")
include("genera_segnale_sinusoidale.jl")
include("get_signal.jl")
include("crea_freqs.jl")
include("build_trapezoidal_pulse.jl")
include("compute_fields_components.jl")
include("computeVs.jl")
include("fft_UAq.jl")
include("genera_punti_circonferenza.jl")
include("get_punti_oss_3D.jl")
include("iter_solver_E_Gaussian_Is_type.jl")
include("complex_matrix_to_float_array_matrix.jl")
include("../format_input_output_solver_functions.jl")
include("../sharedRis/calcola_P.jl")
include("../sharedRis/calcola_Lp.jl")
include("../sharedRis/find_nodes_ports_or_le.jl")
include("../utility.jl")

using MKL
using JSON


function doSolvingElectricFields(incidence_selection, volumi, superfici, nodi_coord, escalings, solverInput, solverAlgoParams, solverType, theta, phi, e_theta, e_phi, baricentro, r_circ, times, signal_type_E, ind_freq_interest, id, aws_config, bucket_name; chan=nothing, commentsEnabled=true)
    #@save "./test/electricFieldsSimulation/test_input_ris6x6.jld2" incidence_selection volumi superfici nodi_coord escalings solverInput solverAlgoParams solverType theta phi e_theta e_phi baricentro r_circ times signal_type_E ind_freq_interest id aws_config bucket_name
    try
        N_circ = 360
        N_circ_3D = 10
        time_delay_vs = 0.0
        f0=0.0
        dev_stand=0.0
        tr=0.0
        power=0.0
        if signal_type_E["type"] == "exponential"
            f0=0
            dev_stand=0;
            time_delay_vs=parse(Float64,signal_type_E["params"]["time_delay_vs"]);
            tw = parse(Float64,signal_type_E["params"]["tw"]);
            power=parse(Float64,signal_type_E["params"]["power"]);
            vs, tr = genera_segnale_esponenziale(tw, power, times, time_delay_vs)
        elseif signal_type_E["type"] == "gaussian_modulated"
            tr=0;
            power=0;
            time_delay_vs=parse(Float64,signal_type_E["params"]["time_delay_vs"]);
            f0=parse(Float64,signal_type_E["params"]["f0"]);
            dev_stand=parse(Float64,signal_type_E["params"]["dev_stand"]);
            vs = genera_segnale_Gaussiano_modulato(f0, dev_stand, times, time_delay_vs)
        elseif signal_type_E["type"] == "sinusoidal"
            tr=0;
            dev_stand=0;
            power=0;
            time_delay_vs=parse(Float64,signal_type_E["params"]["time_delay_vs"]);
            f0=parse(Float64,signal_type_E["params"]["f0"]);
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

        println("reading ports")
        nodi_coord = round.(nodi_coord, digits=8)
        ports, lumped_elements = find_nodes_ports_or_le(inputDict["ports"], inputDict["lumped_elements"], nodi_coord, escal)
        println("reading ports completed")
        GMRES_settings = Dict("Inner_Iter" => solverAlgoParams["innerIteration"], "Outer_Iter" => solverAlgoParams["outerIteration"], "tol" => solverAlgoParams["convergenceThreshold"] * ones((n_freq)))
        ind_low_freq = findall(x -> x < 1e5, freq)
        GMRES_settings["tol"][ind_low_freq] .= 1e-8
        QS_Rcc_FW = solverType # 1 QS, 2 Rcc, 3 Taylor
        use_Zs_in = false

        is_matrix = zeros(size(ports[:port_start], 1), length(times))
        for (index, signal_type) in enumerate(ports[:signals_port])
            if signal_type["type"] != "no_signal"
                vs = getSignalbasedOn(signal_type, times)
                is_matrix[index, :] = vs
            end
        end
        #is=getSignalbasedOn()

        println("P and Lp")
        P_data = @time calcola_P(superfici, escalings, QS_Rcc_FW, id)
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
        Lp_data = @time calcola_Lp(volumi, incidence_selection, escalings, QS_Rcc_FW, id)
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
        if is_stop_requested(id)
            println("Simulazione $(id) interrotta per richiesta stop.")
            return nothing # O un altro valore che indica interruzione
        end
        row_indices, col_indices, nz_values = findnz(incidence_selection[:A])
        A = sparse(row_indices, col_indices, nz_values)
        E,K,H,E_theta_v,E_phi_v = compute_fields_components(theta, phi, e_theta, e_phi)
        Vs= @time  computeVs(times,time_delay_vs,signal_type_E,volumi,nodi_coord,E,K,A,tr,power,dev_stand,f0,ind_freq_interest);

        Is = zeros(ComplexF64, size(ports[:port_nodes], 1), n_freq)
        for k in 1:size(ports[:port_nodes], 1)
            #Trasformata=fft_UAq(times, is_matrix[k, :])
            #Is[k,:]=Trasformata[2, ind_freq_interest]
            #Is[k,:].=0.02+0im
            # println(angle.(Is))
        end
        Is[1,:].=0.02+0im
        r_circ = r_circ*escal
        baricentro = baricentro .* escal
        punti_xy=genera_punti_circonferenza(r_circ,N_circ,baricentro,1);
        punti_zx=genera_punti_circonferenza(r_circ,N_circ,baricentro,2);
        punti_yz=genera_punti_circonferenza(r_circ,N_circ,baricentro,3);

        centri_oss=[punti_xy;punti_zx;punti_yz];
        centri_oss_3D, distanze_3D, theta_vals, x_grid, y_grid, z_grid = @time get_punti_oss_3D(r_circ, N_circ_3D, baricentro);
        println("gmres")
        out = @time iter_solver_E_Gaussian_Is_type(
            freq, escalings, incidence_selection, P_data, Lp_data,
                ports, lumped_elements, GMRES_settings, volumi, superfici, use_Zs_in, QS_Rcc_FW, ports_scatter_value, Vs, Is, centri_oss, centri_oss_3D, id, chan, commentsEnabled
        )
        if (isnothing(out))
            return nothing
        end
        # if test == true
        #     Ex = out["Ex"]
        #     Ey = out["Ey"]
        #     Ez = out["Ez"]
        #     @save "Ex.jld2" Ex
        #     @save "Ey.jld2" Ey
        #     @save "Ez.jld2" Ez
        # end
        # open("Ex_.txt", "w") do io
        # 	JSON.print(io, complex_matrix_to_float_array_matrix(out["Ex"]))
    	# end
		# open("Ey_.txt", "w") do io
        # 	JSON.print(io, complex_matrix_to_float_array_matrix(out["Ey"]))
    	# end
		# open("Ez_.txt", "w") do io
        # 	JSON.print(io, complex_matrix_to_float_array_matrix(out["Ez"]))
    	# end

        # if is_stopped_computation(id, chan)
        #     return false
        # end
        if (commentsEnabled == true)
            println("data publish")
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
                "centri_oss_3D" => JSON.json(transpose(centri_oss_3D)),
                "distanze_3D" => JSON.json(distanze_3D),
                "theta_vals" => JSON.json(theta_vals),
                "x_grid" => JSON.json(transpose(x_grid)),
                "y_grid" => JSON.json(transpose(y_grid)),
                "z_grid" => JSON.json(transpose(z_grid)),
                "baricentro" => JSON.json(baricentro),
                "f" => JSON.json(out["f"])
            )
            # publish_data(dataToReturn, "solver_results", chan)
            filename = id * "_results.json.gz"
            saveOnS3GZippedResults(id, resultsToStoreOnS3, aws_config, bucket_name)
            send_rabbitmq_feedback(Dict("computation_completed" => true, "path" => filename, "id" => id), "solver_feedback")
            #s3_put(aws_config, bucket_name, filename, JSON.json(resultsToStoreOnS3))
            # if !isnothing(chan)
            #     publish_data(Dict("computation_completed" => true, "path" => filename, "id" => id), "solver_feedback", chan)
            # end
        end
        return out
    catch e
        if e isa OutOfMemoryError
            send_rabbitmq_feedback(Dict("error" => "out of memory", "id" => id, isStopped => false, partial: false), "solver_feedback")
        else
            error_msg = sprint(showerror, e)
            st = sprint((io,v) -> show(io, "text/plain", v), stacktrace(catch_backtrace()))
            @warn "Trouble doing things:\n$(error_msg)\n$(st)"
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
