

aws_access_key_id = ENV["AWS_ACCESS_KEY_ID"]
aws_secret_access_key = ENV["AWS_SECRET_ACCESS_KEY"]
aws_region = ENV["AWS_DEFAULT_REGION"]
aws_bucket_name = ENV["AWS_BUCKET_NAME"]
creds = AWSCredentials(aws_access_key_id, aws_secret_access_key)
aws = global_aws_config(; region=aws_region, creds=creds)

const CORS_HEADERS = [
    "Access-Control-Allow-Origin" => "*",
    "Access-Control-Allow-Headers" => "*",
    "Access-Control-Allow-Methods" => "POST, GET, OPTIONS"
]

# https://juliaweb.github.io/HTTP.jl/stable/examples/#Cors-Server
function CorsMiddleware(handler)
    return function (req::HTTP.Request)
        # determine if this is a pre-flight request from the browser
        if HTTP.method(req) == "OPTIONS"
            return HTTP.Response(200, CORS_HEADERS)
        else
            return handler(req) # passes the request to the Application
        end
    end
end


const VIRTUALHOST = "/"
#const HOST = "rabbitmq"
const HOST = "127.0.0.1"
const PORT = 8001

# ==============================================================================
# Variabili condivise per lo stato del server e delle simulazioni
# Sarà necessario usare Locks per proteggere l'accesso a queste variabili
# se più thread/tasks le modificano contemporaneamente.
# In questo scenario, le modifiche provengono principalmente dai Tasks delle simulazioni
# e dalle API di Oxygen.
# ==============================================================================
const solver_overall_status = Ref("ready") # ready, busy, error
const active_simulations = Dict{String,Dict{String,Any}}() # ID simulazione -> {status, progress, start_time, etc.}
const simulations_lock = ReentrantLock() # Lock per proteggere `active_simulations`

const commentsEnabled = []

# ==============================================================================
# Funzione per inviare feedback su RabbitMQ (connessione on-demand)
# ==============================================================================
function send_rabbitmq_feedback(data::Dict, routing_key::String, virtualhost::String=VIRTUALHOST, host::String=HOST)
    try
        # 1. Create a connection to RabbitMQ (on-demand)
        connection(; virtualhost=VIRTUALHOST, host=HOST) do conn
            # 2. Create a channel to send messages
            AMQPClient.channel(conn, AMQPClient.UNUSED_CHANNEL, true) do chan
                # 3. Publish the message (make it persistent if it's critical)
                publish_data(data, routing_key, chan)
                println("Feedback RabbitMQ inviato a $(routing_key)")
            end # Channel is closed here
        end # Connection is closed here
    catch e
        println("Errore durante l'invio del feedback RabbitMQ: $(e)")
        # Implementa qui una logica di retry o di logging più sofisticata se necessario
        # (es. scrivere i messaggi non inviati in un log file per ritentarli dopo)
    end
end

# ==============================================================================
# Funzioni per lo stato del solver e la gestione delle simulazioni
# ==============================================================================

# ==============================================================================
# Funzione eseguita sul worker Distributed - usa aws/bucket del worker
# ==============================================================================
function run_solver_on_worker(solver_type_key::String, args...)
    if solver_type_key == "fft"
        doSolvingFFT(args..., aws, aws_bucket_name)
    elseif solver_type_key == "ris"
        doSolvingRis(args..., aws, aws_bucket_name)
    elseif solver_type_key == "aca"
        doSolvingACA(args..., aws, aws_bucket_name)
    elseif solver_type_key == "electric_fields"
        doSolvingElectricFields(args..., aws, aws_bucket_name)
    else
        error("Tipo di simulazione sconosciuto: $(solver_type_key)")
    end
end

# ==============================================================================
# Spawn di una simulazione su un worker Distributed
# ==============================================================================
function force_kill_worker(w_id::Int)
    try
        if w_id != 1
            w = Distributed.worker_from_id(w_id)
            pid = w.config.ospid
            println("Tentativo terminazione forzata Worker $(w_id) con PID: $pid")
            if Sys.iswindows()
                run(`taskkill /F /PID $pid`)
            else
                run(`kill -9 $pid`)
            end
            println("Worker $(w_id) terminato con successo.")
        end
    catch e
        if !(e isa Distributed.ProcessExitedException)
            println("Errore nel clean-up del worker $(w_id) (possibile PID inesistente): ", e)
        end
    end
end

function spawn_worker_simulation(simulation_id::String, simulation_type::String, solver_type_key::String, args...)
    # Crea un worker process dedicato
    worker_id = addprocs(1)[1]
    println("Worker $(worker_id) creato per simulazione $(simulation_id)")

    # Carica il modulo Solver sul worker (remotecall_eval instead of @everywhere, which is toplevel-only)
    Distributed.remotecall_eval(Main, worker_id, :(using Solver))
    Distributed.remotecall_eval(Main, worker_id, :(using Solver.Checkpointing))

    # Traccia la simulazione
    is_already_stopped = lock(simulations_lock) do
        # Salva il worker_id così può essere terminato se arriva un comando nel frattempo
        if haskey(active_simulations, simulation_id)
            active_simulations[simulation_id]["worker_id"] = worker_id
            active_simulations[simulation_id]["status"] = "running"
            active_simulations[simulation_id]["start_time"] = time()
            active_simulations[simulation_id]["progress"] = 0
            return get(active_simulations[simulation_id], "stopped", false)
        end
        return false
    end

    if is_already_stopped
        println("Simulazione $(simulation_id) fermata prima ancora di iniziare. Rimuovo worker $(worker_id).")
        force_kill_worker(worker_id)
        lock(simulations_lock) do
            if haskey(active_simulations, simulation_id)
                active_simulations[simulation_id]["status"] = "stopped"
                active_simulations[simulation_id]["end_time"] = time()
            end
        end
        send_rabbitmq_feedback(Dict("id" => simulation_id, "isStopped" => true), "solver_feedback")
        return
    end

    send_rabbitmq_feedback(Dict("id" => simulation_id, "status" => "running", "type" => simulation_type), "solver_results")

    # Lancia la computazione sul worker
    future = remotecall(run_solver_on_worker, worker_id, solver_type_key, args...)

    # Monitora in background
    Threads.@spawn monitor_worker_simulation(simulation_id, simulation_type, worker_id, future)
end

# ==============================================================================
# Monitoraggio del worker: gestisce completamento, errore, stop
# ==============================================================================
function monitor_worker_simulation(simulation_id::String, simulation_type::String, worker_id::Int, future::Future)
    try
        res = fetch(future)
        if isnothing(res)
            throw(ErrorException("Simulazione interrotta per un errore interno (out of memory o crash/eccezione del solver). Consultare i log per i dettagli."))
        end
        # Completata con successo
        clear_all_checkpoints(simulation_id)
        lock(simulations_lock) do
            if haskey(active_simulations, simulation_id)
                active_simulations[simulation_id]["status"] = "completed"
                active_simulations[simulation_id]["progress"] = 100
                active_simulations[simulation_id]["end_time"] = time()
            end
        end
        send_rabbitmq_feedback(Dict("id" => simulation_id, "status" => "completed", "type" => simulation_type), "solver_results")
    catch e
        was_stopped = lock(simulations_lock) do
            haskey(active_simulations, simulation_id) && get(active_simulations[simulation_id], "stopped", false)
        end
        if was_stopped || e isa ProcessExitedException
            lock(simulations_lock) do
                if haskey(active_simulations, simulation_id)
                    active_simulations[simulation_id]["status"] = "stopped"
                    active_simulations[simulation_id]["end_time"] = time()
                end
            end
            println("Simulazione $(simulation_id) fermata dall'utente.")
            send_rabbitmq_feedback(Dict("id" => simulation_id, "isStopped" => true), "solver_feedback")
        else
            println("Errore nella simulazione $(simulation_id): $(e)")
            lock(simulations_lock) do
                if haskey(active_simulations, simulation_id)
                    active_simulations[simulation_id]["status"] = "failed"
                    active_simulations[simulation_id]["error_message"] = string(e)
                    active_simulations[simulation_id]["end_time"] = time()
                end
            end
            send_rabbitmq_feedback(Dict("id" => simulation_id, "status" => "failed", "type" => simulation_type, "message" => string(e)), "solver_results")
        end
    finally
        # Rimuovi il worker
        try
            force_kill_worker(worker_id)
            println("Worker $(worker_id) terminato.")
        catch
        end
        # Pulizia dopo un delay
        Threads.@spawn begin
            sleep(60)
            lock(simulations_lock) do
                if haskey(active_simulations, simulation_id) && active_simulations[simulation_id]["status"] in ["completed", "failed", "stopped"]
                    delete!(active_simulations, simulation_id)
                end
            end
        end
        lock(simulations_lock) do
            running = any(v -> v["status"] == "running", values(active_simulations))
            if !running
                solver_overall_status[] = "ready"
                send_rabbitmq_feedback(Dict("target" => "solver", "status" => solver_overall_status[]), "server_init")
            end
        end
    end
end

# ==============================================================================
# Funzioni Oxygen.jl (API Web)
# ==============================================================================

function setup_oxygen_routes()
    # Endpoint per lo stato generale del server: NON USATA PER IL MOMENTO
    get("/status") do
        lock(simulations_lock) do
            json(Dict(
                "solver_overall_status" => solver_overall_status[],
                "active_simulations_count" => length(active_simulations),
                "active_simulations" => active_simulations # Potresti voler limitare i dati inviati qui
            ))
        end
    end


    # Endpoint per avviare simulazioni
    @post "/solve" function (req)
        try
            req_data = Oxygen.json(req) # Assume JSON body
            simulation_id = get(req_data, "id", "randomid") # Genera ID se non fornito
            simulation_type = get(req_data, "simulationType", "matrix") # 'matrix', 'ris', 'electric fields'
            mesher = get(req_data, "mesher", "standard")
            lock(simulations_lock) do
                if haskey(active_simulations, simulation_id)
                    #return JSON.json(Dict("error" => "Simulazione con ID $simulation_id già in corso"))
                    HTTP.Response(500, CORS_HEADERS)
                end
                active_simulations[simulation_id] = Dict(
                    "status" => "pending",
                    "progress" => 0,
                    "type" => simulation_type
                )
            end
            # solver_overall_status[] = "busy"
            # send_rabbitmq_feedback(Dict("target" => "solver", "status" => solver_overall_status[]), "server_init")


            # Avvia la simulazione in un task/thread separato
            if simulation_type == "Matrix" && mesher == "standard"
                mesher_file_id = req_data["mesherFileId"]
                mesherOutput = download_json_gz(aws, aws_bucket_name, mesher_file_id)
                # Calculate current checksum for validation
                current_checksum = string(hash(req_data))
                verify_or_clear_checkpoints(simulation_id, current_checksum)
                
                spawn_worker_simulation(
                    simulation_id, "matrix", "fft",
                    mesherOutput,
                    req_data["solverInput"],
                    req_data["solverAlgoParams"],
                    req_data["solverType"],
                    simulation_id
                )
            elseif simulation_type == "Matrix" && mesher == "ris"
                mesher_file_id = req_data["mesherFileId"]
                surface_file_id = req_data["surfaceFileId"]

                mesherOutput = if req_data["mesherType"] === "backend"
                    download_serialized_data(aws, aws_bucket_name, mesher_file_id)
                else
                    m = get_solverInput_from_s3(aws, aws_bucket_name, mesher_file_id, req_data["mesherType"])
                    m["incidence_selection"]["Gamma"] = convertSparseMatrixFromJavascriptToJulia(m["incidence_selection"]["Gamma"])
                    m["incidence_selection"]["A"] = convertSparseMatrixFromJavascriptToJulia(m["incidence_selection"]["A"])
                    m["nodi_coord"] = transpose(hcat(m["nodi_coord"]...))
                    m["volumi"]["coordinate"] = transpose(hcat(m["volumi"]["coordinate"]...))
                    deep_symbolize_keys(m)
                end
                surface = download_json_gz(aws, aws_bucket_name, surface_file_id) # O get_solverInput_from_s3 a seconda del tipo
                surface["sigma"] = Float64.(surface["sigma"])
                surface["S"] = Float64.(surface["S"])
                surface["normale"] = map(inner -> map(Float64, inner), surface["normale"])
                surface[" materials"] = String.(surface["materials"]) # <- typo from older version, kept logic just fixed parsing
                surface["epsr"] = Float64.(surface["epsr"])
                surface["centri"] = map(inner -> map(Float64, inner), surface["centri"])
                
                verify_or_clear_checkpoints(simulation_id, string(hash((
                    mesherOutput[:incidence_selection], mesherOutput[:volumi], surface, mesherOutput[:nodi_coord], mesherOutput[:escalings],
                    req_data["solverInput"], req_data["solverAlgoParams"], req_data["solverType"]
                ))))

                spawn_worker_simulation(
                    simulation_id, "ris", "ris",
                    mesherOutput[:incidence_selection], mesherOutput[:volumi], surface, mesherOutput[:nodi_coord], mesherOutput[:escalings],
                    req_data["solverInput"], req_data["solverAlgoParams"], req_data["solverType"],
                    simulation_id
                )
            elseif simulation_type == "Matrix_ACA" && mesher == "ris"
                mesher_file_id = req_data["mesherFileId"]
                surface_file_id = req_data["surfaceFileId"]

                # Gestione acaSelectedPorts con fallback
                aca_ports = get(req_data, "acaSelectedPorts", [0])
                ports_to_excite = if isa(aca_ports, AbstractArray)
                    [Int(p) + 1 for p in aca_ports]  # Convert to 1-based indexing
                else
                    [1]  # Default fallback
                end
                println(req_data)
                mesherOutput = if req_data["mesherType"] === "backend"
                    download_serialized_data(aws, aws_bucket_name, mesher_file_id)
                else
                    m = get_solverInput_from_s3(aws, aws_bucket_name, mesher_file_id, req_data["mesherType"])
                    m["incidence_selection"]["Gamma"] = convertSparseMatrixFromJavascriptToJulia(m["incidence_selection"]["Gamma"])
                    m["incidence_selection"]["A"] = convertSparseMatrixFromJavascriptToJulia(m["incidence_selection"]["A"])
                    m["nodi_coord"] = transpose(hcat(m["nodi_coord"]...))
                    m["volumi"]["coordinate"] = transpose(hcat(m["volumi"]["coordinate"]...))
                    deep_symbolize_keys(m)
                end
                surface = download_json_gz(aws, aws_bucket_name, surface_file_id) # O get_solverInput_from_s3 a seconda del tipo
                surface["sigma"] = Float64.(surface["sigma"])
                surface["S"] = Float64.(surface["S"])
                surface["normale"] = map(inner -> map(Float64, inner), surface["normale"])
                surface["materials"] = String.(surface["materials"])
                surface["epsr"] = Float64.(surface["epsr"])
                surface["centri"] = map(inner -> map(Float64, inner), surface["centri"])
                
                verify_or_clear_checkpoints(simulation_id, string(hash((
                    mesherOutput[:incidence_selection], mesherOutput[:volumi], surface, mesherOutput[:nodi_coord], mesherOutput[:escalings],
                    req_data["solverInput"], req_data["solverAlgoParams"], req_data["solverType"], ports_to_excite
                ))))

                spawn_worker_simulation(
                    simulation_id, "ris", "aca",
                    mesherOutput[:incidence_selection], mesherOutput[:volumi], surface, mesherOutput[:nodi_coord], mesherOutput[:escalings],
                    req_data["solverInput"], req_data["solverAlgoParams"], req_data["solverType"], ports_to_excite,
                    simulation_id
                )
            elseif simulation_type == "Electric Fields"
                mesher_file_id = req_data["mesherFileId"]
                surface_file_id = req_data["surfaceFileId"]

                mesherOutput = if req_data["mesherType"] === "backend"
                    download_serialized_data(aws, aws_bucket_name, mesher_file_id)
                else
                    m = get_solverInput_from_s3(aws, aws_bucket_name, mesher_file_id, req_data["mesherType"])
                    m["incidence_selection"]["Gamma"] = convertSparseMatrixFromJavascriptToJulia(m["incidence_selection"]["Gamma"])
                    m["incidence_selection"]["A"] = convertSparseMatrixFromJavascriptToJulia(m["incidence_selection"]["A"])
                    m["nodi_coord"] = transpose(hcat(m["nodi_coord"]...))
                    m["volumi"]["coordinate"] = transpose(hcat(m["volumi"]["coordinate"]...))
                    deep_symbolize_keys(m)
                end
                surface = download_json_gz(aws, aws_bucket_name, surface_file_id) # O get_solverInput_from_s3
                surface["sigma"] = Float64.(surface["sigma"])
                surface["S"] = Float64.(surface["S"])
                surface["normale"] = map(inner -> map(Float64, inner), surface["normale"])
                surface["materials"] = String.(surface["materials"])
                surface["epsr"] = Float64.(surface["epsr"])
                surface["centri"] = map(inner -> map(Float64, inner), surface["centri"])
                
                verify_or_clear_checkpoints(simulation_id, string(hash((
                    mesherOutput[:incidence_selection], mesherOutput[:volumi], surface, mesherOutput[:nodi_coord], mesherOutput[:escalings],
                    req_data["solverInput"], req_data["solverAlgoParams"], req_data["solverType"],
                    req_data["theta"], req_data["phi"], req_data["e_theta"], req_data["e_phi"],
                    req_data["baricentro"], req_data["r_circ"], req_data["times"],
                    req_data["signal_type_E"], req_data["ind_freq_interest"]
                ))))
                
                spawn_worker_simulation(
                    simulation_id, "electric fields", "electric_fields",
                    mesherOutput[:incidence_selection], mesherOutput[:volumi], surface, mesherOutput[:nodi_coord], mesherOutput[:escalings],
                    req_data["solverInput"], req_data["solverAlgoParams"], req_data["solverType"],
                    req_data["theta"], req_data["phi"], req_data["e_theta"], req_data["e_phi"],
                    req_data["baricentro"], req_data["r_circ"], req_data["times"],
                    req_data["signal_type_E"], req_data["ind_freq_interest"],
                    simulation_id
                )
            else
                #return JSON.json(Dict("error" => "Unsupported simulation type: $(simulation_type)"))
                return HTTP.Response(500, CORS_HEADERS)
            end
            return HTTP.Response(200, CORS_HEADERS)
            #JSON.json(Dict("message" => "Simulation started", "id" => simulation_id, "status" => "accepted"))
        catch e
            println("Errore nell'avvio della simulazione: $(e)")
            #JSON.json(Dict("error" => "Failed to start simulation: $(e)"))
            return HTTP.Response(500, CORS_HEADERS)
        end
    end

    @post "/get_results_electric_fields" function (req)
        file_id = queryparams(req)["file_id"]
        freq_index = Oxygen.json(req)["freq_index"]
        id = Oxygen.json(req)["id"]
        try
            # Scarica il file grezzo o il JSON gz compresso da S3
            res = download_json_gz(aws, aws_bucket_name, file_id)
            resultsToPublish = Dict(
                "Vp" => res["Vp"],
                "Ex" => JSON.json(JSON.parse(res["Ex"])[freq_index]),
                "Ey" => JSON.json(JSON.parse(res["Ey"])[freq_index]),
                "Ez" => JSON.json(JSON.parse(res["Ez"])[freq_index]),
                "Ex_3D" => JSON.json(JSON.parse(res["Ex_3D"])[freq_index]),
                "Ey_3D" => JSON.json(JSON.parse(res["Ey_3D"])[freq_index]),
                "Ez_3D" => JSON.json(JSON.parse(res["Ez_3D"])[freq_index]),
                "Hx_3D" => JSON.json(JSON.parse(res["Hx_3D"])[freq_index]),
                "Hy_3D" => JSON.json(JSON.parse(res["Hy_3D"])[freq_index]),
                "Hz_3D" => JSON.json(JSON.parse(res["Hz_3D"])[freq_index]),
                "centri_oss_3D" => res["centri_oss_3D"],
                "distanze_3D" => res["distanze_3D"],
                "theta_vals" => res["theta_vals"],
                "x_grid" => res["x_grid"],
                "y_grid" => res["y_grid"],
                "z_grid" => res["z_grid"],
                "baricentro" => res["baricentro"],
                "f" => res["f"]
            )
            dataToReturn = Dict(
                "results" => resultsToPublish,
                "simulationType" => "electric fields",
                "id" => id
            )
            send_rabbitmq_feedback(dataToReturn, "solver_results")
            #JSON.json("dati restituiti")
            return HTTP.Response(200, CORS_HEADERS)
        catch e
            println("Errore nel recupero dei risultati per $(file_id): $(e)")
            #JSON.json(Dict("error" => "Could not retrieve results for $(file_id): $(e)"))
            return HTTP.Response(500, CORS_HEADERS)
        end
    end

    @post "/get_results_matrix" function (req)
        file_id = queryparams(req)["file_id"]
        port_index = Oxygen.json(req)["port_index"]
        try
            res = download_json_gz(aws, aws_bucket_name, file_id)
            matrixZ = JSON.parse(res["matrices"]["matrix_Z"])
            matrixS = JSON.parse(res["matrices"]["matrix_S"])
            matrixY = JSON.parse(res["matrices"]["matrix_Y"])
            dataToReturn = Dict(
                "portIndex" => port_index,
                "results" => Dict(
                    "matrixZ" => matrixZ[port_index+1],
                    "matrixS" => matrixS[port_index+1],
                    "matrixY" => matrixY[port_index+1],
                ),
                "simulationType" => "matrix"
            )
            #println(dataToReturn)
            send_rabbitmq_feedback(dataToReturn, "solver_results")
            #JSON.json("dati restituiti")
            return HTTP.Response(200, CORS_HEADERS)
        catch e
            println("Errore nel recupero dei risultati per $(file_id): $(e)")
            #JSON.json(Dict("error" => "Could not retrieve results for $(file_id): $(e)"))
            return HTTP.Response(500, CORS_HEADERS)
        end
    end

    # Endpoint per fermare una simulazione - termina il worker process
    @post "/stop_computation" function (req)
        sim_id = queryparams(req)["sim_id"]

        worker_id = lock(simulations_lock) do
            if haskey(active_simulations, sim_id)
                status = active_simulations[sim_id]["status"]
                if status in ["running", "pending"]
                    active_simulations[sim_id]["stopped"] = true
                    val = get(active_simulations[sim_id], "worker_id", nothing)
                    println("Richiesta stop per simulazione $(sim_id) (stato: $(status), worker: $(val))")
                    return val
                end
            end
            return nothing
        end

        if !isnothing(worker_id)
            try
                println("Terminazione worker $(worker_id) per simulazione $(sim_id)...")
                force_kill_worker(worker_id)
                println("Worker $(worker_id) terminato.")
            catch e
                println("Errore nella terminazione del worker $(worker_id): $(e)")
            end
            return HTTP.Response(200, CORS_HEADERS)
        else
            # Se worker_id è nothing ma abbiamo settato stopped = true (caso pending), 
            # lo spawn_worker_simulation se ne accorgerà e pulirà tutto.
            println("Simulazione $(sim_id) non ancora avviata su un worker o già completata. Flag stop impostato.")
            return HTTP.Response(200, CORS_HEADERS)
        end
    end
end

# ==============================================================================
# Main execution flow
# ==============================================================================

function force_kill_workers()
    println("Terminazione forzata di tutti i worker attivi...")
    try
        for w_id in Distributed.workers()
            if w_id != 1
                try
                    w = Distributed.worker_from_id(w_id)
                    pid = w.config.ospid
                    if Sys.iswindows()
                        run(`taskkill /F /PID $pid`)
                    else
                        run(`kill -9 $pid`)
                    end
                catch e
                    if !(e isa Distributed.ProcessExitedException)
                        println("Errore terminazione forzata worker $(w_id): ", e)
                    end
                end
            end
        end
    catch e
        if !(e isa Distributed.ProcessExitedException)
            println("Errore nel clean-up workers: ", e)
        end
    end
end

function julia_main()
    # Invia lo stato iniziale del solver tramite RabbitMQ
    is_building_app = get(ENV, "JULIA_APP_BUILD", "false") == "true"
    if !is_building_app
        send_rabbitmq_feedback(Dict("target" => "solver", "status" => "starting"), "server_init")
        solver_overall_status[] = "starting"
        println("Configurazione delle rotte Oxygen...")
        setup_oxygen_routes()

        println("Avvio del server Oxygen...")
    end

    if !is_building_app
        try
            #up(8001, async = true) #con async a true non blocca il thread principale
            #serve(middleware=[CorsMiddleware], host="0.0.0.0", port=8001, async=true)
            serve(middleware=[CorsMiddleware], host="127.0.0.1", port=8001, async=true)
            # Invia lo stato "ready" dopo aver avviato Oxygen e precompilato
            send_rabbitmq_feedback(Dict("target" => "solver", "status" => "ready"), "server_init")
            solver_overall_status[] = "ready"
            while true
                sleep(1)
            end
        catch ex
            if ex isa InterruptException
                println("Server Oxygen interrotto da Ctrl-C.")
            else
                println("Eccezione durante l'esecuzione del server Oxygen: $(ex)")
            end
        finally
            println("Server Oxygen sta per spegnersi. Invio stato 'idle' a RabbitMQ.")
            send_rabbitmq_feedback(Dict("target" => "solver", "status" => "idle"), "server_init")
            solver_overall_status[] = "idle"
            exit() # Chiude il processo Julia
        end
    else
        # Se siamo nel processo di PackageCompiler.jl, facciamo solo le configurazioni
        # e poi la funzione `main` terminerà naturalmente.
        println("Processo di PackageCompiler.jl in corso (generazione output). Il server non verrà avviato.")
        # Non è necessario aggiungere qui chiamate a `Pkg.precompile()` perché
        # `PackageCompiler.jl` lo gestisce già.
        # Evita chiamate a funzioni che hanno effetti collaterali esterni (es. network I/O)
        # o che richiedono un ambiente di runtime completo.
    end
end

# Punto di ingresso principale del tuo server
# Avvia il server solo sul processo principale (myid()==1).
# Quando un worker Distributed carica il modulo, NON deve avviare il server.
if myid() == 1
    Base.exit_on_sigint(false)
    try
        julia_main()
    catch ex
        if ex isa InterruptException
            println("Catturato Ctrl-C nel blocco principale. Chiusura pulita.")
            println("Killo tutti i worker ancora attivi...")
            force_kill_workers()
            if get(ENV, "JULIA_APP_BUILD", "false") != "true"
                send_rabbitmq_feedback(Dict("target" => "solver", "status" => "idle"), "server_init")
            end
            exit()
        else
            println("Eccezione non gestita nel server principale: $(ex)")
            println("Killo tutti i worker ancora attivi...")
            force_kill_workers()
            if get(ENV, "JULIA_APP_BUILD", "false") != "true"
                send_rabbitmq_feedback(Dict("target" => "solver", "status" => "error", "message" => string(ex)), "server_init")
            end
            exit()
        end
    end
end


# DotEnv.load!(joinpath(@__DIR__, "..", ".env"))

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