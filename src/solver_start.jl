DotEnv.load!()

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
    return function(req::HTTP.Request)
        # determine if this is a pre-flight request from the browser
        if HTTP.method(req)=="OPTIONS"
            return HTTP.Response(200, CORS_HEADERS)  
        else 
            return handler(req) # passes the request to the Application
        end
    end
end


const VIRTUALHOST = "/"
const HOST = "127.0.0.1"

# ==============================================================================
# Variabili condivise per lo stato del server e delle simulazioni
# Sarà necessario usare Locks per proteggere l'accesso a queste variabili
# se più thread/tasks le modificano contemporaneamente.
# In questo scenario, le modifiche provengono principalmente dai Tasks delle simulazioni
# e dalle API di Oxygen.
# ==============================================================================
const solver_overall_status = Ref("ready") # ready, busy, error
const active_simulations = Dict{String, Dict{String, Any}}() # ID simulazione -> {status, progress, start_time, etc.}
const simulations_lock = ReentrantLock() # Lock per proteggere `active_simulations`
# const stopComputation = []
const stopComputation = Dict{String, Ref{Bool}}() # ID simulazione -> Ref{Bool} per il flag di stop
const stop_computation_lock = ReentrantLock() # Aggiungi un lock per proteggere stopComputation
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
                # Dichiara la coda dei risultati come durevole, se non lo è già
                # Questo è importante se il client RabbitMQ si aspetta messaggi durevoli
                # AMQPClient.queue_declare(chan, "solver_results"; durable=true)
                # AMQPClient.queue_declare(chan, "server_init"; durable=true) # Per stati generali

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

# Funzione wrapper per le tue funzioni di solving originali
# Questa funzione gestirà il ciclo di vita di una simulazione
# e invierà feedback su RabbitMQ.
function run_simulation_task(
    simulation_id::String,
    solver_function::Function, # es. doSolvingFFT, doSolvingRis, doSolvingElectricFields
    args...; # Argomenti specifici per la funzione solver
    simulation_type::String
)
    lock(simulations_lock) do
        active_simulations[simulation_id] = Dict(
            "status" => "running",
            "progress" => 0,
            "start_time" => time(),
            "type" => simulation_type
        )
    end
    send_rabbitmq_feedback(Dict("id" => simulation_id, "status" => "running", "type" => simulation_type), "solver_results")

    try
        # Precompila, se non lo hai già fatto in fase di avvio del server
        # force_compile2() # Potresti volerlo fare una volta all'avvio del server Julia

        # Esegui la simulazione
        # La funzione solver_function DOVRA' essere modificata per accettare
        # un callback o un canale per il progresso, e per i feedback intermedi.
        # Per ora, si assume che pubblichi solo il risultato finale.
        results = solver_function(args...)

        # Simulazione completata
        lock(simulations_lock) do
            active_simulations[simulation_id]["status"] = "completed"
            active_simulations[simulation_id]["progress"] = 100
            active_simulations[simulation_id]["end_time"] = time()
        end
        send_rabbitmq_feedback(Dict("id" => simulation_id, "status" => "completed", "type" => simulation_type, "results" => results), "solver_results")

    catch e
        println("Errore critico nella simulazione $(simulation_id): $(e)")
        lock(simulations_lock) do
            active_simulations[simulation_id]["status"] = "failed"
            active_simulations[simulation_id]["error_message"] = string(e)
            active_simulations[simulation_id]["end_time"] = time()
        end
        send_rabbitmq_feedback(Dict("id" => simulation_id, "status" => "failed", "type" => simulation_type, "message" => string(e)), "solver_results")
    finally
        # Rimuovi la simulazione dalla lista delle attive dopo un po'
        # o sposta in una lista di "simulate terminate"
        Threads.@spawn begin
            sleep(60) # Mantieni i risultati per 1 minut0
            lock(simulations_lock) do
                if haskey(active_simulations, simulation_id) && active_simulations[simulation_id]["status"] in ["completed", "failed"]
                    delete!(active_simulations, simulation_id)
                    println("Simulazione $(simulation_id) rimossa dalla lista attiva.")
                end
            end
        end
        # Aggiorna lo stato generale del solver se non ci sono altre simulazioni attive
        lock(simulations_lock) do
            if isempty(active_simulations)
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
    @post "/solve" function(req)
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
                Threads.@spawn run_simulation_task(
                    simulation_id,
                    doSolvingFFT,
                    mesherOutput,
                    req_data["solverInput"],
                    req_data["solverAlgoParams"],
                    req_data["solverType"],
                    simulation_id, # ID della simulazione passato al solver
                    aws, aws_bucket_name;
                    simulation_type="matrix"
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
                surface["materials"] = String.(surface["materials"])
                surface["epsr"] = Float64.(surface["epsr"])
                surface["centri"] = map(inner -> map(Float64, inner), surface["centri"])
                Threads.@spawn run_simulation_task(
                    simulation_id,
                    doSolvingRis,
                    mesherOutput[:incidence_selection], mesherOutput[:volumi], surface, mesherOutput[:nodi_coord], mesherOutput[:escalings],
                    req_data["solverInput"], req_data["solverAlgoParams"], req_data["solverType"],
                    simulation_id, # ID della simulazione
                    aws, aws_bucket_name;
                    simulation_type="ris"
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
                Threads.@spawn run_simulation_task(
                    simulation_id,
                    doSolvingElectricFields,
                    mesherOutput[:incidence_selection], mesherOutput[:volumi], surface, mesherOutput[:nodi_coord], mesherOutput[:escalings],
                    req_data["solverInput"], req_data["solverAlgoParams"], req_data["solverType"],
                    req_data["theta"], req_data["phi"], req_data["e_theta"], req_data["e_phi"],
                    req_data["baricentro"], req_data["r_circ"], req_data["times"],
                    req_data["signal_type_E"], req_data["ind_freq_interest"],
                    simulation_id, # ID della simulazione
                    aws, aws_bucket_name;
                    simulation_type="electric fields"
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

    @post "/get_results_electric_fields" function(req)
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

    @post "/get_results_matrix" function(req)
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

    # Endpoint per fermare una simulazione (solo se supportato dal solver)
    @post "/stop_computation" function(req)
        sim_id = queryparams(req)["sim_id"]
    
        lock(stop_computation_lock) do
            if haskey(active_simulations, sim_id)
                if !haskey(stopComputation, sim_id) # Crea il Ref{Bool} se non esiste
                    stopComputation[sim_id] = Ref(false)
                end
                stopComputation[sim_id][] = true # Imposta il flag di stop su true
                println("Richiesta di stop per simulazione $(sim_id) ricevuta. Flag impostato su $(stopComputation[sim_id][]).")
                
                # Opzionale: Invia un feedback immediato al client via RabbitMQ che la richiesta è stata accettata
                #send_rabbitmq_feedback(Dict("id" => sim_id, "status" => "stopping", "type" => active_simulations[sim_id]["type"]), "solver_results")

                #return JSON.json(Dict("message" => "Stop request for simulation $(sim_id) acknowledged.", "status" => "stopping"))
                return HTTP.Response(200, CORS_HEADERS)
            else
                println("Richiesta di stop per simulazione $(sim_id) ma la simulazione non è attiva.")
                #return JSON.json(Dict("error" => "Simulation $sim_id not found or already completed/stopped."))
                return HTTP.Response(200, CORS_HEADERS)
            end
        end
    end
end

function is_stop_requested(sim_id::String)
    lock(stop_computation_lock) do
        return haskey(stopComputation, sim_id) && stopComputation[sim_id][]
    end
end

# ==============================================================================
# Main execution flow
# ==============================================================================

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
            serve(middleware=[CorsMiddleware],port=8001, async=true)
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
Base.exit_on_sigint(false) # Non uscire su Ctrl-C immediatamente
try
    julia_main()
catch ex
    if ex isa InterruptException
        println("Catturato Ctrl-C nel blocco principale. Chiusura pulita.")
        if get(ENV, "JULIA_APP_BUILD", "false") != "true" # Invia solo se non è il compilatore
            send_rabbitmq_feedback(Dict("target" => "solver", "status" => "idle"), "server_init")
        end
        exit()
    else
        println("Eccezione non gestita nel server principale: $(ex)")
        if get(ENV, "JULIA_APP_BUILD", "false") != "true" # Invia solo se non è il compilatore
            send_rabbitmq_feedback(Dict("target" => "solver", "status" => "error", "message" => string(ex)), "server_init")
        end
        exit()
    end
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