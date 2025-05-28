using Base.Threads, AMQPClient, JSON, AWS, AWSS3, DotEnv
using Genie, Genie.Router, Genie.Renderer.Json, Genie.Requests

include("./lib/solve3.jl") # Contiene doSolving, doSolvingRis, doSolvingElectricFields
include("./lib/utility.jl") # Contiene download_json_gz, get_solverInput_from_s3, ecc.

DotEnv.load!()

aws_access_key_id = ENV["AWS_ACCESS_KEY_ID"]
aws_secret_access_key = ENV["AWS_SECRET_ACCESS_KEY"]
aws_region = ENV["AWS_DEFAULT_REGION"]
aws_bucket_name = ENV["AWS_BUCKET_NAME"]
creds = AWSCredentials(aws_access_key_id, aws_secret_access_key)
aws = global_aws_config(; region=aws_region, creds=creds)

Genie.config.run_as_server = true
Genie.config.cors_headers["Access-Control-Allow-Origin"] = "*"
# This has to be this way - you should not include ".../*"
Genie.config.cors_headers["Access-Control-Allow-Headers"] = "Content-Type"
Genie.config.cors_headers["Access-Control-Allow-Methods"] ="GET,POST,PUT,DELETE,OPTIONS" 
Genie.config.cors_allowed_origins = ["*"]

const VIRTUALHOST = "/"
const HOST = "127.0.0.1"

# ==============================================================================
# Variabili condivise per lo stato del server e delle simulazioni
# Sarà necessario usare Locks per proteggere l'accesso a queste variabili
# se più thread/tasks le modificano contemporaneamente.
# In questo scenario, le modifiche provengono principalmente dai Tasks delle simulazioni
# e dalle API di Genie.
# ==============================================================================
const solver_overall_status = Ref("ready") # ready, busy, error
const active_simulations = Dict{String, Dict{String, Any}}() # ID simulazione -> {status, progress, start_time, etc.}
const simulations_lock = ReentrantLock() # Lock per proteggere `active_simulations`
const stopComputation = []
const commentsEnabled = []

# ==============================================================================
# Funzione per inviare feedback su RabbitMQ (connessione on-demand)
# ==============================================================================
function send_rabbitmq_feedback(data::Dict, routing_key::String)
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
    solver_function::Function, # es. doSolving, doSolvingRis, doSolvingElectricFields
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
            sleep(300) # Mantieni i risultati per 5 minuti
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
# Funzioni Genie.jl (API Web)
# ==============================================================================

function setup_genie_routes()
    # Endpoint per lo stato generale del server
    route("/status") do
        lock(simulations_lock) do
            json(Dict(
                "solver_overall_status" => solver_overall_status[],
                "active_simulations_count" => length(active_simulations),
                "active_simulations" => active_simulations # Potresti voler limitare i dati inviati qui
            ))
        end
    end

    # Endpoint per caricare file
    # route("/upload", method = POST) do
    #     try
    #         payload = Genie.Requests.filespayload()
    #         if isempty(payload)
    #             return  JSON.json(Dict("error" => "No file uploaded"), status = 400)
    #         end

    #         file = first(values(payload)) # Prende il primo file caricato
    #         filename = file.name
    #         file_data = file.data # Questo è un Vector{UInt8}

    #         # Genera un ID univoco per il file, es. usando UUIDs
    #         file_id = string(Base.UUIDs.uuid4()) * "_" * filename
    #         s3_key = "uploads/" * file_id # La chiave per S3

    #         # Carica il file su S3
    #         AWSS3.s3_put(aws, aws_bucket_name, s3_key, file_data)
    #         println("File $(filename) caricato su S3 come $(s3_key)")

    #         json(Dict("message" => "File uploaded successfully", "file_id" => s3_key))
    #     catch e
    #         println("Errore nell'upload del file: $(e)")
    #         json(Dict("error" => "Failed to upload file: $(e)"), status = 500)
    #     end
    # end

    # Endpoint per avviare simulazioni
    route("/solve", method="POST" ) do
        try
            req_data = jsonpayload() # Assume JSON body
            simulation_id = get(req_data, "id", "randomid") # Genera ID se non fornito
            simulation_type = get(req_data, "simulationType", "matrix") # 'matrix', 'ris', 'electric fields'
            println(simulation_type)
            lock(simulations_lock) do
                if haskey(active_simulations, simulation_id)
                    return JSON.json(Dict("error" => "Simulazione con ID $simulation_id già in corso"))
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
            if simulation_type == "matrix"
                mesher_file_id = req_data["mesherFileId"]
                mesherOutput = download_json_gz(aws, aws_bucket_name, mesher_file_id)
                Threads.@spawn run_simulation_task(
                    simulation_id,
                    doSolving,
                    mesherOutput,
                    req_data["solverInput"],
                    req_data["solverAlgoParams"],
                    req_data["solverType"],
                    simulation_id, # ID della simulazione passato al solver
                    aws, aws_bucket_name;
                    simulation_type="matrix"
                )
            elseif simulation_type == "ris"
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
                println("qui")
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
                return JSON.json(Dict("error" => "Unsupported simulation type: $(simulation_type)"))
            end

            json(Dict("message" => "Simulation started", "id" => simulation_id, "status" => "accepted"))
        catch e
            println("Errore nell'avvio della simulazione: $(e)")
            JSON.json(Dict("error" => "Failed to start simulation: $(e)"))
        end
    end

    # Endpoint per ottenere lo stato di una specifica simulazione
    route("/simulation_status/:id") do
        sim_id = params(:id)
        lock(simulations_lock) do
            if haskey(active_simulations, sim_id)
                return json(active_simulations[sim_id])
            else
                return JSON.json(Dict("error" => "Simulation $sim_id not found or completed"))
            end
        end
    end

    # Endpoint per ottenere i risultati finali (potrebbero essere molto grandi)
    # Questa API riceverebbe un ID di un file da S3 e lo restituirebbe.
    route("/get_results_electric_fields", method="POST") do
        file_id = params(:file_id)
        freq_index = jsonpayload()["freq_index"]
        id = jsonpayload()["id"]
        try
            # Scarica il file grezzo o il JSON gz compresso da S3
            # Assumendo che download_json_gz restituisca un dict
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
            JSON.json("dati restituiti")
        catch e
            println("Errore nel recupero dei risultati per $(file_id): $(e)")
            JSON.json(Dict("error" => "Could not retrieve results for $(file_id): $(e)"))
        end
    end

    # Endpoint per fermare una simulazione (solo se supportato dal solver)
    route("/stop_computation/:id", method = POST) do
        sim_id = params(:id)
        # Implementa qui la logica per fermare la simulazione.
        # Questo tipicamente imposta un flag che il solver controllerà periodicamente.
        # Ad esempio, potresti aggiungere un `stopComputation` Dict mappando ID a boolean.
        # Se `doSolving` e le sue varianti controllano un flag `stopComputation[sim_id]`,
        # allora lo imposterai qui.
        println("Richiesta di stop per simulazione $(sim_id).")
        # Esempio: aggiungi l'ID alla lista di quelle da stoppare
        push!(stopComputation, sim_id) # Assumendo che stopComputation sia ancora globale
        json(Dict("message" => "Stop request for simulation $(sim_id) acknowledged."))
    end
end

# ==============================================================================
# Main execution flow
# ==============================================================================

function main()
    # Invia lo stato iniziale del solver tramite RabbitMQ
    send_rabbitmq_feedback(Dict("target" => "solver", "status" => "starting"), "server_init")
    solver_overall_status[] = "starting"

    # Precompilazione del solver (se lunga, farla qui prima di servire richieste)
    # force_compile2() # Scommenta se vuoi precompilare all'avvio

    println("Configurazione delle rotte Genie...")
    setup_genie_routes()

    println("Avvio del server Genie...")

    # Invia lo stato "ready" dopo aver avviato Genie e precompilato
    send_rabbitmq_feedback(Dict("target" => "solver", "status" => "ready"), "server_init")
    solver_overall_status[] = "ready"

    try
        up(8001, async = false) # Blocca il thread principale
    catch ex
        if ex isa InterruptException
            println("Server Genie interrotto da Ctrl-C.")
        else
            println("Eccezione durante l'esecuzione del server Genie: $(ex)")
        end
    finally
        println("Server Genie sta per spegnersi. Invio stato 'idle' a RabbitMQ.")
        send_rabbitmq_feedback(Dict("target" => "solver", "status" => "idle"), "server_init")
        solver_overall_status[] = "idle"
        exit() # Chiude il processo Julia
    end
end

# Punto di ingresso principale del tuo server
Base.exit_on_sigint(false) # Non uscire su Ctrl-C immediatamente
try
    main()
catch ex
    if ex isa InterruptException
        println("Catturato Ctrl-C nel blocco principale. Chiusura pulita.")
        send_rabbitmq_feedback(Dict("target" => "solver", "status" => "idle"), "server_init")
        exit()
    else
        println("Eccezione non gestita nel server principale: $(ex)")
        send_rabbitmq_feedback(Dict("target" => "solver", "status" => "error", "message" => string(ex)), "server_init")
        exit()
    end
end