using Base.Threads, AMQPClient, JSON, AWS, AWSS3, DotEnv, BSON
include("./lib/solve.jl")
include("./lib/utility.jl")
include("./lib/solve2.jl")

DotEnv.load!()

aws_access_key_id = ENV["AWS_ACCESS_KEY_ID"]
aws_secret_access_key = ENV["AWS_SECRET_ACCESS_KEY"]
aws_region = ENV["AWS_DEFAULT_REGION"]
aws_bucket_name = ENV["AWS_BUCKET_NAME"]
creds = AWSCredentials(aws_access_key_id, aws_secret_access_key)
aws = global_aws_config(; region=aws_region, creds=creds)

#server = WebsocketServer()

#Threads.@spawn serve(server; verbose=false)

const stopComputation = []
const commentsEnabled = []

function dict_to_state(data::Dict)
  return ElectricFieldsSolvingState(
    data["id"], data["incidence_selection"], data["volumi"], data["superfici"],
    data["nodi_coord"], data["escalings"], data["solverInput"],
    data["solverAlgoParams"], data["solverType"],
    data["theta"], data["phi"], data["e_theta"], data["e_phi"],
    data["baricentro"], data["r_circ"], data["times"], data["signal_type_E"],
    Int.(data["ind_freq_interest"]), data["currentState"], data["currentStateComputeFields"],
    get(data, "progress", 0.0), get(data, "intermediateResults", Dict())
  )
end

function deserialize_with_symbols_recursive(data::Any)
  if data isa Dict
      restored = Dict()
      for (key, value) in data
          if key isa String && startswith(key, "_sym_")
            key = Symbol(replace(key, "_sym_" => ""))
          end
          restored[key] = deserialize_with_symbols_recursive(value)
      end
      return restored
  else
      return data
  end
end


function load_task_state(task_id::String, aws_config, bucket_name, task_type="electric_fields")
    s3_key = "solver_state/$(task_type)_$(task_id).json.gz"
    fileName = "$(task_type)_$(task_id).json.gz"
    try
        if s3_exists(aws_config, bucket_name, s3_key)
            serialized_data = download_state_json_gz(aws_config, bucket_name, s3_key)
            return dict_to_state(deserialize_with_symbols_recursive(serialized_data))
        else
            return nothing
        end
    catch e
        println("Errore nel caricamento dello stato per $task_id: $(e)")
        return nothing
    end
end

function stopSolver(VIRTUALHOST, HOST)
  #nuova connessione con il broker per avvisare il client che il solver è stato stoppato
  connection(; virtualhost=VIRTUALHOST, host=HOST) do conn
    AMQPClient.channel(conn, AMQPClient.UNUSED_CHANNEL, true) do chan
      publish_data(Dict("target" => "solver", "status" => "idle"), "server_init", chan)
    end
  end
  println("Shutdown initiated. The 'idle' status should have been published.")
  exit()
end

function reRunDoSolvingElectricFields(task_id, chan, stopComputation, aws, aws_bucket_name, VIRTUALHOST, HOST, exception)
  existing_state = load_task_state(task_id, aws, aws_bucket_name, "electric_fields")
  if !isnothing(existing_state) && existing_state.currentState != "finished" && existing_state.currentState != "failed"
    println("[$task_id] Ripresa task 'electric fields' dallo stato: $(existing_state.currentState)")
    Threads.@spawn doSolvingElectricFields2(existing_state, aws, aws_bucket_name; chan, stopComputation)
  else
    println("Relevated Exception on Solving electric fields: $(exception)")
    stopSolver(VIRTUALHOST, HOST)
  end
end

function getInitialState(task_id, mesherOutput, surface, data)
  return ElectricFieldsSolvingState(
    task_id,
    mesherOutput[:incidence_selection],
    mesherOutput[:volumi], 
    surface, 
    mesherOutput[:nodi_coord], 
    mesherOutput[:escalings],
    data["body"]["solverInput"],
    data["body"]["solverAlgoParams"],
    data["body"]["solverType"],
    data["body"]["theta"],
    data["body"]["phi"],
    data["body"]["e_theta"],
    data["body"]["e_phi"],
    data["body"]["baricentro"],
    data["body"]["r_circ"],
    data["body"]["times"],
    data["body"]["signal_type_E"],
    Int.(data["body"]["ind_freq_interest"]),
    "initial",
    "compute hc",
    0.0,
    Dict()
)
end


function force_compile2()
  println("------ Precompiling routes...wait for solver to be ready ---------")
  data = open(JSON.parse, "first_run_data.json")
  doSolving(data["mesherOutput"], data["solverInput"], data["solverAlgoParams"], data["solverType"], "init", aws, aws_bucket_name; commentsEnabled=false)
  println("SOLVER READY")
end

const VIRTUALHOST = "/"
const HOST = "127.0.0.1"
const stop_condition = Ref{Float64}(0.0)

function receive()
  connection(; virtualhost=VIRTUALHOST, host=HOST) do conn
      # 2. Create a channel to send messages
      AMQPClient.channel(conn, AMQPClient.UNUSED_CHANNEL, true) do chan
          publish_data(Dict("target" => "solver", "status" => "starting"), "server_init", chan)
          #force_compile2()
          println(" [*] Waiting for messages. To exit press CTRL+C")
          # 3. Declare a queue
          management_queue = "management_solver"

          # 4. Setup function to receive message
          on_receive_management = (msg) -> begin
              #basic_ack(chan, msg.delivery_tag)
              data = JSON.parse(String(msg.data))
              println(data["message"])
              try
                if data["message"] == "ping"
                  publish_data(Dict("target" => "solver", "status" => "ready"), "server_init", chan)
                  basic_ack(chan, msg.delivery_tag)
                end
              catch e
                println("Rilevated Exception on ping: $(e)")
                stopSolver(VIRTUALHOST, HOST)
              end

              try
                if data["message"] == "solving electric fields"
                  if length(stopComputation) > 0
                    pop!(stopComputation)
                  end
                  task_id = data["body"]["id"]
                  existing_state = load_task_state(task_id, aws, aws_bucket_name, "electric_fields")
                  if data["body"]["mesherType"] === "backend"
                    mesherOutput = download_serialized_data(aws, aws_bucket_name, data["body"]["mesherFileId"])
                    surface = download_json_gz(aws, aws_bucket_name, data["body"]["surfaceFileId"])
                    println("[$task_id] Inizio nuovo task 'electric fields'.")
                    initial_state = getInitialState(task_id, mesherOutput, surface, data)
                    Threads.@spawn doSolvingElectricFields2(initial_state, aws, aws_bucket_name; chan, stopComputation)
                    basic_ack(chan, msg.delivery_tag)
                  else
                    mesherOutput = get_solverInput_from_s3(aws, aws_bucket_name, data["body"]["mesherFileId"], data["body"]["mesherType"])
                    surface = get_solverInput_from_s3(aws, aws_bucket_name, data["body"]["surfaceFileId"], data["body"]["mesherType"])
                    mesherOutput["incidence_selection"]["Gamma"] = convertSparseMatrixFromJavascriptToJulia(mesherOutput["incidence_selection"]["Gamma"])
                    mesherOutput["incidence_selection"]["A"] = convertSparseMatrixFromJavascriptToJulia(mesherOutput["incidence_selection"]["A"])
                    mesherOutput["nodi_coord"] = transpose(hcat(mesherOutput["nodi_coord"]...))
                    mesherOutput["volumi"]["coordinate"] = transpose(hcat(mesherOutput["volumi"]["coordinate"]...))
                    mesherOutput = deep_symbolize_keys(mesherOutput)
                    if !isnothing(existing_state) && existing_state.currentState != "finished" && existing_state.currentState != "failed"
                      println("[$task_id] Ripresa task 'electric fields' dallo stato: $(existing_state.currentState)")
                      Threads.@spawn doSolvingElectricFields2(existing_state, aws, aws_bucket_name; chan, stopComputation)
                    else
                      println("[$task_id] Inizio nuovo task 'electric fields'.")
                      initial_state = getInitialState(task_id, mesherOutput, surface, data)
                      Threads.@spawn doSolvingElectricFields2(initial_state, aws, aws_bucket_name; chan, stopComputation)
                    end
                    basic_ack(chan, msg.delivery_tag)
                  end
                end
              catch e
                task_id = data["body"]["id"]
                reRunDoSolvingElectricFields(task_id, chan, stopComputation, aws, aws_bucket_name, VIRTUALHOST, HOST, e)
                basic_ack(chan, msg.delivery_tag)
              end
              if data["message"] == "stop"
                stop_condition[] = 1.0
                basic_ack(chan, msg.delivery_tag)
              end
              if data["message"] == "stop_computation"
                push!(stopComputation, data["id"])
                basic_ack(chan, msg.delivery_tag)
              end
              if data["message"] == "get results"
                res = download_json_gz(aws, aws_bucket_name, data["body"]["fileId"])
                matrixZ = JSON.parse(res["matrices"]["matrix_Z"])
                matrixS = JSON.parse(res["matrices"]["matrix_S"])
                matrixY = JSON.parse(res["matrices"]["matrix_Y"])
                dataToReturn = Dict(
                  "portIndex" => data["body"]["portIndex"],
                  "results" => Dict(
                    "matrixZ" => matrixZ[data["body"]["portIndex"]+1],
                    "matrixS" => matrixS[data["body"]["portIndex"]+1],
                    "matrixY" => matrixY[data["body"]["portIndex"]+1],
                  ),
                  "simulationType" => "matrix"
                )
                publish_data(dataToReturn, "solver_results", chan)
                basic_ack(chan, msg.delivery_tag)
              end
              if data["message"] == "get results electric fields"
                freq_index = data["body"]["freq_index"]
                res = download_json_gz(aws, aws_bucket_name, data["body"]["fileId"])
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
                  "id" => data["body"]["id"]
                )
                publish_data(dataToReturn, "solver_results", chan)
                basic_ack(chan, msg.delivery_tag)
              end
          end

          # 5. Configure Quality of Service
          basic_qos(chan, 0, 1, false)
          success_management, consumer_tag = basic_consume(chan, management_queue, on_receive_management)

          @assert success_management == true
          publish_data(Dict("target" => "solver", "status" => "ready"), "server_init", chan)
          while stop_condition[] != 1.0
              sleep(1)
          end
          # 5. Close the connection
          publish_data(Dict("target" => "solver", "status" => "idle"), "server_init", chan)
          sleep(3)
      end
  end
end


# Don't exit on Ctrl-C
Base.exit_on_sigint(false)
try
  receive()
catch ex
  if ex isa InterruptException
    stopSolver(VIRTUALHOST, HOST)
  else
      println("Exception: $ex")
  end
end