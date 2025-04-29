using Base.Threads, AMQPClient, JSON, AWS, AWSS3, DotEnv
include("./lib/solve.jl")
include("./lib/utility.jl")

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
  # 1. Create a connection to the localhost or 127.0.0.1 of virtualhost '/'
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
              basic_ack(chan, msg.delivery_tag)
              data = JSON.parse(String(msg.data))
              println(data["message"])
              if (data["message"] == "solving")
                mesherOutput = download_json_gz(aws, aws_bucket_name, data["body"]["mesherFileId"])
                Threads.@spawn doSolving(mesherOutput, data["body"]["solverInput"], data["body"]["solverAlgoParams"], data["body"]["solverType"], data["body"]["id"], aws, aws_bucket_name; chan)
              elseif data["message"] == "solving ris"
                if data["body"]["mesherType"] === "backend"
                  mesherOutput = download_serialized_data(aws, aws_bucket_name, data["body"]["mesherFileId"])
                  println(typeof(mesherOutput))
                  surface = download_json_gz(aws, aws_bucket_name, data["body"]["surfaceFileId"])
                  println(typeof(surface))
                  Threads.@spawn doSolvingRis(mesherOutput[:incidence_selection], mesherOutput[:volumi], surface, mesherOutput[:nodi_coord], mesherOutput[:escalings], data["body"]["solverInput"], data["body"]["solverAlgoParams"], data["body"]["solverType"], data["body"]["id"], aws, aws_bucket_name; chan)
                else
                  mesherOutput = get_solverInput_from_s3(aws, aws_bucket_name, data["body"]["mesherFileId"], data["body"]["mesherType"])
                  surface = get_solverInput_from_s3(aws, aws_bucket_name, data["body"]["surfaceFileId"], data["body"]["mesherType"])
                  mesherOutput["incidence_selection"]["Gamma"] = convertSparseMatrixFromJavascriptToJulia(mesherOutput["incidence_selection"]["Gamma"])
                  mesherOutput["incidence_selection"]["A"] = convertSparseMatrixFromJavascriptToJulia(mesherOutput["incidence_selection"]["A"])
                  mesherOutput["nodi_coord"] = transpose(hcat(mesherOutput["nodi_coord"]...))
                  mesherOutput["volumi"]["coordinate"] = transpose(hcat(mesherOutput["volumi"]["coordinate"]...))
                  mesherOutput = deep_symbolize_keys(mesherOutput)
                  Threads.@spawn doSolvingRis(mesherOutput[:incidence_selection], mesherOutput[:volumi], surface, mesherOutput[:nodi_coord], mesherOutput[:escalings], data["body"]["solverInput"], data["body"]["solverAlgoParams"], data["body"]["solverType"], data["body"]["id"], aws, aws_bucket_name; chan)
                end
              elseif data["message"] == "solving electric fields"
                if data["body"]["mesherType"] === "backend"
                  mesherOutput = download_serialized_data(aws, aws_bucket_name, data["body"]["mesherFileId"])
                  surface = download_json_gz(aws, aws_bucket_name, data["body"]["surfaceFileId"])
                  Threads.@spawn doSolvingElectricFields(mesherOutput[:incidence_selection], mesherOutput[:volumi], surface, mesherOutput[:nodi_coord], mesherOutput[:escalings], data["body"]["solverInput"], data["body"]["solverAlgoParams"], data["body"]["solverType"], data["body"]["theta"], data["body"]["phi"], data["body"]["e_theta"], data["body"]["e_phi"], data["body"]["baricentro"], data["body"]["r_circ"], data["body"]["times"], data["body"]["signal_type_E"], data["body"]["ind_freq_interest"], data["body"]["id"], aws, aws_bucket_name; chan)
                else
                  mesherOutput = get_solverInput_from_s3(aws, aws_bucket_name, data["body"]["mesherFileId"], data["body"]["mesherType"])
                  surface = get_solverInput_from_s3(aws, aws_bucket_name, data["body"]["surfaceFileId"], data["body"]["mesherType"])
                  mesherOutput["incidence_selection"]["Gamma"] = convertSparseMatrixFromJavascriptToJulia(mesherOutput["incidence_selection"]["Gamma"])
                  mesherOutput["incidence_selection"]["A"] = convertSparseMatrixFromJavascriptToJulia(mesherOutput["incidence_selection"]["A"])
                  mesherOutput["nodi_coord"] = transpose(hcat(mesherOutput["nodi_coord"]...))
                  mesherOutput["volumi"]["coordinate"] = transpose(hcat(mesherOutput["volumi"]["coordinate"]...))
                  mesherOutput = deep_symbolize_keys(mesherOutput)
                  Threads.@spawn doSolvingElectricFields(mesherOutput[:incidence_selection], mesherOutput[:volumi], surface, mesherOutput[:nodi_coord], mesherOutput[:escalings], data["body"]["solverInput"], data["body"]["solverAlgoParams"], data["body"]["solverType"], data["body"]["theta"], data["body"]["phi"], data["body"]["e_theta"], data["body"]["e_phi"], data["body"]["baricentro"], data["body"]["r_circ"], data["body"]["times"], data["body"]["signal_type_E"], data["body"]["ind_freq_interest"], data["body"]["id"], aws, aws_bucket_name; chan)
                end
              elseif data["message"] == "stop"
                stop_condition[] = 1.0
              elseif data["message"] == "stop_computation"
                push!(stopComputation, data["id"])
              elseif data["message"] == "get results"
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
                #println(dataToReturn)
                publish_data(dataToReturn, "solver_results", chan)
              elseif data["message"] == "get results electric fields"
                freq_index = data["body"]["freq_index"]
                res = download_json_gz(aws, aws_bucket_name, data["body"]["fileId"])
                println(size(JSON.parse(res["Ex"])))
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
      println("Interrupted")
  else
      println("Exception: $ex")
  end
end