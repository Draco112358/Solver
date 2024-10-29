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
          force_compile2()
          # EXCG_DIRECT = "MyDirectExcg"
          # @assert exchange_declare(chan1, EXCG_DIRECT, EXCHANGE_TYPE_DIRECT)
          println(" [*] Waiting for messages. To exit press CTRL+C")
          # 3. Declare a queue
          management_queue = "management_solver"
          #queue_bind(chan, "mesher_results", EXCG_DIRECT, "mesher_results")

          # 4. Setup function to receive message
          on_receive_management = (msg) -> begin
              basic_ack(chan, msg.delivery_tag)
              data = JSON.parse(String(msg.data))
              #data = String(msg.data)
              println(data["message"])
              if (data["message"] == "solving")
                # mesherOutput = JSON.parsefile(data["body"]["mesherFileId"])
                mesherOutput = download_json_gz(aws, aws_bucket_name, data["body"]["mesherFileId"])
                Threads.@spawn doSolving(mesherOutput, data["body"]["solverInput"], data["body"]["solverAlgoParams"], data["body"]["solverType"], data["body"]["id"], aws, aws_bucket_name; chan)
              elseif data["message"] == "stop"
                stop_condition[] = 1.0
              elseif data["message"] == "stop_computation"
                push!(stopComputation, data["id"])
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