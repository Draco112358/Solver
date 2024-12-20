using AMQPClient, JSON

const VIRTUALHOST = "/"
const HOST = "35.224.53.85"

function receive()
  # 1. Create a connection to the localhost or 127.0.0.1 of virtualhost '/'
  connection(; virtualhost=VIRTUALHOST, host=HOST) do conn
      # 2. Create a channel to send messages
      AMQPClient.channel(conn, AMQPClient.UNUSED_CHANNEL, true) do chan
        # 4. Setup function to receive message
        msg1 = basic_get(chan, "solver_results", false)
        msg2 = basic_get(chan, "solver_results", false)
        msg3 = basic_get(chan, "solver_results", false)

        # check if we got a message
        if msg1 !== nothing
            # process msg...
            json_msg1 = JSON.parse(String(msg1.data))
            println(json_msg1["partial"])
            # if (json_msg["partial"] == false)
            #     println("msg partial")
            # end
            #println(JSON.parse(String(msg.data)))
            # acknowledge receipt
            #basic_ack(chan1, msg.delivery_tag)
        end
        if msg2 !== nothing
            # process msg...
            json_msg2 = JSON.parse(String(msg2.data))
            println(json_msg2["partial"])
            # if (json_msg["partial"] == false)
            #     println("msg partial")
            # end
            #println(JSON.parse(String(msg.data)))
            # acknowledge receipt
            #basic_ack(chan1, msg.delivery_tag)
        end
        if msg3 !== nothing
            # process msg...
            json_msg3 = JSON.parse(String(msg3.data))
            println(json_msg3["partial"])
            open("message.json", "w") do f
                JSON.print(f, json_msg3)
            end
            # if (json_msg["partial"] == false)
            #     println("msg partial")
            # end
            #println(JSON.parse(String(msg.data)))
            # acknowledge receipt
            #basic_ack(chan1, msg.delivery_tag)
        end
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