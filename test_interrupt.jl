using Base.Threads

# Mock shared state
const active_tasks = Dict{String,Task}()
const simulation_lock = ReentrantLock()

function heavy_computation(id)
    println("[$id] Starting heavy computation... (NO POLLING)")
    try
        # Simulate heavy work without manual checks
        for i in 1:100
            sleep(0.1) # Yield point
            # No if is_stopped check here!
            if i % 10 == 0
                println("[$id] Step $i/100")
            end
        end
        println("[$id] Finished naturally.")
    catch e
        if e isa InterruptException
            println("[$id] CAUGHT INTERRUPT! Cleaning up...")
            printstyled("[$id] Simulation ABORTED successfully.\n", color=:green)
        else
            rethrow(e)
        end
    end
end

function start_server()
    println("Server started. Press Ctrl+C to send SIGINT, or wait for simulated stop.")

    # Start a simulation
    sim_id = "sim_1"
    t = Threads.@spawn heavy_computation(sim_id)

    lock(simulation_lock) do
        active_tasks[sim_id] = t
    end

    # Simulate an external "Stop" command after 2 seconds
    Threads.@spawn begin
        sleep(2)
        println("\n>>> TRIGGERING ABORT for $sim_id <<<\n")
        lock(simulation_lock) do
            if haskey(active_tasks, sim_id)
                task = active_tasks[sim_id]
                schedule(task, InterruptException(), error=true)
            end
        end
    end

    # Keep main thread alive
    try
        wait(t)
    catch e
        # Task might fail, but we wait for it
    end
    println("Server finished.")
end

start_server()
