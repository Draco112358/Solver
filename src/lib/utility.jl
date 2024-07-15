using AMQPClient, JSON, GZip

function publish_data(result::Dict, queue::String, chan)
    data = convert(Vector{UInt8}, codeunits(JSON.json(result)))
    message = Message(data, content_type="application/json", delivery_mode=PERSISTENT)
    basic_publish(chan, message; exchange="", routing_key=queue)
end

function is_stopped_computation(id::String, chan)
    if length(filter(i -> i == id, stopComputation)) > 0
        filter!(i -> i != id, stopComputation)
        publish_data(Dict("id" => id, "isStopped" => true), "solver_results", chan)
        return true
    end
    return false
end

function read_gzipped_json_file(file_path::String) :: Dict
    mesh_file = gzopen(file_path)
    s = IOBuffer()
    while !eof(mesh_file)
        write(s, readline(mesh_file))
    end
    close(mesh_file)
    data2 = String(take!(s))
    return JSON.parse(data2)
end