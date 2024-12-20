using AMQPClient, JSON, GZip, AWSS3, AWS, CodecZlib

function publish_data(result::Dict, queue::String, chan)
    data = convert(Vector{UInt8}, codeunits(JSON.json(result)))
    message = Message(data, content_type="application/json", delivery_mode=PERSISTENT)
    basic_publish(chan, message; exchange="", routing_key=queue)
end

function is_stopped_computation(id::String, chan)
    if length(filter(i -> i == id, stopComputation)) > 0
        filter!(i -> i != id, stopComputation)
        #publish_data(Dict("id" => id, "isStopped" => true), "solver_results", chan)
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

function download_json_gz(aws_config, bucket, key)
    response = s3_get(aws_config, bucket, key)
    content = transcode(GzipDecompressor, response)
    GZip.open(key*".tmp.gz", "w") do f
      write(f, content)
    end
    s = IOBuffer()
    file = gzopen(key*".tmp.gz")
    while !eof(file)
      write(s, readline(file))
    end
    close(file)
    Base.Filesystem.rm(key*".tmp.gz", force=true)
    data2 = String(take!(s))
    return JSON.parse(data2)
  end

  function get_solverInput_from_s3(aws, aws_bucket_name::String, fileName::String)
    resposnse_dict = Dict()
    if (s3_exists(aws, aws_bucket_name, fileName))
        response = s3_get(aws, aws_bucket_name, fileName)
        resposnse_dict = to_standard_dict(response)
    end
    return resposnse_dict
end

function to_standard_dict(data)
    if isa(data, OrderedCollections.LittleDict)
        # Convert the LittleDict to Dict, applying the function recursively
        return Dict(k => to_standard_dict(v) for (k, v) in data)
    elseif isa(data, AbstractArray)
        # If it's an array, apply the function to each element
        return map(to_standard_dict, data)
    else
        # For other types, return as is
        return data
    end
end