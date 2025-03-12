using AMQPClient, JSON, GZip, AWSS3, AWS, CodecZlib, Serialization

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

# function download_json_gz(aws_config, bucket, key)
#     response = s3_get(aws_config, bucket, key)
#     content = transcode(GzipDecompressor, response)
#     GZip.open(key*".tmp.gz", "w") do f
#       write(f, content)
#     end
#     s = IOBuffer()
#     file = gzopen(key*".tmp.gz")
#     while !eof(file)
#       write(s, readline(file))
#     end
#     close(file)
#     Base.Filesystem.rm(key*".tmp.gz", force=true)
#     data2 = String(take!(s))
#     return JSON.parse(data2)
#   end

function download_json_gz(aws_config, bucket, key)
    # Get the gzipped response from S3
    response = s3_get(aws_config, bucket, key)
    
    # Write the response directly to a temporary file.
    tmp_filename = key * ".tmp.gz"
    open(tmp_filename, "w") do f
      write(f, response)
    end

    # Open the temporary file for reading with gzip decompression.
    s = IOBuffer()
    file = gzopen(tmp_filename)
    while !eof(file)
      println("Reading chunk")  
      write(s, readline(file))
    end
    close(file)
    rm(tmp_filename, force=true)
    println("Reading chunk completed")
    data2 = String(take!(s))
    return JSON.parse(data2)
end

function download_serialized_data(aws_config, bucket, key)
    # Get the gzipped response from S3
    response = s3_get(aws_config, bucket, key)
    # Wrap the byte data in an IOBuffer
    io = IOBuffer(response)
    # Deserialize the variable
    variable = Dict(Serialization.deserialize(io))
    return variable
end


function get_solverInput_from_s3(aws, aws_bucket_name::String, fileName::String, mesherType)
    resposnse_dict = Dict()
    if (s3_exists(aws, aws_bucket_name, fileName))
        response = s3_get(aws, aws_bucket_name, fileName)
        if mesherType == "backend"
            resposnse_dict = to_standard_dict(response)
        else
            resposnse_dict = JSON.parse(String(response))
        end
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

function convertSparseMatrixFromJavascriptToJulia(sparseMatrixJavascript)
    values_js = Float64.(sparseMatrixJavascript["_values"])
    index_js  = sparseMatrixJavascript["_index"]
    ptr_js    = sparseMatrixJavascript["_ptr"]
    m, n    = sparseMatrixJavascript["_size"]
    # Julia uses 1-based indexing so adjust the indices:
    row_julia    = [i + 1 for i in index_js]
    col_ptr_julia = [p + 1 for p in ptr_js]
    sparseMatrixJulia = SparseMatrixCSC(m, n, col_ptr_julia, row_julia, values_js)
    return sparseMatrixJulia
end

function deep_symbolize_keys(x)
    if x isa AbstractDict
        # Create a new Dict with Symbol keys and recursively converted values.
        return Dict(Symbol(k) => deep_symbolize_keys(v) for (k, v) in x)
    elseif x isa AbstractVector
        # Recursively process arrays in case they contain dictionaries.
        return [deep_symbolize_keys(item) for item in x]
    else
        # Return any other value unchanged.
        return x
    end
end

function saveOnS3GZippedResults(fileName::String, data::Dict, aws_config, bucket_name)
    res_id = fileName*"_results.json.gz"
    if(s3_exists(aws_config, bucket_name, res_id))
      s3_delete(aws_config, bucket_name, res_id)
    end
    upload_json_gz(aws_config, bucket_name, res_id, data)
    return res_id
  end
  
  function upload_json_gz(aws_config, bucket_name, file_name, data_to_save)
    dato_compresso = transcode(GzipCompressor, JSON.json(data_to_save))
    s3_put(aws_config, bucket_name, file_name, dato_compresso)
  end