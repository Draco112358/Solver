function distfcm(center, data)
    # DISTFCM Distance measure in fuzzy c-mean clustering.
    # OUT = DISTFCM(CENTER, DATA) calculates the Euclidean distance
    # between each row in CENTER and each row in DATA, and returns a
    # distance matrix OUT of size M by N, where M and N are row
    # dimensions of CENTER and DATA, respectively, and OUT(I, J) is
    # the distance between CENTER(I,:) and DATA(J,:).
    
    # Initialize the output matrix
    out = zeros(size(center[1], 2), size(data, 1))
    # If the center has more than one dimension, calculate Euclidean distance
    if size(center[1],1) > 1
        for k in range(1, size(center[1], 2))
            out[k, :] = transpose(sqrt.(sum((data .- ones(size(data, 1)) * transpose(center[k])) .^ 2, dims=2)))
        end
    else  # For 1-D data, calculate absolute distance
        for k in range(1,size(center[1], 2))
            out[k, :] = transpose(abs.(center[k] .- data))
        end
    end
    
    return out
end
