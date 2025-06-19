using MAT

function saveComplexMatrix(filename::String, A::AbstractArray; varname::String="A")
    matwrite(filename, Dict(varname => A))
end