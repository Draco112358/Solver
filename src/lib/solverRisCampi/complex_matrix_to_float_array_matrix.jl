using MKL
function complex_matrix_to_float_array_matrix(complex_matrix::Matrix{ComplexF64})
    rows, cols = size(complex_matrix)
    float_array_matrix = Array{Array{Float64, 1}, 2}(undef, rows, cols)
    for i in 1:rows
        for j in 1:cols
            float_array_matrix[i, j] = [real(complex_matrix[i, j]), imag(complex_matrix[i, j])]
        end
    end
    return float_array_matrix
end