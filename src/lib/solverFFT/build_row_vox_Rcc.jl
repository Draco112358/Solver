using MKL

function build_row_vox_Rcc(mat_QS, Rcc, freq)
    # Constants
    mu0 = 4 * π * 1e-7
    eps0 = 8.854187816997944e-12
    beta = 2π * freq * sqrt(eps0 * mu0)

    # Element-wise multiplication of mat_QS with the exponential factor
    mat_out = mat_QS .* exp.(-1im * beta * Rcc)

    return mat_out
end
