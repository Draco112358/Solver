function ComputeMatrixVectorRis(x, w, incidence_selection, P_data, Lp_data, Z_self, Yle, invZ, invP, F, resProd)
    m = size(incidence_selection[:A], 1)
    ns = size(incidence_selection[:Gamma], 2)

    # Extract vectors from x
    I = @view x[1:m]
    Q = @view x[m+1:m+ns]
    Phi = @view x[m+ns+1:end]

    mx = incidence_selection[:mx]
    my = incidence_selection[:my]
    mz = incidence_selection[:mz]

    # Initialize Y1
    Y1 = zeros(ComplexF64, m)
    @views Y1[1:mx] .= Lp_data[:Lp_x] * I[1:mx]
    @views Y1[mx+1:mx+my] .= Lp_data[:Lp_y] * I[mx+1:mx+my]
    @views Y1[mx+my+1:mx+my+mz] .= Lp_data[:Lp_z] * I[mx+my+1:mx+my+mz]

    # Update Y1 with additional terms
    Y1 .= 1im * w * Y1 .+ Z_self .* I .+ incidence_selection[:A] * Phi

    # Compute Y2
    Y2 = P_data[:P] * Q .- transpose(incidence_selection[:Gamma]) * Phi

    # Compute Y3
    Y3 = -transpose(incidence_selection[:A]) * I .+ Yle * Phi .+ 1im * w * (incidence_selection[:Gamma] * Q)

    # Combine results
    MatrixVector = precond_3_3_vector(F, invZ, invP, incidence_selection[:A], incidence_selection[:Gamma], w, Y1, Y2, Y3)

    return MatrixVector
end

function precond_3_3_vector(lu, invZ, invP, A, Gamma, w, X1, X2, X3)
    n1 = length(X1)
    n2 = length(X2)
    n3 = length(X3)

    i1 = 1:n1
    i2 = n1 + 1:n1 + n2
    i3 = n1 + n2 + 1:n1 + n2 + n3

    Y = zeros(ComplexF64, n1 + n2 + n3)

    M1 = prod_real_complex_Ris(invZ, X1)
    M2 = lu \ prod_real_complex_Ris(transpose(A), M1)
    M3 = prod_real_complex_Ris(invP, X2)
    M4 = lu \ prod_real_complex_Ris(Gamma, M3)
    M5 = lu \ X3

    @views Y[i1] .= M1 .- prod_real_complex_Ris(invZ, prod_real_complex_Ris(A, M2)) .+ 1im * w * prod_real_complex_Ris(invZ, prod_real_complex_Ris(A, M4)) .- prod_real_complex_Ris(invZ, prod_real_complex_Ris(A, M5))
    @views Y[i2] .= prod_real_complex_Ris(invP, prod_real_complex_Ris(transpose(Gamma), M2)) .+ M3 .- 1im * w * prod_real_complex_Ris(invP, prod_real_complex_Ris(transpose(Gamma), M4)) .+ prod_real_complex_Ris(invP, prod_real_complex_Ris(transpose(Gamma), M5))
    @views Y[i3] .= M2 .- 1im * w * M4 .+ M5

    return Y
end

function prod_real_complex_Ris(A, x)
    # A is a N x N real matrix and x is a complex vector
    N = size(A, 1)
    y = zeros(ComplexF64, N)
    y .= A * real(x) .+ 1im * (A * imag(x))
    return y
end