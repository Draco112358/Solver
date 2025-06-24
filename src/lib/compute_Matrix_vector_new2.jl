using MKL, LinearAlgebra

function ComputeMatrixVectorNew2(x, w, incidence_selection, P_data, Lp_data, Z_self, Yle, invZ, invP, F, resProd)
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

function precond_3_3_vector(lu, invZ, invP, A, A_t, Gamma, Gamma_t, w, X1, X2, X3)
    n1 = length(X1)
    n2 = length(X2)
    n3 = length(X3)

    i1 = 1:n1
    i2 = n1 + 1:n1 + n2
    i3 = n1 + n2 + 1:n1 + n2 + n3

    Y = zeros(ComplexF64, n1 + n2 + n3)

    M1 = invZ*X1
    M2 = lu \ (A_t*M1)
    M3 = invP*X2
    M4 = lu \ (Gamma*M3)
    M5 = lu \ X3

    @views Y[i1] .= M1 .- invZ*(A*M2) .+ 1im * w * invZ*(A*M4) .- invZ*(A*M5)
    @views Y[i2] .= invP*(Gamma_t*M2) .+ M3 .- 1im * w * invP*(Gamma_t*M4) .+ invP*(Gamma_t*M5)
    @views Y[i3] .= M2 .- 1im * w * M4 .+ M5

    return Y
end

"""
    ComputeMatrixVectorNew2!(y, x, w, isc, P, Lp, Zs, Yle,
                             invZ, invP, luF, work)
Versione in-place senza allocazioni.  
`y`     : vettore risultato già allocato  
`work`  : named‐tuple con i buffer Y1, Y2, Y3 (anch’essi preallocati)
"""
function ComputeMatrixVectorNew2!(y, work, x, w,
                                  isc, P, Lp, Zs, Yle,
                                  invZ, invP, luF)
    A      = isc[:A]
    Γ      = isc[:Gamma]
    At     = isc[:A_t]      # ⇒ pre‐calcola fuori: At = A'
    Γt     = isc[:Gamma_t]    # idem         : Γt = Γ'
    mx, my, mz = isc[:mx], isc[:my], isc[:mz]
    m  = size(A,1)
    ns = size(Γ,2)
    # ───────── Slice di x (solo view, nessuna allocazione)
    I   = @view x[1:m]
    Q   = @view x[m+1:m+ns]
    Φ   = @view x[m+ns+1:end]
    Y1, Y2, Y3 = work.Y1, work.Y2, work.Y3   # alias locali
    # ───────── Y1 = (Lp*I)  poi   Y1 = jω·Y1 + Zs⊙I + A*Φ
    @views begin
        mul!(Y1[1:mx],       Lp[:Lp_x], I[1:mx])        # Lp_x*Ix
        mul!(Y1[mx+1:mx+my], Lp[:Lp_y], I[mx+1:mx+my])  # Lp_y*Iy
        mul!(Y1[mx+my+1:end],Lp[:Lp_z], I[mx+my+1:end]) # Lp_z*Iz
    end
    @. Y1 = im*w*Y1 + Zs*I                              # Zs è vettore diag
    mul!(Y1, A, Φ, 1, 1)                                # Y1 += A*Φ
    # ───────── Y2 = P*Q  - Γᵀ*Φ
    mul!(Y2, P[:P], Q)              # Y2 = P*Q
    mul!(Y2, Γt, Φ, -1, 1)          # Y2 .-= Γᵀ*Φ
    # ───────── Y3 = -Aᵀ*I + Yle*Φ + jω*(Γ*Q)
    mul!(Y3, At, I, -1, 0)          # Y3 = -Aᵀ*I
    mul!(Y3, Yle, Φ, 1, 1)          # Y3 += Yle*Φ
    mul!(Y3, Γ,  Q, im*w, 1)        # Y3 += jω Γ*Q
    # ───────── Precondizionatore in-place
    y = precond_3_3_vector(luF, invZ, invP, A, At, Γ, Γt, w, Y1, Y2, Y3)
    return y
end