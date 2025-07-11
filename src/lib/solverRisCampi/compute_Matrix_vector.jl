struct PCWork{T}
    b1::Vector{T}      # n3
    b2::Vector{T}      # n3
    b3::Vector{T}      # n3
    z1::Vector{T}      # n1   (M1)
    ztmp::Vector{T}    # n1   scratch invZ
    ztmp2::Vector{T}
    p1::Vector{T}      # n2   (M3)
    ptmp::Vector{T}    # n2   scratch invP
    ptmp2::Vector{T}
end

PCWork(T::Type, n1::Int, n2::Int, n3::Int) = PCWork(
    Vector{T}(undef, n3),
    Vector{T}(undef, n3),
    Vector{T}(undef, n3),
    Vector{T}(undef, n1),
    Vector{T}(undef, n1),
    Vector{T}(undef, n1),
    Vector{T}(undef, n2),
    Vector{T}(undef, n2),
    Vector{T}(undef, n2)
)

function precond_3_3_vector!(
        Y::Vector{ComplexF64},
        lu::SparseArrays.UMFPACK.UmfpackLU{ComplexF64,Int},
        invZ::SparseMatrixCSC{ComplexF64,Int},
        invP::SparseMatrixCSC{Float64,Int},
        A::SparseMatrixCSC{Float64,Int},
        A_t::Transpose{Float64,SparseMatrixCSC{Float64,Int}},
        Γ::SparseMatrixCSC{Float64,Int},
        Γ_t::Transpose{Float64,SparseMatrixCSC{Float64,Int}},
        w::Float64,
        X1::Vector{ComplexF64},
        X2::Vector{ComplexF64},
        X3::Vector{ComplexF64},
        ws::PCWork{ComplexF64})

    n1 = length(X1);  n2 = length(X2);  n3 = length(X3)
    i1 = 1:n1
    i2 = n1+1 : n1+n2
    i3 = n1+n2+1 : n1+n2+n3

    b1, b2, b3 = ws.b1, ws.b2, ws.b3
    M1, ztmp, ztmp2   = ws.z1, ws.ztmp, ws.ztmp2
    M3, ptmp, ptmp2   = ws.p1, ws.ptmp, ws.ptmp2

    # ---- M1 = invZ * X1  (n1) -------------------------------------------
    mul!(M1, invZ, X1)

    # ---- M2 -------------------------------------------------------------
    mul!(b1, A_t, M1)          # RHS
    ldiv!(b2, lu, b1)          # M2  (b2)

    # ---- M3 = invP * X2  (n2) ------------------------------------------
    mul!(M3, invP, X2)

    # ---- M4 -------------------------------------------------------------
    mul!(b1, Γ, M3)           # RHS
    ldiv!(b3, lu, b1)    # M4  (b3)

    # ---- M5 -------------------------------------------------------------
    ldiv!(b1, lu, X3)          # M5  (b1)

    # ---------- blocco Y[i1] --------------------------------------------
    mul!(ztmp, A, b2);  mul!(ztmp2, invZ, ztmp)
    @views Y[i1] .= M1 .- ztmp2

    mul!(ztmp, A, b3);  mul!(ztmp2, invZ, ztmp)
    @views Y[i1] .+= 1im*w .* ztmp2

    mul!(ztmp, A, b1);  mul!(ztmp2, invZ, ztmp)
    @views Y[i1] .-= ztmp2

    # ---------- blocco Y[i2] --------------------------------------------
    mul!(ptmp, Γ_t, b2);  mul!(ptmp2, invP, ptmp)
    @views Y[i2] .= ptmp2 .+ M3

    mul!(ptmp, Γ_t, b3);  mul!(ptmp2, invP, ptmp)
    @views Y[i2] .-= 1im*w .* ptmp2

    mul!(ptmp, Γ_t, b1);  mul!(ptmp2, invP, ptmp)
    @views Y[i2] .+= ptmp2

    # ---------- blocco Y[i3] --------------------------------------------
    @views Y[i3] .= b2 .- 1im*w .* b3 .+ b1

    return Y
end


function ComputeMatrixVector!(y, work, pc_work, x, w,isc, P, Lp, Zs, Yle,invZ, invP, luF)
    @time begin
        A      = isc[:A]
        Γ      = isc[:Gamma]
        At     = isc[:A_t]      # ⇒ pre‐calcola fuori: At = A'
        Γt     = isc[:Gamma_t]    # idem         : Γt = Γ'
        mx, my, mz = isc[:mx], isc[:my], isc[:mz]
       
        m  = size(A,1) #2 allocs
        ns = size(Γ,2)
       
        
        # ───────── Slice di x (solo view, nessuna allocazione)
        
        I   = @view x[1:m] #2 allocs
        Q   = @view x[m+1: m+ns] #10 allocs
        Φ   = @view x[m+ns+1:end]
    
        
        Y1, Y2, Y3 = work.Y1, work.Y2, work.Y3   # alias locali
        # ───────── Y1 = (Lp*I)  poi   Y1 = jω·Y1 + Zs⊙I + A*Φ
        @views begin
            mul!(Y1[1:mx],       Lp[:Lp_x], I[1:mx])        # Lp_x*Ix
            mul!(Y1[mx+1:mx+my], Lp[:Lp_y], I[mx+1:mx+my])  # Lp_y*Iy
            mul!(Y1[mx+my+1:end],Lp[:Lp_z], I[mx+my+1:end]) # Lp_z*Iz
        end

        @. Y1 = im*w*Y1 + Zs*I                          # Zs è vettore diag
        mul!(Y1, A, Φ, 1, 1)                                # Y1 += A*Φ
        # ───────── Y2 = P*Q  - Γᵀ*Φ
        mul!(Y2, P[:P], Q)              # Y2 = P*Q
        mul!(Y2, Γt, Φ, -1, 1)          # Y2 .-= Γᵀ*Φ
        # ───────── Y3 = -Aᵀ*I + Yle*Φ + jω*(Γ*Q)
        mul!(Y3, At, I, -1, 0)          # Y3 = -Aᵀ*I
        mul!(Y3, Yle, Φ, 1, 1)       # Y3 += Yle*Φ
        mul!(Y3, Γ,  Q, im*w, 1)        # Y3 += jω Γ*Q
        # ───────── Precondizionatore in-place
        # precond_3_3_vector!(ComplexF64.(y), luF, invZ, invP, A, At, Γ, Γt, w, Y1, Y2, Y3,
        #                         M1_pre, M2_pre, M3_pre, M4_pre, M5_pre,
        #                         temp_res_n1, temp_res_n2, temp_res_n3
        #     )
        println("time precond")
        @time precond_3_3_vector!(y, luF, invZ, invP, A, At, Γ, Γt, w, Y1, Y2, Y3, pc_work)
        println("time compute matrix vector")
    end
    
end