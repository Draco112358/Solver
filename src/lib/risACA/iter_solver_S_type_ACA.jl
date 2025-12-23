function iter_solver_S_type_ACA(freq, escalings, incidence_selection, P_data, Lp_data, ports, lumped_elements, GMRES_settings, volumi, superfici, use_Zs_in, QS_Rcc_FW, ACA_thres, id, chan, commentsEnabled)

    # PER LA RUGOSITA'
    # rugosity_lines;
    # warning off; (Julia doesn't need this explicit warning off usually)

    freq = freq .* escalings[:freq]

    # GMRES settings ----------------------------
    Inner_Iter::Int64 = GMRES_settings["Inner_Iter"]
    Outer_Iter::Int64 = GMRES_settings["Outer_Iter"] # Note: gmres_custom uses this differently or expects it in a specific way
    # -------------------------------------------

    mx = incidence_selection[:mx]
    my = incidence_selection[:my]
    mz = incidence_selection[:mz]
    m = mx + my + mz
    n::Int64 = size(incidence_selection[:A], 2)
    ns::Int64 = size(incidence_selection[:Gamma], 2)

    w = 2 .* pi .* freq

    nfreq = length(w)

    is = zeros(ComplexF64, n) # Should be complex or float based on usage. MATLAB inits as zeros(n,1) which is double.

    S = zeros(ComplexF64, size(ports[:port_nodes], 1), size(ports[:port_nodes], 1), nfreq)

    Vrest = zeros(ComplexF64, m + n + ns, size(ports[:port_nodes], 1))

    R_chiusura = 50.0

    invCd = zeros(ComplexF64, m)
    not_switched = true

    # Pre-allocation for GMRES product
    resProd = zeros(ComplexF64, m + n + ns) # Estimate size, verify logic later

    local P_sp, Lp_sp
    if QS_Rcc_FW < 2
        P_sp, Lp_sp = prepare_sparse_mats_prec_ACA(P_data, Lp_data)
    end

    for k = 1:nfreq

        if freq[k] / escalings[:freq] >= 1e8 && not_switched

            not_switched = false

            freq = freq ./ escalings[:freq]
            w = 2 .* pi .* freq

            volumi[:R] = volumi[:R] ./ escalings[:R]
            volumi[:Zs_part] = volumi[:Zs_part] ./ escalings[:R]

            if !isempty(volumi[:indici_dielettrici])
                volumi[:Cd] = volumi[:Cd] ./ escalings[:Cd]
            end

            remove_escalings_comp_mats!(Lp_data[:Lpx], Lp_data[:Lpx_comp], escalings[:Lp])
            remove_escalings_comp_mats!(Lp_data[:Lpy], Lp_data[:Lpy_comp], escalings[:Lp])
            remove_escalings_comp_mats!(Lp_data[:Lpz], Lp_data[:Lpz_comp], escalings[:Lp])
            remove_escalings_comp_mats!(P_data, P_data[:P_comp], escalings[:P])

            if QS_Rcc_FW < 2
                P_sp = P_sp ./ escalings[:P]
                Lp_sp = Lp_sp ./ escalings[:Lp]
            end

            Vrest[1:m, :] = Vrest[1:m, :] ./ escalings[:Is]
            if ns > 0
                Vrest[m+1:m+ns, :] = Vrest[m+1:m+ns, :] ./ escalings[:Cd]
            end

            escalings[:Lp] = 1.0
            escalings[:P] = 1.0
            escalings[:R] = 1.0
            escalings[:Cd] = 1.0
            escalings[:Is] = 1.0
            escalings[:Yle] = 1.0
            escalings[:freq] = 1.0
            escalings[:time] = 1.0

        end

        println("Freq n=$(k) - Freq Tot=$(nfreq)")

        if QS_Rcc_FW == 2
            mu0 = 4 * pi * 1e-7
            eps0 = 8.854187816997944e-12
            beta = 2 * pi * freq[k] / escalings[:freq] * sqrt(eps0 * mu0)

            add_delays_P!(P_data, beta, superfici, ACA_thres, escalings, id)
            add_delays_Lp!(Lp_data, volumi, incidence_selection, beta, ACA_thres, escalings, id)

            P_sp, Lp_sp = prepare_sparse_mats_prec_ACA(P_data, Lp_data)
        end

        Yle = build_Yle_S(lumped_elements, [], ports, escalings, n, w[k] / escalings[:freq], R_chiusura,
            lumped_elements[:type], lumped_elements[:R], lumped_elements[:L], lumped_elements[:C])

        if !isempty(volumi[:indici_dielettrici])
            invCd[Int64.(volumi[:indici_dielettrici])] .= 1.0 ./ (1im * w[k] * volumi[:Cd][Int64.(volumi[:indici_dielettrici])])
        end

        local Z_self
        if use_Zs_in == 1
            Zs = real.(sqrt(1im * w[k] / escalings[:freq]) .* volumi[:Zs_part])
            indR = findall(x -> volumi[:R][x] > Zs[x], 1:length(volumi[:R]))
            indZs = setdiff(1:m, indR)
            Z_self = zeros(ComplexF64, m)
            Z_self[indR] = volumi[:R][indR]
            Z_self[indZs] = Zs[indZs]
            Z_self = Z_self .+ invCd
        else
            Z_self = volumi[:R] .+ invCd
        end

        # --------------------- preconditioner ------------------------

        # MatPrec = [sparse(1:m, 1:m, Z_self, m, m)+1im*w[k]*Lp_sp spzeros(m, ns) incidence_selection[:A];
        #     spzeros(ns, m) P_sp -transpose(incidence_selection[:Gamma]);
        #     -transpose(incidence_selection[:A]) 1im*w[k]*incidence_selection[:Gamma] Yle]

        # Explicit construction to avoid block matrix (Matrix of Matrices) and ensure SparseMatrixCSC{ComplexF64, Int64}
        # Stripping all wrappers (Transpose, Adjoint) with copy and SparseMatrixCSC conversion
        block11 = sparse(1:m, 1:m, ComplexF64.(Z_self), m, m) + 1im * w[k] * SparseMatrixCSC{ComplexF64,Int64}(Lp_sp)
        block12 = spzeros(ComplexF64, m, ns)
        block13 = SparseMatrixCSC{ComplexF64,Int64}(incidence_selection[:A])

        block21 = spzeros(ComplexF64, ns, m)
        block22 = SparseMatrixCSC{ComplexF64,Int64}(P_sp)
        block23 = -SparseMatrixCSC{ComplexF64,Int64}(copy(transpose(incidence_selection[:Gamma])))

        block31 = -SparseMatrixCSC{ComplexF64,Int64}(copy(transpose(incidence_selection[:A])))
        block32 = 1im * w[k] * SparseMatrixCSC{ComplexF64,Int64}(incidence_selection[:Gamma])
        block33 = SparseMatrixCSC{ComplexF64,Int64}(Yle)

        row1 = hcat(block11, block12, block13)
        row2 = hcat(block21, block22, block23)
        row3 = hcat(block31, block32, block33)

        MatPrec_interm = vcat(row1, row2, row3)

        # CRITICAL: Reconstruct from triplets to ensure sorted, valid CSC and no structural issues
        # This is the most robust way to fix "invalid matrix" in UMFPACK
        I_tri, J_tri, V_tri = findnz(MatPrec_interm)
        MatPrec::SparseMatrixCSC{ComplexF64,Int64} = sparse(I_tri, J_tri, V_tri, size(MatPrec_interm)...)

        F = lu(MatPrec)

        # time_Lu = toc
        # println("time LU $time_Lu")

        for c1 = range(1, length(ports_to_excite))

            n1 = Int64(ports[:port_nodes][c1, 1])
            n2 = Int64(ports[:port_nodes][c1, 2])

            is[n1] = 1.0 * escalings[:Is]
            is[n2] = -1.0 * escalings[:Is]

            # tn = Q1*(U1\(L1\(P1*[zeros(m+ns,1);is])))
            rhs = [zeros(ComplexF64, m + ns); is]
            tn = F \ rhs

            # gmres call
            # GMRES_settings.tol(k) -> GMRES_settings["tol"][k]

            V, flag, relres, iter, resvec = gmres_ACA(tn, false, GMRES_settings["tol"][k], Inner_Iter, Vrest[:, c1], w[k], incidence_selection,
                P_data, Lp_data, Z_self, Yle, F, resProd, id, chan, c1)

            tot_iter_number = (iter[1] - 1) * Inner_Iter + iter[2] + 1
            if tot_iter_number < 1
                tot_iter_number = 1
            end
            println("Flag $flag - Number of iterations = $tot_iter_number")

            Vrest[:, c1] = V

            is[n1] = 0
            is[n2] = 0

            for c2 = c1:size(ports[:port_nodes], 1)
                n3 = Int64(ports[:port_nodes][c2, 1])
                n4 = Int64(ports[:port_nodes][c2, 2])

                if c1 == c2
                    S[c1, c2, k] = (2 * (V[m+ns+n3] - V[m+ns+n4]) - R_chiusura) / R_chiusura
                else
                    S[c1, c2, k] = (2 * (V[m+ns+n3] - V[m+ns+n4])) / R_chiusura
                end
                S[c2, c1, k] = S[c1, c2, k]
            end
        end

        if QS_Rcc_FW == 2
            restore_diag_P!(P_data, beta, superfici)
            restore_diag_Lp!(Lp_data, volumi, beta)
        end

    end

    out = Dict()
    out[:S] = S
    out[:Z] = s2z_S_ACA(S, 50.0) # Assuming R_chiusura usage for conversion is similar to refernece, usually 50 or port impedance
    out[:Y] = s2y_S_ACA(S, 50.0)
    out[:f] = freq ./ escalings[:freq]

    return out

end

function ComputeMatrixVectorACA(x, w, incidence_selection, P_data, Lp_data, Z_self, Yle, F, resProd)
    # x corresponds to [I; Q; Phi]
    # m, ns, n sizes

    m = size(incidence_selection[:A], 1)
    ns = size(incidence_selection[:Gamma], 2)

    mx = Lp_data[:mx]
    my = Lp_data[:my]
    mz = Lp_data[:mz]

    I = view(x, 1:m)
    Q = view(x, m+1:m+ns)
    Phi = view(x, m+ns+1:length(x))

    # Y1 = 1j*w*Lp*I + Z_self.*I + A*Phi
    # But Lp*I is done via compression.

    # We can use resProd as temp buffer if needed, but for clarity let's just compute.
    # To reduce allocs, we should prealloc Y1, Y2, Y3 or use resProd segments.
    # But for now, let's implement logic.

    Y1 = zeros(ComplexF64, m)

    # compress_matrix_times_vector
    Y1[1:mx] = compress_matrix_times_vector(Lp_data[:Lpx], Lp_data[:Lpx_comp], I[1:mx])
    Y1[mx+1:mx+my] = compress_matrix_times_vector(Lp_data[:Lpy], Lp_data[:Lpy_comp], I[mx+1:mx+my])
    Y1[mx+my+1:mx+my+mz] = compress_matrix_times_vector(Lp_data[:Lpz], Lp_data[:Lpz_comp], I[mx+my+1:mx+my+mz])

    Y1 .= 1im .* w .* Y1 .+ (Z_self .* I) .+ (incidence_selection[:A] * Phi)

    # Y2 = (P*Q) - Gamma'*Phi
    Y2_part = compress_matrix_times_vector(P_data, P_data[:P_comp], Q)
    Y2 = Y2_part .- (transpose(incidence_selection[:Gamma]) * Phi)

    # Y3 = -A'*I + Yle*Phi + 1j*w*Gamma*Q;
    Y3 = -(transpose(incidence_selection[:A]) * I) .+ (Yle * Phi) .+ (1im * w .* (incidence_selection[:Gamma] * Q))

    # MatrixVector = F \ [Y1; Y2; Y3]
    RHS = [Y1; Y2; Y3]
    MatrixVector = F \ RHS

    return MatrixVector
end


function compress_matrix_times_vector(mat_compress, comp, V)
    num_blocks = length(comp)
    res = zeros(ComplexF64, length(V))

    for c1 in 1:num_blocks
        idx = comp[c1]
        if !isempty(idx)
            res[idx] .+= mat_compress[:D][c1] * V[idx]
        end
    end

    for c1 in 1:num_blocks-1
        for c2 in c1+1:num_blocks
            idx1, idx2 = comp[c1], comp[c2]
            if !isempty(idx1) && !isempty(idx2)
                # Use isassigned to avoid UndefRefError with #undef entries
                if isassigned(mat_compress[:U], c1, c2) && isassigned(mat_compress[:V], c1, c2)
                    U_blk = mat_compress[:U][c1, c2]
                    V_blk = mat_compress[:V][c1, c2]

                    if !isnothing(U_blk) && !isnothing(V_blk)
                        res[idx1] .+= U_blk * (V_blk * V[idx2])
                        # Symmetric part (transpose)
                        res[idx2] .+= transpose((transpose(V[idx1]) * U_blk) * V_blk)
                    end
                end
            end
        end
    end
    return res
end

function remove_escalings_comp_mats!(mat_compress, comp, esca)
    num_blocks = length(comp)
    for c1 = 1:num_blocks
        if !isempty(comp[c1])
            mat_compress[:D][c1] ./= esca
        end
    end
    for c1 = 1:num_blocks-1
        for c2 = c1+1:num_blocks
            if !isempty(comp[c1]) && !isempty(comp[c2])
                mat_compress[:U][c1, c2] ./= esca
            end
        end
    end
end

function convert_comp_struct_to_complex!(comp_struct)
    # Check D
    if haskey(comp_struct, :D) && isa(comp_struct[:D], Vector{Matrix{Float64}})
        comp_struct[:D] = [ComplexF64.(d) for d in comp_struct[:D]]
    end

    # Check U and V
    for key in [:U, :V]
        if haskey(comp_struct, key)
            mat = comp_struct[key]
            if isa(mat, Matrix{Matrix{Float64}})
                rows, cols = size(mat)
                new_mat = Matrix{Matrix{ComplexF64}}(undef, rows, cols)
                for i in eachindex(mat)
                    if isassigned(mat, i)
                        new_mat[i] = ComplexF64.(mat[i])
                    end
                end
                comp_struct[key] = new_mat
            elseif isa(mat, Matrix{Union{Nothing,Matrix{Float64}}})
                rows, cols = size(mat)
                new_mat = Matrix{Union{Nothing,Matrix{ComplexF64}}}(nothing, rows, cols)
                for i in eachindex(mat)
                    if isassigned(mat, i)
                        val = mat[i]
                        if !isnothing(val)
                            new_mat[i] = ComplexF64.(val)
                        else
                            new_mat[i] = nothing
                        end
                    end
                end
                comp_struct[key] = new_mat
            end
        end
    end
end

function add_delays_P!(P_data, beta, superfici, ACA_thres, escalings, id)
    convert_comp_struct_to_complex!(P_data)

    use_suppression = 1
    num_blocks = length(P_data[:P_comp])

    for c1 = 1:num_blocks
        m = P_data[:P_comp][c1]
        if !isempty(m)
            # pdist2 in MATLAB returns distance matrix.
            # We assume a pdist2 equivalent exists or we implement it.
            # generic Julia: pairwise(Euclidean(), ...) from Distances.jl
            # OR simple broadcast logic.
            # Assuming superfici.centri is a matrix.
            dist_mat = pdist2_ACA(superfici["centri"][m, :], superfici["centri"][m, :])
            P_data[:D][c1] .*= exp.(-1im * beta * dist_mat)
        end
    end

    for c1 = 1:num_blocks-1
        for c2 = c1+1:num_blocks
            if !isempty(P_data[:P_comp][c1]) && !isempty(P_data[:P_comp][c2])
                U, V = ACA_P_delays(ACA_thres, P_data[:P_comp][c1], P_data[:P_comp][c2],
                    superfici, escalings, id, beta, use_suppression)
                P_data[:U][c1, c2] = U
                P_data[:V][c1, c2] = V
            end
        end
    end
end

function add_delays_Lp!(Lp_data, volumi, incidence_selection, beta, ACA_thres, escalings, id)
    use_suppression = 1
    mx = Lp_data[:mx]
    my = Lp_data[:my]
    mz = Lp_data[:mz]

    # X part
    centri = volumi[:centri][1:mx, :]
    process_Lp_delays!(Lp_data[:Lpx], Lp_data[:Lpx_comp], volumi, incidence_selection, centri, beta, ACA_thres, escalings, use_suppression, "x", id)

    # Y part
    centri = volumi[:centri][mx+1:mx+my, :]
    process_Lp_delays!(Lp_data[:Lpy], Lp_data[:Lpy_comp], volumi, incidence_selection, centri, beta, ACA_thres, escalings, use_suppression, "y", id)

    # Z part
    centri = volumi[:centri][mx+my+1:mx+my+mz, :]
    process_Lp_delays!(Lp_data[:Lpz], Lp_data[:Lpz_comp], volumi, incidence_selection, centri, beta, ACA_thres, escalings, use_suppression, "z", id)
end

function process_Lp_delays!(Lp_comp_struct, comp, volumi, incidence_selection, centri, beta, ACA_thres, escalings, use_suppression, component, id)
    convert_comp_struct_to_complex!(Lp_comp_struct)

    num_blocks = length(comp)
    for c1 = 1:num_blocks
        m = comp[c1]
        if !isempty(m)
            dist_mat = pdist2_ACA(centri[m, :], centri[m, :])
            Lp_comp_struct[:D][c1] .*= exp.(-1im * beta * dist_mat)
        end
    end

    for c1 = 1:num_blocks-1
        for c2 = c1+1:num_blocks
            if !isempty(comp[c1]) && !isempty(comp[c2])
                U, V = ACA_Lp_delays(ACA_thres, comp[c1], comp[c2],
                    volumi, incidence_selection, escalings, id, beta, use_suppression, component)
                Lp_comp_struct[:U][c1, c2] = U
                Lp_comp_struct[:V][c1, c2] = V
            end
        end
    end
end

function restore_diag_P!(P_data, beta, superfici)
    num_blocks = length(P_data[:P_comp])
    for c1 = 1:num_blocks
        m = P_data[:P_comp][c1]
        if !isempty(m)
            dist_mat = pdist2_ACA(superfici["centri"][m, :], superfici["centri"][m, :])
            P_data[:D][c1] ./= exp.(-1im * beta * dist_mat)
        end
    end
end

function restore_diag_Lp!(Lp_data, volumi, beta)
    mx = Lp_data[:mx]
    my = Lp_data[:my]
    mz = Lp_data[:mz]

    centri = volumi[:centri][1:mx, :]
    restore_diag_component!(Lp_data[:Lpx], Lp_data[:Lpx_comp], centri, beta)

    centri = volumi[:centri][mx+1:mx+my, :]
    restore_diag_component!(Lp_data[:Lpy], Lp_data[:Lpy_comp], centri, beta)

    centri = volumi[:centri][mx+my+1:mx+my+mz, :]
    restore_diag_component!(Lp_data[:Lpz], Lp_data[:Lpz_comp], centri, beta)
end

function restore_diag_component!(Lp_comp_struct, comp, centri, beta)
    num_blocks = length(comp)
    for c1 = 1:num_blocks
        m = comp[c1]
        if !isempty(m)
            dist_mat = pdist2_ACA(centri[m, :], centri[m, :])
            Lp_comp_struct[:D][c1] ./= exp.(-1im * beta * dist_mat)
        end
    end
end

function pdist2_ACA(A, B)
    # Simple Euclidean distance matrix between rows of A and B
    # A: n x d, B: m x d
    n = size(A, 1)
    m = size(B, 1)
    D = zeros(Float64, n, m)
    for i = 1:n
        for j = 1:m
            D[i, j] = norm(A[i, :] - B[j, :])
        end
    end
    return D
end

function s2z_S_ACA(S, Zo)
    num_ports = size(S)[1]
    nfreq = size(S)[3]
    Z = zeros(ComplexF64, num_ports, num_ports, nfreq)
    Id = Matrix{Float64}(I, num_ports, num_ports)
    for cont = 1:nfreq
        Z[:, :, cont] = Zo * ((Id - S[:, :, cont]) \ (Id + S[:, :, cont]))
    end
    return Z
end

function s2y_S_ACA(S, Zo)
    num_ports = size(S)[1]
    nfreq = size(S)[3]
    Y = zeros(ComplexF64, num_ports, num_ports, nfreq)
    Id = Matrix{Float64}(I, num_ports, num_ports)
    for cont = 1:nfreq
        Y[:, :, cont] = (1 / Zo) * ((Id + S[:, :, cont]) \ (Id - S[:, :, cont]))
    end
    return Y
end
