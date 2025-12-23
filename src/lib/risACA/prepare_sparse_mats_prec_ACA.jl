using SparseArrays

"""
    prepare_sparse_mats_prec_ACA(P_data::Dict, Lp_data::Dict)

Translates the ACA compressed data structures into sparse preconditioner matrices
by extracting the diagonal blocks.
"""
function prepare_sparse_mats_prec_ACA(P_data::Dict, Lp_data::Dict)

    # --- P_sp Construction ---
    ns = P_data[:ns]

    # Calculate total number of elements to preallocate
    nElemTot_P = sum(length(idx)^2 for idx in P_data[:P_comp])

    # Determine the type for Vi based on the first non-empty diagonal block
    # Or just use promote_type of all diagonal blocks. For simplicity, we can look at the first non-empty block.
    # Actually, better yet, use eltype(P_data[:D][findfirst(!isempty, P_data[:P_comp])])

    first_idx = findfirst(!isempty, P_data[:P_comp])
    T_P = !isnothing(first_idx) ? eltype(P_data[:D][first_idx]) : Float64

    Ii = zeros(Int, nElemTot_P)
    Ji = zeros(Int, nElemTot_P)
    Vi = zeros(T_P, nElemTot_P)

    pos = 0
    for kcon in 1:length(P_data[:P_comp])
        idx = P_data[:P_comp][kcon]
        n_blk = length(idx)
        if n_blk == 0
            continue
        end

        # Diagonal block matrix
        blk_mat = P_data[:D][kcon]

        # Fill coordinates
        for j in 1:n_blk
            for i in 1:n_blk
                pos += 1
                Ii[pos] = idx[i]
                Ji[pos] = idx[j]
                Vi[pos] = blk_mat[i, j]
            end
        end
    end

    P_sp = sparse(Ii, Ji, Vi, ns, ns)

    # --- Lp_sp Construction ---
    mx = Lp_data[:mx]
    my = Lp_data[:my]
    mz = Lp_data[:mz]

    # Helper to collect sparse triplets for Lp components
    function collect_triplets(comp_indices, blocks_data, offset)
        nElem = sum(length(idx)^2 for idx in comp_indices)

        # Determine type
        first_idx = findfirst(!isempty, comp_indices)
        T_Lp = !isnothing(first_idx) ? eltype(blocks_data[first_idx]) : Float64

        I = zeros(Int, nElem)
        J = zeros(Int, nElem)
        V = zeros(T_Lp, nElem)

        p = 0
        for k in 1:length(comp_indices)
            idx = comp_indices[k]
            n_blk = length(idx)
            if n_blk == 0
                continue
            end

            blk_mat = blocks_data[k]

            for j in 1:n_blk
                for i in 1:n_blk
                    p += 1
                    I[p] = idx[i] + offset
                    J[p] = idx[j] + offset
                    V[p] = blk_mat[i, j]
                end
            end
        end
        return I[1:p], J[1:p], V[1:p]
    end

    Iix, Jix, Vix = collect_triplets(Lp_data[:Lpx_comp], Lp_data[:Lpx][:D], 0)
    Iiy, Jiy, Viy = collect_triplets(Lp_data[:Lpy_comp], Lp_data[:Lpy][:D], mx)
    Iiz, Jiz, Viz = collect_triplets(Lp_data[:Lpz_comp], Lp_data[:Lpz][:D], mx + my)

    # Combine components
    I_all = [Iix; Iiy; Iiz]
    J_all = [Jix; Jiy; Jiz]
    V_all = [Vix; Viy; Viz]

    Lp_sp = sparse(I_all, J_all, V_all, mx + my + mz, mx + my + mz)

    return P_sp, Lp_sp
end

"""
    crea_vect_s(sizemat, comp, mats)

Helper matching MATLAB implementation (defined but unused in original main body).
"""
function crea_vect_s(sizemat, comp, mats)
    nElemTot = sum(length(x)^2 for x in comp)

    Ii = zeros(Int, nElemTot)
    Ji = zeros(Int, nElemTot)
    Vi = zeros(Float64, nElemTot)

    pos = 0
    for k in 1:length(comp)
        idx = comp[k]
        n = length(idx)
        for j in 1:n
            for i in 1:n
                pos += 1
                Ii[pos] = idx[i]
                Ji[pos] = idx[j]
                # Assuming mats is a dense matrix being sampled
                Vi[pos] = mats[idx[i], idx[j]]
            end
        end
    end

    return Ii, Ji, Vi
end
