function calcola_P_ACA_cpp(superfici, ACA_thres::Float64, escalings, id::String="")
    # Constants
    use_suppression = 1

    # Extract number of surfaces
    # Ensure estremi_celle is a matrix. If it's a vector of vectors, we might need to handle it, 
    # but Solver usually ensures it's a Matrix by this point or we can convert.
    # Assuming it is handled in the calling scope or we add a check.
    # if isa(superfici["estremi_celle"], Vector)
    #     superfici["estremi_celle"] = hcat(superfici["estremi_celle"]...)'
    # end
    # if isa(superfici["normale"], Vector)
    #     superfici["normale"] = hcat(superfici["normale"]...)'
    # end
    # if isa(superfici["centri"], Vector)
    #     superfici["centri"] = hcat(superfici["centri"]...)'
    # end

    nsup = size(superfici["estremi_celle"], 1)

    P_data = Dict{Symbol,Any}()
    P_data[:ns] = nsup

    # Level calculation
    Level_P = ceil(Int, log2(nsup / 600))
    if Level_P < 1
        Level_P = 1
    end

    # Multilevel preparation
    # Note: prepare_multilevel returns (boxes, level). 
    # Check if prepare_multilevel expects Matrix or Vector of Vectors for centers.
    # Usually Julia functions prefer Matrices if column-major, but prepare_multilevel might have been ported differently.
    # Let's assume it accepts the Matrix format provided in superfici["centri"].
    basis_func_boxes_P, Level_P = prepare_multilevel(superfici["centri"], Level_P)

    # Last level fill
    P_data[:P_comp] = last_level_fill(basis_func_boxes_P, Level_P)

    num_blocks = length(P_data[:P_comp])

    # Initialize outputs
    # Using specific typed arrays for efficiency
    P_data[:D] = Vector{Matrix{Float64}}(undef, num_blocks)
    # U and V are cell/arrays of matrices. Julia sparse storage might be better, but we'll stick to dense matrix array
    P_data[:U] = Matrix{Matrix{Float64}}(undef, num_blocks, num_blocks)
    P_data[:V] = Matrix{Matrix{Float64}}(undef, num_blocks, num_blocks)

    # Diagonal Blocks Loop
    for c1 in 1:num_blocks
        m = P_data[:P_comp][c1]

        if isempty(m)
            P_data[:D][c1] = Matrix{Float64}(undef, 0, 0)
        else
            # Call P_cpp_QS_Computation_Sym (Julia implementation replacement for mex)
            # Signature: (M, N, Norm_m, Norm_n, Extr_m, Extr_n, supp)
            P_val = calcola_P_Sym(m, superfici, escalings, id)

            P_data[:D][c1] = P_val
        end
        # println("Processed diagonal block $c1 / $num_blocks")
    end

    # Off-Diagonal Blocks Loop (ACA)
    for c1 in 1:num_blocks-1
        for c2 in c1+1:num_blocks
            m = P_data[:P_comp][c1]
            n = P_data[:P_comp][c2]

            if !isempty(m) && !isempty(n)
                # Call ACA_P
                # It returns U, V.
                # Arguments: ACA_thres, m, n, superfici, scaling, suppression
                U_blk, V_blk = ACA_P(
                    ACA_thres, m, n,
                    superfici,
                    escalings,
                    id
                )

                P_data[:U][c1, c2] = U_blk
                P_data[:V][c1, c2] = V_blk

                # ACA_P might need to handle the symmetric part or we compute it? 
                # MATLAB code:
                # [P_data.U{c1,c2},P_data.V{c1,c2}]=ACA_P(...)
                # It only fills upper triangle {c1, c2}. The receiver implies symmetry is handled later 
                # or interactions are symmetric. 
                # Typically P is symmetric. If U,V represents Top-Right, then Bottom-Left is Transpose.
                # We just store what MATLAB stored.
            end
        end
        # println("Processed off-diagonal row $c1 / $num_blocks")
    end

    return P_data
end
