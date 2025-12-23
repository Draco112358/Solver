"""
    ACA_P(ACA_thres, m, n, superfici, escalingP, use_suppression)

Adaptive Cross Approximation (ACA) algorithm for P matrix compression.

This function computes a low-rank approximation U*V of a matrix block using
the ACA algorithm, which is particularly efficient for compressing dense
matrices arising from boundary element methods.

# Arguments
- `ACA_thres`: Convergence threshold for ACA compression
- `m`: Vector of row indices for the matrix block
- `n`: Vector of column indices for the matrix block
- `superfici`: Dictionary/NamedTuple containing surface data with fields:
  - `normale`: Matrix of surface normals
  - `estremi_celle`: Matrix of cell extremes
- `escalings`: Scaling factor for P matrix
- `id`: simulation id

# Returns
- `U`: Left factor of the low-rank approximation (M×k matrix)
- `V`: Right factor of the low-rank approximation (k×N matrix)
"""
function ACA_P(ACA_thres, m, n, superfici, escalings, id)
    M = length(m)
    N = length(n)

    # If P is a vector, there is nothing to compress
    if M == 1
        U = fill(escalingP, 1, 1)
        # Using the overloaded P_cpp_QS_Computation_Sym for matrix slice inputs
        V = calcola_P_noSym(m, n, superfici, escalings, id)
        return U, V
    elseif N == 1
        U = calcola_P_noSym(m, n, superfici, escalings, id)
        V = fill(escalingP, 1, 1)
        return U, V
    else
        J = zeros(Int, N)
        I = zeros(Int, M)
        ii = collect(2:M)
        jj = collect(1:N)

        # Initialization

        # Initialize the 1st row index I[1] = 1
        I[1] = 1

        # Initialize the 1st row of the approximate error matrix
        # Passing slices for row I[1]
        Rik = calcola_P_noSym(m[I[1]:I[1]], n, superfici, escalings, id)

        # Find the 1st column index J[1]
        # Rik is 1xN matrix, so access as vector or reshape
        col = findall(abs.(Rik[jj]) .== maximum(abs.(Rik[jj])))
        J[1] = jj[col[1]]
        jj = remove_index(jj, J[1])

        size_max = max(length(m), length(n))

        U = zeros(length(m), size_max)
        V = zeros(size_max, length(n))

        # First row of V
        V[1, :] = Rik / (Rik[J[1]] + eps())

        # First column of U
        # Passing slices for col J[1]
        U[:, 1] = calcola_P_noSym(m, n[J[1]:J[1]], superfici, escalings, id)

        cont_cU = 1
        cont_rV = 1

        normP = norm(U[:, cont_cU])^2 * norm(V[cont_rV, :])^2

        # Find 2nd row index I[2]
        row = findall(abs.(U[ii, 1]) .== maximum(abs.(U[ii, 1])))
        I[2] = ii[row[1]]
        ii = remove_index(ii, I[2])

        # Iteration
        for k in 2:min(M, N)

            # Compute R(I[k], :)
            row_ik_P = calcola_P_noSym(m[I[k]:I[k]], n, superfici, escalings, id)

            # U[I[k], 1:cont_cU]' * V[1:cont_rV, :] is (cU x 1)' * (rV x N) -> (1 x cU) * (rV x N). 
            # If cU == rV, it is (1 x k-1) * (k-1 x N) -> 1 x N. Correct.
            approx_row = (U[I[k], 1:cont_cU]' * V[1:cont_rV, :])
            # approx_row will be a RowVector or Matrix. row_ik_P is Matrix (1xN).
            # Transpose of approx_row if needed? No, standard matrix multiply gives row vector.

            Rik = row_ik_P - approx_row

            col = findall(abs.(Rik[jj]) .== maximum(abs.(Rik[jj])))
            J[k] = jj[col[1]]
            jj = remove_index(jj, J[k])

            # Terminate if R(I[k], J[k]) == 0
            if Rik[J[k]] == 0
                break
            end

            # Set k-th row of V equal to normalized error
            Vk = Rik / (Rik[J[k]] + eps())

            # Set k-th column of U equal to updated error
            col_jk_P = calcola_P_noSym(m, n[J[k]:J[k]], superfici, escalings, id)

            approx_col = U[:, 1:cont_cU] * V[1:cont_rV, J[k]]
            Uk = col_jk_P - approx_col

            # Norm of approximate P (ACA convergence check)
            normP = normP + 2 * sum(real.((U[:, 1:cont_cU]' * Uk) .* (Vk * V[1:cont_rV, :]')')) +
                    norm(Uk)^2 * norm(Vk)^2

            # Update U and V
            cont_cU += 1
            U[:, cont_cU] = Uk
            cont_rV += 1
            V[cont_rV, :] = Vk

            # Check convergence
            if norm(Uk) * norm(Vk) <= ACA_thres * sqrt(normP)
                break
            end

            if k == min(M, N)
                break
            end

            # Find next row index
            row = findall(abs.(Uk[ii]) .== maximum(abs.(Uk[ii])))
            I[k+1] = ii[row[1]]
            ii = remove_index(ii, I[k+1])
        end

        U = U[:, 1:cont_cU]
        V = V[1:cont_rV, :]
    end

    return U, V
end


