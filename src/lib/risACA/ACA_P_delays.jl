using LinearAlgebra

"""
    ACA_P_delays(ACA_thres, m, n, superfici, escalings, id, beta, use_suppression)

Adaptive Cross Approximation (ACA) algorithm for P matrix compression with delays (Full-Wave).
"""
function ACA_P_delays(ACA_thres, m, n, superfici, escalings, id, beta, use_suppression)
    M = length(m)
    N = length(n)

    # If P is a vector, there is nothing to compress
    if M == 1
        U = fill(1.0, 1, 1)
        # Compute base values
        V_base = calcola_P_noSym(m, n, superfici, escalings, id)
        # Apply delays
        dist = [norm(superfici["centri"][m[1], :] - superfici["centri"][n[j], :]) for j in 1:N]
        V = V_base .* exp.(-1im * beta * dist')
        return U, V
    elseif N == 1
        # Compute base values
        U_base = calcola_P_noSym(m, n, superfici, escalings, id)
        # Apply delays
        dist = [norm(superfici["centri"][m[i], :] - superfici["centri"][n[1], :]) for i in 1:M]
        U = U_base .* exp.(-1im * beta * dist)
        V = fill(1.0, 1, 1)
        return U, V
    else
        J = zeros(Int, N)
        I = zeros(Int, M)
        ii = collect(2:M)
        jj = collect(1:N)

        I[1] = 1

        # Initial row with delays
        Rik_base = calcola_P_noSym(m[I[1]:I[1]], n, superfici, escalings, id)
        dist_row = [norm(superfici["centri"][m[I[1]], :] - superfici["centri"][n[j], :]) for j in 1:N]
        Rik = Rik_base .* exp.(-1im * beta * dist_row')

        col_idx = findmax(abs.(Rik[jj]))[2]
        J[1] = jj[col_idx]
        jj = filter(x -> x != J[1], jj)

        size_max = max(M, N)
        U = zeros(ComplexF64, M, size_max)
        V = zeros(ComplexF64, size_max, N)

        V[1, :] = Rik / (Rik[J[1]] + eps())

        # Initial column with delays
        U_col_base = calcola_P_noSym(m, n[J[1]:J[1]], superfici, escalings, id)
        dist_col = [norm(superfici["centri"][m[i], :] - superfici["centri"][n[J[1]], :]) for i in 1:M]
        U[:, 1] = U_col_base .* exp.(-1im * beta * dist_col)

        cont_cU = 1
        cont_rV = 1

        normP = norm(U[:, 1])^2 * norm(V[1, :])^2

        row_idx = findmax(abs.(U[ii, 1]))[2]
        I[2] = ii[row_idx]
        ii = filter(x -> x != I[2], ii)

        for k in 2:min(M, N)
            # Row k with delays
            row_k_base = calcola_P_noSym(m[I[k]:I[k]], n, superfici, escalings, id)
            dist_k = [norm(superfici["centri"][m[I[k]], :] - superfici["centri"][n[j], :]) for j in 1:N]
            row_k = row_k_base .* exp.(-1im * beta * dist_k')

            # Note: Rik should be 1xN.
            Rik = row_k - reshape(U[I[k], 1:cont_cU], 1, cont_cU) * V[1:cont_rV, :]

            col_idx = findmax(abs.(Rik[jj]))[2]
            J[k] = jj[col_idx]
            jj = filter(x -> x != J[k], jj)

            if Rik[J[k]] == 0
                break
            end

            Vk = Rik / (Rik[J[k]] + eps())

            # Column k with delays
            col_k_base = calcola_P_noSym(m, n[J[k]:J[k]], superfici, escalings, id)
            dist_ck = [norm(superfici["centri"][m[i], :] - superfici["centri"][n[J[k]], :]) for i in 1:M]
            col_k = col_k_base .* exp.(-1im * beta * dist_ck)

            Uk = col_k - U[:, 1:cont_cU] * V[1:cont_rV, J[k]]

            # Update normP

            normP += 2 * sum(real.((U[:, 1:cont_cU]' * Uk) .* (Vk * V[1:cont_rV, :]')')) + norm(Uk)^2 * norm(Vk)^2

            cont_cU += 1
            U[:, cont_cU] = Uk
            cont_rV += 1
            V[cont_rV, :] = Vk

            if norm(Uk) * norm(Vk) <= ACA_thres * sqrt(norm(normP)) # MATLAB uses sqrt(normP)
                break
            end

            if k == min(M, N)
                break
            end

            row_idx = findmax(abs.(Uk[ii]))[2]
            I[k+1] = ii[row_idx]
            ii = filter(x -> x != I[k+1], ii)
        end

        U = U[:, 1:cont_cU]
        V = V[1:cont_rV, :]
    end

    return U, V
end
