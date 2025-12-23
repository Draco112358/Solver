"""
    ACA_Lp(ACA_thres, m, n, estremi_celle, versori, centri, escalings_Lp, use_suppression)

Adaptive Cross Approximation (ACA) algorithm for Lp matrix compression.
"""
function ACA_Lp(ACA_thres, m, n, volumi, incidence_selection, escalings, id, use_suppression, component::String="all")
    M = length(m)
    N = length(n)

    escalings_Lp = escalings[:Lp]

    if M == 1
        U = fill(Float64(escalings_Lp), 1, 1)
        V = calcola_Lp_noSym(m, n, volumi, incidence_selection, escalings, id, component)
        return U, V
    elseif N == 1
        U = calcola_Lp_noSym(m, n, volumi, incidence_selection, escalings, id, component)
        V = fill(Float64(escalings_Lp), 1, 1)
        return U, V
    else
        J = zeros(Int, N)
        I = zeros(Int, M)
        ii = collect(2:M)
        jj = collect(1:N)

        I[1] = 1

        # R(I[1], :)
        Rik = calcola_Lp_noSym(m[I[1]:I[1]], n, volumi, incidence_selection, escalings, id, component)

        col = findall(abs.(Rik[jj]) .== maximum(abs.(Rik[jj])))
        J[1] = jj[col[1]]
        jj = remove_index(jj, J[1])

        size_max = max(M, N)
        U = zeros(M, size_max)
        V = zeros(size_max, N)

        V[1, :] = Rik / (Rik[J[1]] + eps())

        # U(:, 1)
        U[:, 1] = calcola_Lp_noSym(m, n[J[1]:J[1]], volumi, incidence_selection, escalings, id, component)

        cont_cU = 1
        cont_rV = 1
        normP = norm(U[:, 1])^2 * norm(V[1, :])^2

        row = findall(abs.(U[ii, 1]) .== maximum(abs.(U[ii, 1])))
        I[2] = ii[row[1]]
        ii = remove_index(ii, I[2])

        for k in 2:min(M, N)
            row_ik_P = calcola_Lp_noSym(m[I[k]:I[k]], n, volumi, incidence_selection, escalings, id, component)

            approx_row = (U[I[k], 1:cont_cU]' * V[1:cont_rV, :])
            Rik = row_ik_P - approx_row

            col = findall(abs.(Rik[jj]) .== maximum(abs.(Rik[jj])))
            J[k] = jj[col[1]]
            jj = remove_index(jj, J[k])

            if Rik[J[k]] == 0
                break
            end

            Vk = Rik / (Rik[J[k]] + eps())

            col_jk_P = calcola_Lp_noSym(m, n[J[k]:J[k]], volumi, incidence_selection, escalings, id, component)

            approx_col = U[:, 1:cont_cU] * V[1:cont_rV, J[k]]
            Uk = col_jk_P - approx_col

            normP = normP + 2 * sum(real.((U[:, 1:cont_cU]' * Uk) .* (Vk * V[1:cont_rV, :]')')) +
                    norm(Uk)^2 * norm(Vk)^2

            cont_cU += 1
            U[:, cont_cU] = Uk
            cont_rV += 1
            V[cont_rV, :] = Vk

            if norm(Uk) * norm(Vk) <= ACA_thres * sqrt(normP)
                break
            end

            if k == min(M, N)
                break
            end

            row = findall(abs.(Uk[ii]) .== maximum(abs.(Uk[ii])))
            I[k+1] = ii[row[1]]
            ii = remove_index(ii, I[k+1])
        end

        U = U[:, 1:cont_cU]
        V = V[1:cont_rV, :]
    end
    return U, V
end

function remove_index(ii, e)
    filter!(x -> x != e, ii)
    return ii
end
