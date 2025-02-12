include("Compute_Lp_Self2.jl")
include("Song_improved_Ivana_strategy2.jl")

function calcola_Lp2(volumi, incidence_selection, escalings, QS_Rcc_FW)
    epsilon1 = 5e-3
    epsilon2 = 1e-3
    epsilon3 = 1e-3
    epsilon4 = 3e-1

    use_suppression = true

    mx = incidence_selection[:mx]
    my = incidence_selection[:my]
    mz = incidence_selection[:mz]

    Rx = nothing
    Ry = nothing
    Rz = nothing

    if QS_Rcc_FW >= 2
        Rx = Matrix{Float64}(undef, mx, mx)
        Ry = Matrix{Float64}(undef, my, my)
        Rz = Matrix{Float64}(undef, mz, mz)

        Threads.@threads for m in 1:mx
            for n in m:mx
                dist = norm(volumi[:centri][m, :] .- volumi[:centri][n, :])
                Rx[m, n] = dist
                Rx[n, m] = dist
            end
        end

        Threads.@threads for m in 1:my
            for n in m:my
                dist = norm(volumi[:centri][m + mx, :] .- volumi[:centri][n + mx, :])
                Ry[m, n] = dist
                Ry[n, m] = dist
            end
        end

        Threads.@threads for m in 1:mz
            for n in m:mz
                dist = norm(volumi[:centri][m + mx + my, :] .- volumi[:centri][n + mx + my, :])
                Rz[m, n] = dist
                Rz[n, m] = dist
            end
        end
    end

    Lp_x = Matrix{Float64}(undef, mx, mx)
    Lp_y = Matrix{Float64}(undef, my, my)
    Lp_z = Matrix{Float64}(undef, mz, mz)

    Threads.@threads for m in 1:mx
        Lp_x[m, m] = Compute_Lp_Self2(volumi[:coordinate][m, :], 1) * escalings[:Lp]

        for n in m+1:mx
            integ, _ = Song_improved_Ivana_strategy2(
                volumi[:coordinate][m, :], volumi[:coordinate][n, :],
                epsilon1, epsilon2, epsilon3, epsilon4, use_suppression
            )
            value = 1e-7 / (volumi[:S][m] * volumi[:S][n]) * integ * escalings[:Lp]
            Lp_x[m, n] = value
            Lp_x[n, m] = value
        end
    end

    Threads.@threads for m in 1:my
        Lp_y[m, m] = Compute_Lp_Self2(volumi[:coordinate][m + mx, :], 2) * escalings[:Lp]

        for n in m+1:my
            integ, _ = Song_improved_Ivana_strategy2(
                volumi[:coordinate][m + mx, :], volumi[:coordinate][n + mx, :],
                epsilon1, epsilon2, epsilon3, epsilon4, use_suppression
            )
            value = 1e-7 / (volumi[:S][m + mx] * volumi[:S][n + mx]) * integ * escalings[:Lp]
            Lp_y[m, n] = value
            Lp_y[n, m] = value
        end
    end

    Threads.@threads for m in 1:mz
        Lp_z[m, m] = Compute_Lp_Self2(volumi[:coordinate][m + mx + my, :], 3) * escalings[:Lp]

        for n in m+1:mz
            integ, _ = Song_improved_Ivana_strategy2(
                volumi[:coordinate][m + mx + my, :], volumi[:coordinate][n + mx + my, :],
                epsilon1, epsilon2, epsilon3, epsilon4, use_suppression
            )
            value = 1e-7 / (volumi[:S][m + mx + my] * volumi[:S][n + mx + my]) * integ * escalings[:Lp]
            Lp_z[m, n] = value
            Lp_z[n, m] = value
        end
    end
    
    return Dict(
        :Lp_x => Lp_x,
        :Lp_y => Lp_y,
        :Lp_z => Lp_z,
        :Rx => Rx,
        :Ry => Ry,
        :Rz => Rz
    )
end
