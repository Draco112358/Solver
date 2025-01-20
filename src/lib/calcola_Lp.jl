include("Compute_Lp_Self.jl")
include("Song_improved_Ivana_strategy.jl")

function calcola_Lp(volumi, incidence_selection, escalings, QS_Rcc_FW)
    epsilon1 = 5e-3
    epsilon2 = 1e-3
    epsilon3 = 1e-3
    epsilon4 = 3e-1

    use_suppression = true

    mx = incidence_selection["mx"]
    my = incidence_selection["my"]
    mz = incidence_selection["mz"]

    Rx = nothing
    Ry = nothing
    Rz = nothing

    volumi["centri"] = hcat(volumi["centri"]...)
    volumi["coordinate"] = hcat(volumi["coordinate"]...)

    if QS_Rcc_FW >= 2
        Rx = zeros(mx, mx)
        Threads.@threads for m in 1:mx
            for n in m:mx
                Rx[m, n] = norm(volumi["centri"][m, :] .- volumi["centri"][n, :])
                Rx[n, m] = Rx[m, n]
            end
        end

        Ry = zeros(my, my)
        Threads.@threads for m in 1:my
            for n in m:my
                Ry[m, n] = norm(volumi["centri"][m + mx, :] .- volumi["centri"][n + mx, :])
                Ry[n, m] = Ry[m, n]
            end
        end

        Rz = zeros(mz, mz)
        Threads.@threads for m in 1:mz
            for n in m:mz
                Rz[m, n] = norm(volumi["centri"][m + mx + my, :] .- volumi["centri"][n + mx + my, :])
                Rz[n, m] = Rz[m, n]
            end
        end
    end

    Lp_x = zeros(mx, mx)

    Threads.@threads for m in 1:mx
        Lp_x[m, m] = Compute_Lp_Self(volumi["coordinate"][m, :], 1) * escalings["Lp"]

        for n in m+1:mx
            integ, _ = Song_improved_Ivana_strategy(
                volumi["coordinate"][m, :], volumi["coordinate"][n, :],
                epsilon1, epsilon2, epsilon3, epsilon4, use_suppression
            )
            Lp_x[m, n] = 1e-7 / (volumi["S"][m] * volumi["S"][n]) * integ * escalings["Lp"]
            Lp_x[n, m] = Lp_x[m, n]
        end
    end

    Lp_y = zeros(my, my)

    Threads.@threads for m in 1:my
        Lp_y[m, m] = Compute_Lp_Self(volumi["coordinate"][m + mx, :], 2) * escalings["Lp"]

        for n in m+1:my
            integ, _ = Song_improved_Ivana_strategy(
                volumi["coordinate"][m + mx, :], volumi["coordinate"][n + mx, :],
                epsilon1, epsilon2, epsilon3, epsilon4, use_suppression
            )
            Lp_y[m, n] = 1e-7 / (volumi["S"][m + mx] * volumi["S"][n + mx]) * integ * escalings["Lp"]
            Lp_y[n, m] = Lp_y[m, n]
        end
    end

    Lp_z = zeros(mz, mz)

    Threads.@threads for m in 1:mz
        Lp_z[m, m] = Compute_Lp_Self(volumi["coordinate"][m + mx + my, :], 3) * escalings["Lp"]

        for n in m+1:mz
            integ, _ = Song_improved_Ivana_strategy(
                volumi["coordinate"][m + mx + my, :], volumi["coordinate"][n + mx + my, :],
                epsilon1, epsilon2, epsilon3, epsilon4, use_suppression
            )
            Lp_z[m, n] = 1e-7 / (volumi["S"][m + mx + my] * volumi["S"][n + mx + my]) * integ * escalings["Lp"]
            Lp_z[n, m] = Lp_z[m, n]
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
