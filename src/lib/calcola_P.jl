include("Song_P_improved_Ivana_strategy.jl")

function calcola_P(superfici, escalings, QS_Rcc_FW)

    eps0 = 8.854187816997944e-12

    epsilon1 = 5e-3
    epsilon2 = 1e-3
    epsilon3 = 1e-3
    epsilon4 = 3e-1

    use_suppression = 1
    superfici["estremi_celle"] = hcat(superfici["estremi_celle"]...)
    superfici["centri"] = hcat(superfici["centri"]...)
    nsup = size(superfici["estremi_celle"], 1)
    R_cc = []
    if QS_Rcc_FW >= 2
        R_cc = zeros(nsup, nsup)
        for m in 1:nsup
            for n in m:nsup
                R_cc[m, n] = norm(superfici["centri"][m, :] - superfici["centri"][n, :])
                R_cc[n, m] = R_cc[m, n]
            end
        end
    end

    P = zeros(nsup, nsup)

    for m in 1:nsup
        for n in m:nsup
            integ, _ = Song_P_improved_Ivana_strategy(superfici["estremi_celle"][m, :], superfici["estremi_celle"][n, :], epsilon1, epsilon2, epsilon3, epsilon4, use_suppression)
            P[m, n] = 1 / (4 * π * eps0 * superfici["S"][m] * superfici["S"][n]) * integ * escalings["P"]
            P[n, m] = P[m, n]
        end
    end

    return Dict(
        :P => P,
        :R_cc => R_cc
    )
end