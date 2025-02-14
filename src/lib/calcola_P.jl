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
    println("Calcolo P initialization")
    if QS_Rcc_FW >= 2
        R_cc = zeros(nsup, nsup)
        block_size1 = 500  # ad esempio, 10 iterazioni per blocco

        Threads.@threads for m_block in 1:block_size1:nsup
            m_end = min(m_block + block_size1 - 1, nsup)
            for m in m_block:m_end
                for n in m:nsup
                    #norm = norm(superfici["centri"][m, :] - superfici["centri"][n, :])
                    dist = norm(view(superfici["centri"], m, :) .- view(superfici["centri"], n, :))
                    R_cc[m, n] = dist
                    R_cc[n, m] = dist
                end
            end
            # Dopo aver processato un blocco, cedi il controllo per permettere la gestione dei heartbeat e altre operazioni
            sleep(0)
            println("block: ", round(m_end/block_size1), " / ", round(nsup/block_size1))
        end
        # Threads.@threads for m in 1:nsup
        #     for n in m:nsup
        #         R_cc[m, n] = norm(superfici["centri"][m, :] - superfici["centri"][n, :])
        #         R_cc[n, m] = R_cc[m, n]
        #     end
        # end
    end
    P = zeros(nsup, nsup)
    println("Calcolo P Song_P_improved_Ivana_strategy")

    # Scegli una dimensione di blocco adatta (da regolare in base alle tue esigenze)
    block_size2 = 500  # ad esempio, processa 1000 iterazioni per blocco

    Threads.@threads for m_block in 1:block_size2:nsup
        m_end = min(m_block + block_size2 - 1, nsup)
        for m in m_block:m_end
            for n in m:nsup
                integ, _ = Song_P_improved_Ivana_strategy(
                    superfici["estremi_celle"][m, :],
                    superfici["estremi_celle"][n, :],
                    epsilon1, epsilon2, epsilon3, epsilon4,
                    use_suppression)
                P[m, n] = 1 / (4 * π * eps0 * superfici["S"][m] * superfici["S"][n]) *
                        integ * escalings[:P]
                P[n, m] = P[m, n]
            end
        end
        # Al termine di ogni blocco, cedi il controllo per consentire l'invio dei heartbeat
        sleep(0)
        println("block: ", round(m_end/block_size2), " / ", round(nsup/block_size1))
    end


    # Threads.@threads for m in 1:nsup
    #     for n in m:nsup
    #         println(" m:Nsup ", m:nsup)
    #         integ, _ = Song_P_improved_Ivana_strategy(superfici["estremi_celle"][m, :], superfici["estremi_celle"][n, :], epsilon1, epsilon2, epsilon3, epsilon4, use_suppression)
    #         P[m, n] = 1 / (4 * π * eps0 * superfici["S"][m] * superfici["S"][n]) * integ * escalings[:P]
    #         P[n, m] = P[m, n]
    #     end
    # end

    return Dict(
        :P => P,
        :R_cc => R_cc
    )
end
