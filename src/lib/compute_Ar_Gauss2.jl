using LinearAlgebra
using StaticArrays # Essenziale per prestazioni elevate con array di piccole dimensioni
using Base.Threads
using Printf
using StaticArrays

function compute_Ar_Gauss(barre::Transpose{Float64, Matrix{Float64}},
                                    centriOss::Matrix{Float64}, ordine::Int, beta::ComplexF64,
                                    simulation_id=nothing, chan=nothing)

    numCentri = size(centriOss, 1)
    numBarre = size(barre, 1)

    ha = zeros(ComplexF64, numBarre, numCentri)

    # --- Pre-calcoli globali (una volta per la funzione) ---

    # 1. qrule values (già fuori, ottimo)
    rootkx, wekx = qrule(ordine)
    rootky, weky = qrule(ordine)
    rootkz, wekz = qrule(ordine)
    nlkx = length(wekx)
    nlky = length(weky)
    nlkz = length(wekz)

    # 2. Pre-calcola le coordinate x_o, y_o, z_o per tutti i centri di osservazione
    # Convertiamo direttamente in SVector{3,Float64} per comodità
    centriOss_svec = Vector{SVector{3,Float64}}(undef, numCentri)
    @inbounds Threads.@threads for cc = 1:numCentri
        centriOss_svec[cc] = SVector{3,Float64}(centriOss[cc, 1], centriOss[cc, 2], centriOss[cc, 3])
    end

    # Questi array conterranno i 3 componenti (x,y,z) di ogni 'r' vettore per ciascuna barra
    rmi_arr = Vector{SVector{3,Float64}}(undef, numBarre)
    rai_arr = Vector{SVector{3,Float64}}(undef, numBarre)
    rbi_arr = Vector{SVector{3,Float64}}(undef, numBarre)
    rci_arr = Vector{SVector{3,Float64}}(undef, numBarre)
    rabi_arr = Vector{SVector{3,Float64}}(undef, numBarre)
    rbci_arr = Vector{SVector{3,Float64}}(undef, numBarre)
    raci_arr = Vector{SVector{3,Float64}}(undef, numBarre)
    rabci_arr = Vector{SVector{3,Float64}}(undef, numBarre)

    @inbounds Threads.@threads for cont = 1:numBarre
        barra_svec_24 = SVector{24, Float64}(@view barre[cont, :])
        # Calcola e memorizza tutti gli 8 vettori 'r' (rmi, rai, ..., rabci)
        # Re-inseriamo i calcoli completi basati sulle espressioni originali.
        # Definiamo i vertici in modo strutturato per chiarezza.
        # r1, r2, ..., r8 corrispondono a barra_svec_24 elementi.
        # Vertici: (x,y,z)
        v1 = SVector{3,Float64}(barra_svec_24[1], barra_svec_24[2], barra_svec_24[3])
        v2 = SVector{3,Float64}(barra_svec_24[4], barra_svec_24[5], barra_svec_24[6])
        v3 = SVector{3,Float64}(barra_svec_24[7], barra_svec_24[8], barra_svec_24[9])
        v4 = SVector{3,Float64}(barra_svec_24[10], barra_svec_24[11], barra_svec_24[12])
        v5 = SVector{3,Float64}(barra_svec_24[13], barra_svec_24[14], barra_svec_24[15])
        v6 = SVector{3,Float64}(barra_svec_24[16], barra_svec_24[17], barra_svec_24[18])
        v7 = SVector{3,Float64}(barra_svec_24[19], barra_svec_24[20], barra_svec_24[21])
        v8 = SVector{3,Float64}(barra_svec_24[22], barra_svec_24[23], barra_svec_24[24])

        rmi_arr[cont] = 0.125 * (v1 + v2 + v3 + v4 + v5 + v6 + v7 + v8)
        rai_arr[cont] = 0.125 * (-v1 + v2 - v3 + v4 - v5 + v6 - v7 + v8)
        rbi_arr[cont] = 0.125 * (-v1 - v2 + v3 + v4 - v5 - v6 + v7 + v8)
        rci_arr[cont] = 0.125 * (-v1 - v2 - v3 - v4 + v5 + v6 + v7 + v8)
        rabi_arr[cont] = 0.125 * (v1 - v2 + v3 - v4 + v5 - v6 + v7 - v8)
        rbci_arr[cont] = 0.125 * (v1 + v2 - v3 - v4 - v5 - v6 + v7 + v8)
        raci_arr[cont] = 0.125 * (v1 - v2 - v3 + v4 - v5 + v6 + v7 - v8)
        rabci_arr[cont] = 0.125 * (-v1 + v2 - v3 + v4 + v5 - v6 + v7 - v8)
    end


    # --- Loop di elaborazione parallelizzata ---
    block_size = 200 # Puoi sperimentare con questo valore

    @inbounds for m_block_start in 1:block_size:numBarre
        m_end = min(m_block_start + block_size - 1, numBarre)
        Threads.@threads for cont in m_block_start:m_end

            # Accedi ai coefficienti 'r' pre-calcolati per la barra corrente
            current_rmi = rmi_arr[cont]
            current_rai = rai_arr[cont]
            current_rbi = rbi_arr[cont]
            current_rci = rci_arr[cont]
            current_rabi = rabi_arr[cont]
            current_rbci = rbci_arr[cont]
            current_raci = raci_arr[cont]
            current_rabci = rabci_arr[cont]

            for cc in 1:numCentri
                # Accedi alle coordinate pre-calcolate del centro di osservazione (già SVector)
                origin_vec = centriOss_svec[cc]

                ha[cont, cc] = compute_ha_optimized_corrected(
                    origin_vec, # Passiamo l'SVector direttamente
                    current_rmi, # Passiamo l'SVector direttamente
                    current_rai,
                    current_rbi,
                    current_rci,
                    current_rabi,
                    current_rbci,
                    current_raci,
                    current_rabci,
                    beta, rootkx, wekx, rootky, weky, rootkz, wekz,
                    nlkx, nlky, nlkz
                )
            end
        end
        #Assicurati che @printf sia disponibile (di solito con using Printf)
        @printf "Block Ar : %.0f / %.0f\n" round(m_end / block_size) round(numBarre / block_size)
        if is_stop_requested(simulation_id)
            println("Simulazione $(simulation_id) interrotta per richiesta stop.")
            return nothing
        end
    end

    return ha
end

@inline function compute_ha_optimized_corrected(
    orig::SVector{3, Float64}, # xo, yo, zo
    rmi::SVector{3, Float64},
    rai::SVector{3, Float64},
    rbi::SVector{3, Float64},
    rci::SVector{3, Float64},
    rabi::SVector{3, Float64},
    rbci::SVector{3, Float64},
    raci::SVector{3, Float64},
    rabci::SVector{3, Float64},
    beta::ComplexF64,
    rootkx::Vector{Float64}, wekx::Vector{Float64},
    rootky::Vector{Float64}, weky::Vector{Float64},
    rootkz::Vector{Float64}, wekz::Vector{Float64},
    nlkx::Int64, nlky::Int64, nlkz::Int64
)
    sum_a1 = zero(ComplexF64)
    @fastmath for a1 in 1:nlkx
        sum_b1 = zero(ComplexF64)
        kx_a1 = rootkx[a1]

        for b1 in 1:nlky
            sum_c1 = zero(ComplexF64)
            ky_b1 = rootky[b1]

            for c1 in 1:nlkz
                kz_c1 = rootkz[c1]

                # Calcolo vettoriale per drai, drbi, drci
                drai = rai + rabi * ky_b1 + raci * kz_c1 + rabci * ky_b1 * kz_c1
                draim = norm(drai) # sqrt(drai.x^2 + drai.y^2 + drai.z^2)

                drbi = rbi + rabi * kx_a1 + rbci * kz_c1 + rabci * kx_a1 * kz_c1
                # Attenzione qui: c'era un errore nel tuo codice originale in drbi_z.
                # Ho assunto che dovesse seguire lo schema degli altri componenti.
                # drbi_z = rbi_z + rabi_z * kx_a1 + rbci_z * kz_c1 + rabci_z * kx_a1 * kz_c1
                drbim = norm(drbi)

                drci = rci + raci * kx_a1 + rbci * ky_b1 + rabci * kx_a1 * ky_b1
                drcim = norm(drci)

                # Calcolo vettoriale per r1
                r1 = rmi + rai * kx_a1 + rbi * ky_b1 + rci * kz_c1 +
                     rabi * kx_a1 * ky_b1 + raci * kx_a1 * kz_c1 + rbci * ky_b1 * kz_c1 +
                     rabci * kx_a1 * ky_b1 * kz_c1

                delta = orig - r1
                R = norm(delta)

                invR = 1.0 / R
                G = invR * exp(-1im * beta * R)

                f = (draim * drbim) * (drcim * G)
                sum_c1 += wekz[c1] * f
            end
            sum_b1 += weky[b1] * sum_c1
        end
        sum_a1 += wekx[a1] * sum_b1
    end

    return sum_a1
end