function compute_lambda_numeric(punti_oss::Matrix{Float64}, volumi::Dict, incidence_selection::Dict{Symbol, Any},
                                          vers_punti_oss::Matrix{Float64}, ordine_int::Int, beta::ComplexF64,
                                          id=nothing, chan=nothing)

    volumi[:coordinate] = transpose(Float64.(volumi[:coordinate].parent))
    N = size(volumi[:coordinate], 1) # Number of volumes
    M = size(punti_oss, 1) # Number of observation points

    # --- Inizializzazioni e pre-calcoli globali (una volta per funzione) ---

    # 1. Pre-allocate vers_celle (Già fuori dal ciclo)
    vers_celle = zeros(Float64, N, 3)
    @inbounds begin
        vers_celle[1:incidence_selection[:mx], 1] .= 1.0
        vers_celle[incidence_selection[:mx]+1:incidence_selection[:mx]+incidence_selection[:my], 2] .= 1.0
        vers_celle[incidence_selection[:mx]+incidence_selection[:my]+1:N, 3] .= 1.0
    end

    # 2. Pre-calculate qrule values once (Già fuori dal ciclo)
    rootk, wek = qrule_Ec(ordine_int)

    # 3. Pre-calcola barra_n_svec_arr e S_n_arr per TUTTI i volumi (N)
    # Spostato FUORI dal Threads.@threads loop.
    # Questi array saranno accessibili da tutti i thread.
    barra_n_svec_arr = Vector{SVector{24, Float64}}(undef, N)
    S_n_arr = Vector{Float64}(undef, N)
    @inbounds Threads.@threads for n = 1:N
        barra_n_svec_arr[n] = SVector{24}(@view volumi[:coordinate][n, :])
        S_n_arr[n] = volumi[:S][n]
    end

    # 4. Pre-calcola punti_oss_m_svec, vers_punti_oss_m_svec, is_pox, is_poy per TUTTI i punti di osservazione (M)
    # Spostato FUORI dal Threads.@threads loop.
    punti_oss_m_svec_arr = Vector{SVector{3, Float64}}(undef, M)
    vers_punti_oss_m_svec_arr = Vector{SVector{3, Float64}}(undef, M)
    is_pox_arr = Vector{Bool}(undef, M)
    is_poy_arr = Vector{Bool}(undef, M)

    @inbounds Threads.@threads for m = 1:M
        punti_oss_m_svec_arr[m] = SVector{3}(@view punti_oss[m, :])
        vers_punti_oss_m_svec_arr[m] = SVector{3}(@view vers_punti_oss[m, :])
        is_pox_arr[m] = abs(vers_punti_oss_m_svec_arr[m][1]) > 1e-5
        is_poy_arr[m] = abs(vers_punti_oss_m_svec_arr[m][2]) > 1e-5
    end

    # Pre-allocazione della matrice finale Lambda
    Lambda = zeros(ComplexF64, M, N)

    # --- Blocco di elaborazione parallelizzata ---
    # L'iterazione esterna parallelizzata ora opera su `m`, ma accede a dati pre-calcolati.

    block_size_M = 10 # Adjust based on your system and problem size

    for m_block_start in 1:block_size_M:M
        m_block_end = min(m_block_start + block_size_M - 1, M)

        # Iterazione interna sui punti di osservazione all'interno del blocco
        @inbounds Threads.@threads for m = m_block_start:m_block_end
            # Accediamo ai dati pre-calcolati per l'attuale 'm'
            punti_oss_m_svec = punti_oss_m_svec_arr[m]
            vers_punti_oss_m_svec = vers_punti_oss_m_svec_arr[m]
            is_pox = is_pox_arr[m]
            is_poy = is_poy_arr[m]

            # Loop sui volumi (N)
            for n = 1:N
                # Accediamo ai dati pre-calcolati per l'attuale 'n'
                barra_n_svec = barra_n_svec_arr[n]
                S_n = S_n_arr[n]

                scelta_val = :vuoto # Using Symbol for `scelta` for faster comparison

                # La logica di `scelta_val` rimane qui, in quanto dipende da `m` (is_pox, is_poy) e `n` (vers_celle[n, :])
                if abs(vers_celle[n, 1]) > 1e-5
                    if is_pox
                        scelta_val = :xx
                    elseif is_poy
                        scelta_val = :yx
                    else
                        scelta_val = :zx
                    end
                elseif abs(vers_celle[n, 2]) > 1e-5
                    if is_pox
                        scelta_val = :xy
                    elseif is_poy
                        scelta_val = :yy
                    else
                        scelta_val = :zy
                    end
                else # Assumiamo che abs(vers_celle[n, 3]) > 1e-5
                    if is_pox
                        scelta_val = :xz
                    elseif is_poy
                        scelta_val = :yz
                    else
                        scelta_val = :zz
                    end
                end

                # Chiama la funzione ottimizzata
                Lambda[m, n] = compute_hi_optimized(barra_n_svec, punti_oss_m_svec, scelta_val, rootk, wek, beta) / S_n
            end
        end
        # MATLAB had sleep and print; comment out for benchmarking, uncomment for progress
        # sleep(0.01)
        @printf "Processing block %d of %d (M: %d-%d)\n" round(Int, m_block_end/block_size_M) round(Int, M/block_size_M) m_block_start m_block_end
        if id !== nothing && is_stop_requested(id) # Check for stop request
           println("Simulation $(id) interrupted by stop request.")
           return nothing
        end
    end

    return Lambda
end


"""
    compute_hi_optimized(barra, centro_oss, scelta, rootk, wek, beta)

Optimized helper function to compute the integral for a single volume and observation point.
Uses StaticArrays for all internal vector operations to maximize performance.
"""
@inline function compute_hi_optimized(barra::SVector{24, Float64}, centro_oss::SVector{3, Float64},
                                      scelta::Symbol, rootk::AbstractVector{Float64}, wek::AbstractVector{Float64},
                                      beta::ComplexF64)

    # Pre-calculate constant inverse
    inv8 = 0.125

    # Calculate rmi, rai, etc. components directly as SVectors
    # This avoids intermediate arrays and leverages StaticArrays' efficiency
    rmi = SVector(
        inv8 * (barra[1] + barra[4] + barra[7] + barra[10] + barra[13] + barra[16] + barra[19] + barra[22]),
        inv8 * (barra[2] + barra[5] + barra[8] + barra[11] + barra[14] + barra[17] + barra[20] + barra[23]),
        inv8 * (barra[3] + barra[6] + barra[9] + barra[12] + barra[15] + barra[18] + barra[21] + barra[24])
    )

    rai = SVector(
        inv8 * (-barra[1] + barra[4] + barra[10] - barra[7] - barra[13] + barra[16] + barra[22] - barra[19]),
        inv8 * (-barra[2] + barra[5] + barra[11] - barra[8] - barra[14] + barra[17] + barra[23] - barra[20]),
        inv8 * (-barra[3] + barra[6] + barra[12] - barra[9] - barra[15] + barra[18] + barra[24] - barra[21]) # Corrected index for barra[21]
    )

    rbi = SVector(
        inv8 * (-barra[1] - barra[4] + barra[10] + barra[7] - barra[13] - barra[16] + barra[22] + barra[19]),
        inv8 * (-barra[2] - barra[5] + barra[11] + barra[8] - barra[14] - barra[17] + barra[23] + barra[20]),
        inv8 * (-barra[3] - barra[6] + barra[12] + barra[9] - barra[15] - barra[18] + barra[24] + barra[21])
    )

    rci = SVector(
        inv8 * (-barra[1] - barra[4] - barra[10] - barra[7] + barra[13] + barra[16] + barra[22] + barra[19]),
        inv8 * (-barra[2] - barra[5] - barra[11] - barra[8] + barra[14] + barra[17] + barra[23] + barra[20]),
        inv8 * (-barra[3] - barra[6] - barra[12] - barra[9] + barra[15] + barra[18] + barra[24] + barra[21])
    )

    rabi = SVector(
        inv8 * (barra[1] - barra[4] + barra[10] - barra[7] + barra[13] - barra[16] + barra[22] - barra[19]),
        inv8 * (barra[2] - barra[5] + barra[11] - barra[8] + barra[14] - barra[17] + barra[23] - barra[20]),
        inv8 * (barra[3] - barra[6] + barra[12] - barra[9] + barra[15] - barra[18] + barra[24] - barra[21])
    )

    rbci = SVector(
        inv8 * (barra[1] + barra[4] - barra[10] - barra[7] - barra[13] - barra[16] + barra[22] + barra[19]),
        inv8 * (barra[2] + barra[5] - barra[11] - barra[8] - barra[14] - barra[17] + barra[23] + barra[20]),
        inv8 * (barra[3] + barra[6] - barra[12] - barra[9] - barra[15] - barra[18] + barra[24] + barra[21])
    )

    raci = SVector(
        inv8 * (barra[1] - barra[4] - barra[10] + barra[7] - barra[13] + barra[16] + barra[22] - barra[19]),
        inv8 * (barra[2] - barra[5] - barra[11] + barra[8] - barra[14] + barra[17] + barra[23] - barra[20]),
        inv8 * (barra[3] - barra[6] - barra[12] + barra[9] - barra[15] + barra[18] + barra[24] - barra[21])
    )

    rabci = SVector(
        inv8 * (-barra[1] + barra[4] - barra[10] + barra[7] + barra[13] - barra[16] + barra[22] - barra[19]),
        inv8 * (-barra[2] + barra[5] - barra[11] + barra[8] + barra[14] - barra[17] + barra[23] - barra[20]),
        inv8 * (-barra[3] + barra[6] - barra[12] + barra[9] + barra[15] - barra[18] + barra[24] - barra[21])
    )

    nlk = length(wek) # Assuming rootkx, wekx, rootky, weky, rootkz, wekz are all the same 'wek'
    sum_a1 = zero(ComplexF64)

    # Destructure centro_oss SVector for direct access
    xo, yo, zo = centro_oss

    @fastmath @inbounds for a1 in 1:nlk
        sum_b1 = zero(ComplexF64)
        k_a1 = rootk[a1]

        for b1 in 1:nlk
            sum_c1 = zero(ComplexF64)
            k_b1 = rootk[b1]

            for c1 in 1:nlk
                k_c1 = rootk[c1]

                # Calculate components of drai, drbi, drci directly as SVectors
                drai = rai + rabi * k_b1 + raci * k_c1 + rabci * k_b1 * k_c1
                drbi = rbi + rabi * k_a1 + rbci * k_c1 + rabci * k_a1 * k_c1
                drci = rci + raci * k_a1 + rbci * k_b1 + rabci * k_a1 * k_b1

                # Use LinearAlgebra.norm() for Euclidean norm (norm(vec, 2))
                draim = norm(drai)
                drbim = norm(drbi)
                drcim = norm(drci)

                # Calculate r1 as an SVector
                r1 = rmi + rai * k_a1 + rbi * k_b1 + rci * k_c1 +
                     rabi * k_a1 * k_b1 + raci * k_a1 * k_c1 + rbci * k_b1 * k_c1 +
                     rabci * k_a1 * k_b1 * k_c1

                # Access components of r1 directly
                delta_x = (xo - r1[1])
                delta_y = (yo - r1[2])
                delta_z = (zo - r1[3])

                R_squared = delta_x*delta_x + delta_y*delta_y + delta_z*delta_z
                R = sqrt(R_squared + 1e-12) # Add a tiny epsilon to prevent division by zero

                # Pre-calculate common terms for G
                invR = 1.0 / R
                invR2 = invR * invR
                invR3 = invR2 * invR

                exp_term = exp(-1im * beta * R)
                coeff_term = (invR3 + 1im * beta * invR2)

                G_val = zero(ComplexF64) # Initialize for type stability
                if scelta == :xx
                    G_val = delta_x * exp_term * coeff_term
                elseif scelta == :xy
                    G_val = delta_z * exp_term * coeff_term
                elseif scelta == :xz
                    G_val = -delta_y * exp_term * coeff_term
                elseif scelta == :yx
                    G_val = -delta_z * exp_term * coeff_term
                elseif scelta == :yy
                    G_val = delta_y * exp_term * coeff_term
                elseif scelta == :yz
                    G_val = delta_x * exp_term * coeff_term
                elseif scelta == :zx
                    G_val = delta_y * exp_term * coeff_term
                elseif scelta == :zy
                    G_val = -delta_x * exp_term * coeff_term
                elseif scelta == :zz # This case was missing in MATLAB, implicitly handled by the final else
                    G_val = delta_z * exp_term * coeff_term
                end
                # MATLAB had an 'else' with G=0 if none matched. If 'zz' is the default and expected
                # when other conditions fail, ensure it's covered.

                f = (draim * drbim) * (drcim * G_val)
                sum_c1 += wek[c1] * f
            end # (c1)
            sum_b1 += wek[b1] * sum_c1
        end  # (b1)
        sum_a1 += wek[a1] * sum_b1
    end # (a1)

    integral = 1e-7 * sum_a1 # Apply 1e-7 factor once at the end
    return integral
end
