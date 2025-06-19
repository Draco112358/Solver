using LinearAlgebra
# using StaticArrays # Essenziale per prestazioni elevate con array di piccole dimensioni
using Base.Threads
using Printf
# using StaticArrays

# function compute_Ar_Gauss(barre::Transpose{Float64, Matrix{Float64}},
#                                     centriOss::Matrix{Float64}, ordine::Int, beta::ComplexF64,
#                                     simulation_id=nothing, chan=nothing)

#     numCentri = size(centriOss, 1)
#     numBarre = size(barre, 1)

#     ha = zeros(ComplexF64, numBarre, numCentri)

#     # --- Pre-calcoli globali (una volta per la funzione) ---

#     # 1. qrule values (già fuori, ottimo)
#     rootkx, wekx = qrule(ordine)
#     rootky, weky = qrule(ordine)
#     rootkz, wekz = qrule(ordine)
#     nlkx = length(wekx)
#     nlky = length(weky)
#     nlkz = length(wekz)

#     # 2. Pre-calcola le coordinate x_o, y_o, z_o per tutti i centri di osservazione
#     # Convertiamo direttamente in SVector{3,Float64} per comodità
#     centriOss_svec = Vector{SVector{3,Float64}}(undef, numCentri)
#     @inbounds Threads.@threads for cc = 1:numCentri
#         centriOss_svec[cc] = SVector{3,Float64}(centriOss[cc, 1], centriOss[cc, 2], centriOss[cc, 3])
#     end

#     # Questi array conterranno i 3 componenti (x,y,z) di ogni 'r' vettore per ciascuna barra
#     rmi_arr = Vector{SVector{3,Float64}}(undef, numBarre)
#     rai_arr = Vector{SVector{3,Float64}}(undef, numBarre)
#     rbi_arr = Vector{SVector{3,Float64}}(undef, numBarre)
#     rci_arr = Vector{SVector{3,Float64}}(undef, numBarre)
#     rabi_arr = Vector{SVector{3,Float64}}(undef, numBarre)
#     rbci_arr = Vector{SVector{3,Float64}}(undef, numBarre)
#     raci_arr = Vector{SVector{3,Float64}}(undef, numBarre)
#     rabci_arr = Vector{SVector{3,Float64}}(undef, numBarre)

#     @inbounds Threads.@threads for cont = 1:numBarre
#         barra_svec_24 = SVector{24, Float64}(@view barre[cont, :])
#         # Calcola e memorizza tutti gli 8 vettori 'r' (rmi, rai, ..., rabci)
#         # Re-inseriamo i calcoli completi basati sulle espressioni originali.
#         # Definiamo i vertici in modo strutturato per chiarezza.
#         # r1, r2, ..., r8 corrispondono a barra_svec_24 elementi.
#         # Vertici: (x,y,z)
#         v1 = SVector{3,Float64}(barra_svec_24[1], barra_svec_24[2], barra_svec_24[3])
#         v2 = SVector{3,Float64}(barra_svec_24[4], barra_svec_24[5], barra_svec_24[6])
#         v3 = SVector{3,Float64}(barra_svec_24[7], barra_svec_24[8], barra_svec_24[9])
#         v4 = SVector{3,Float64}(barra_svec_24[10], barra_svec_24[11], barra_svec_24[12])
#         v5 = SVector{3,Float64}(barra_svec_24[13], barra_svec_24[14], barra_svec_24[15])
#         v6 = SVector{3,Float64}(barra_svec_24[16], barra_svec_24[17], barra_svec_24[18])
#         v7 = SVector{3,Float64}(barra_svec_24[19], barra_svec_24[20], barra_svec_24[21])
#         v8 = SVector{3,Float64}(barra_svec_24[22], barra_svec_24[23], barra_svec_24[24])

#         rmi_arr[cont] = 0.125 * (v1 + v2 + v3 + v4 + v5 + v6 + v7 + v8)
#         rai_arr[cont] = 0.125 * (-v1 + v2 - v3 + v4 - v5 + v6 - v7 + v8)
#         rbi_arr[cont] = 0.125 * (-v1 - v2 + v3 + v4 - v5 - v6 + v7 + v8)
#         rci_arr[cont] = 0.125 * (-v1 - v2 - v3 - v4 + v5 + v6 + v7 + v8)
#         rabi_arr[cont] = 0.125 * (v1 - v2 + v3 - v4 + v5 - v6 + v7 - v8)
#         rbci_arr[cont] = 0.125 * (v1 + v2 - v3 - v4 - v5 - v6 + v7 + v8)
#         raci_arr[cont] = 0.125 * (v1 - v2 - v3 + v4 - v5 + v6 + v7 - v8)
#         rabci_arr[cont] = 0.125 * (-v1 + v2 - v3 + v4 + v5 - v6 + v7 - v8)
#     end


#     # --- Loop di elaborazione parallelizzata ---
#     block_size = 200 # Puoi sperimentare con questo valore

#     @inbounds for m_block_start in 1:block_size:numBarre
#         m_end = min(m_block_start + block_size - 1, numBarre)
#         Threads.@threads for cont in m_block_start:m_end

#             # Accedi ai coefficienti 'r' pre-calcolati per la barra corrente
#             current_rmi = rmi_arr[cont]
#             current_rai = rai_arr[cont]
#             current_rbi = rbi_arr[cont]
#             current_rci = rci_arr[cont]
#             current_rabi = rabi_arr[cont]
#             current_rbci = rbci_arr[cont]
#             current_raci = raci_arr[cont]
#             current_rabci = rabci_arr[cont]

#             for cc in 1:numCentri
#                 # Accedi alle coordinate pre-calcolate del centro di osservazione (già SVector)
#                 origin_vec = centriOss_svec[cc]

#                 ha[cont, cc] = compute_ha_optimized_corrected(
#                     origin_vec, # Passiamo l'SVector direttamente
#                     current_rmi, # Passiamo l'SVector direttamente
#                     current_rai,
#                     current_rbi,
#                     current_rci,
#                     current_rabi,
#                     current_rbci,
#                     current_raci,
#                     current_rabci,
#                     beta, rootkx, wekx, rootky, weky, rootkz, wekz,
#                     nlkx, nlky, nlkz
#                 )
#             end
#         end
#         #Assicurati che @printf sia disponibile (di solito con using Printf)
#         @printf "Block Ar : %.0f / %.0f\n" round(m_end / block_size) round(numBarre / block_size)
#         if is_stop_requested(simulation_id)
#             println("Simulazione $(simulation_id) interrotta per richiesta stop.")
#             return nothing
#         end
#     end

#     return ha
# end

# @inline function compute_ha_optimized_corrected(
#     orig::SVector{3, Float64}, # xo, yo, zo
#     rmi::SVector{3, Float64},
#     rai::SVector{3, Float64},
#     rbi::SVector{3, Float64},
#     rci::SVector{3, Float64},
#     rabi::SVector{3, Float64},
#     rbci::SVector{3, Float64},
#     raci::SVector{3, Float64},
#     rabci::SVector{3, Float64},
#     beta::ComplexF64,
#     rootkx::Vector{Float64}, wekx::Vector{Float64},
#     rootky::Vector{Float64}, weky::Vector{Float64},
#     rootkz::Vector{Float64}, wekz::Vector{Float64},
#     nlkx::Int64, nlky::Int64, nlkz::Int64
# )
#     sum_a1 = zero(ComplexF64)
#     @fastmath for a1 in 1:nlkx
#         sum_b1 = zero(ComplexF64)
#         kx_a1 = rootkx[a1]

#         for b1 in 1:nlky
#             sum_c1 = zero(ComplexF64)
#             ky_b1 = rootky[b1]

#             for c1 in 1:nlkz
#                 kz_c1 = rootkz[c1]

#                 # Calcolo vettoriale per drai, drbi, drci
#                 drai = rai + rabi * ky_b1 + raci * kz_c1 + rabci * ky_b1 * kz_c1
#                 draim = norm(drai) # sqrt(drai.x^2 + drai.y^2 + drai.z^2)

#                 drbi = rbi + rabi * kx_a1 + rbci * kz_c1 + rabci * kx_a1 * kz_c1
#                 # Attenzione qui: c'era un errore nel tuo codice originale in drbi_z.
#                 # Ho assunto che dovesse seguire lo schema degli altri componenti.
#                 # drbi_z = rbi_z + rabi_z * kx_a1 + rbci_z * kz_c1 + rabci_z * kx_a1 * kz_c1
#                 drbim = norm(drbi)

#                 drci = rci + raci * kx_a1 + rbci * ky_b1 + rabci * kx_a1 * ky_b1
#                 drcim = norm(drci)

#                 # Calcolo vettoriale per r1
#                 r1 = rmi + rai * kx_a1 + rbi * ky_b1 + rci * kz_c1 +
#                      rabi * kx_a1 * ky_b1 + raci * kx_a1 * kz_c1 + rbci * ky_b1 * kz_c1 +
#                      rabci * kx_a1 * ky_b1 * kz_c1

#                 delta = orig - r1
#                 R = norm(delta)

#                 invR = 1.0 / R
#                 G = invR * exp(-1im * beta * R)

#                 f = (draim * drbim) * (drcim * G)
#                 sum_c1 += wekz[c1] * f
#             end
#             sum_b1 += weky[b1] * sum_c1
#         end
#         sum_a1 += wekx[a1] * sum_b1
#     end

#     return sum_a1
# end

# # 1-3. xo, yo, zo: Float64
# #      Coordinate di un singolo centro di osservazione.
# xo_single = 1.0
# yo_single = 2.0
# zo_single = 3.0

# # 4-9. rootkx, wekx, rootky, weky, rootkz, wekz: Vector{Float64}
# #      Uguali a quelli usati in `compute_Ar_Gauss_to_para`.
# ordine_ha = 4
# rootkx_ha, wekx_ha = qrule(ordine_ha)
# rootky_ha, weky_ha = qrule(ordine_ha)
# rootkz_ha, wekz_ha = qrule(ordine_ha)

# # 10. beta: ComplexF64
# #     Uguale a quello usato nelle funzioni precedenti.
# beta_ha = 0.5 + 1.2im

# # 11-18. rmi, rai, rbi, rci, rabi, rbci, raci, rabci: Vector{Float64}
# #       Coefficienti derivati dalle coordinate della barra.
# #       Questi vengono calcolati all'interno di `compute_Ar_Gauss_to_para`.
# #       Per un input indipendente, dobbiamo simularli. Sono vettori di 3 elementi.
# rmi_ha = [0.1, 0.2, 0.3]
# rai_ha = [0.05, -0.02, 0.08]
# rbi_ha = [-0.01, 0.03, -0.06]
# rci_ha = [0.07, -0.04, 0.01]
# rabi_ha = [0.02, 0.01, -0.03]
# rbci_ha = [-0.04, 0.05, 0.02]
# raci_ha = [0.03, -0.02, 0.04]
# rabci_ha = [-0.01, 0.03, -0.05]

# # Chiamata della funzione
# # ha_single_result = compute_ha(xo_single, yo_single, zo_single,
# #                               rootkx_ha, wekx_ha, rootky_ha, weky_ha, rootkz_ha, wekz_ha, beta_ha,
# #                               rmi_ha, rai_ha, rbi_ha, rci_ha, rabi_ha, rbci_ha, raci_ha, rabci_ha)
# # println("\nRisultato di compute_ha: ", ha_single_result)

function compute_Ar_Gauss(barre::Transpose{Float64, Matrix{Float64}},
                                    centriOss::Matrix{Float64}, ordine::Int, beta::ComplexF64,
                                    simulation_id=nothing, chan=nothing)
    numCentri = size(centriOss, 1)
    N = size(barre, 1)

    num_punti = ceil(Int, N / 100)
    vect_ind = ceil.(Int, LinRange(1, N, num_punti))

    if length(vect_ind) < 2
        vect_ind = [1, N]
    end
    if vect_ind[end] > N
        vect_ind[end] = N
    end

    Nq = length(vect_ind) - 1
    barrec = Vector{Matrix{Float64}}(undef, Nq) 

    for pos = 1:Nq
        if pos == Nq
            barrec[pos] = barre[vect_ind[pos]:vect_ind[pos+1], :]
        else
            barrec[pos] = barre[vect_ind[pos]:vect_ind[pos+1]-1, :]
        end
    end

    ha_local = Vector{Matrix{ComplexF64}}(undef, Nq)

    rootkx, wekx = qrule(ordine)
    rootky, weky = qrule(ordine)
    rootkz, wekz = qrule(ordine)

    Threads.@threads for c = 1:Nq
        ha_local[c] = compute_Ar_Gauss_to_para(barrec[c], centriOss, beta, rootkx, wekx, rootky, weky, rootkz, wekz)
        println("block Ar : ", c, " / ", Nq) # Commentato per evitare overhead di I/O nel profiler
    end

    ha = zeros(ComplexF64, N, numCentri)

    for c = 1:Nq
        if c == Nq
            ha[vect_ind[c]:vect_ind[c+1], :] = ha_local[c]
        else
            ha[vect_ind[c]:vect_ind[c+1]-1, :] = ha_local[c]
        end
    end

    return ha
end


function compute_Ar_Gauss_to_para(barre::Matrix{Float64}, centriOss::Matrix{Float64}, beta::ComplexF64, rootkx::Vector{Float64}, wekx::Vector{Float64}, rootky::Vector{Float64}, weky::Vector{Float64}, rootkz::Vector{Float64}, wekz::Vector{Float64})
    numCentri = size(centriOss, 1)
    
    ha = zeros(ComplexF64, size(barre, 1), numCentri)
    
    # Pre-allocare `ri` e i vettori dei coefficienti qui per ridurre le allocazioni per barra
    # (anche se il profiler non li ha indicati come hotspot significativi, è buona pratica)
    ri = Matrix{Float64}(undef, 8, 3) 
    rmi = Vector{Float64}(undef, 3)
    rai = Vector{Float64}(undef, 3)
    rbi = Vector{Float64}(undef, 3)
    rci = Vector{Float64}(undef, 3)
    rabi = Vector{Float64}(undef, 3)
    rbci = Vector{Float64}(undef, 3)
    raci = Vector{Float64}(undef, 3)
    rabci = Vector{Float64}(undef, 3)

    for cont = 1:size(barre, 1)
        barra_row = barre[cont, :] # Rinominato per evitare conflitto con la variabile 'barra' interna
        
        xb = [barra_row[1], barra_row[4], barra_row[7], barra_row[10], barra_row[13], barra_row[16], barra_row[19], barra_row[22]]
        yb = [barra_row[2], barra_row[5], barra_row[8], barra_row[11], barra_row[14], barra_row[17], barra_row[20], barra_row[23]]
        zb = [barra_row[3], barra_row[6], barra_row[9], barra_row[12], barra_row[15], barra_row[18], barra_row[21], barra_row[24]]
        
        x_bar = [minimum(xb), maximum(xb)]
        y_bar = [minimum(yb), maximum(yb)]
        z_bar = [minimum(zb), maximum(zb)]

        x1, x2 = x_bar[1], x_bar[2]
        y1, y2 = y_bar[1], y_bar[2]
        z1, z2 = z_bar[1], z_bar[2]
        
        # Define 'barra' (vertices of the 3D box)
        # Nota: questa assegnazione sovrascrive la variabile barra_row. Potrebbe essere meglio rinominare.
        # Ho rinominato la variabile di input a `barra_row` per evitarlo.
        temp_barra_coords = [x1, y1, z1, x2, y1, z1, x1, y2, z1, x2, y2, z1, x1, y1, z2, x2, y1, z2, x1, y2, z2, x2, y2, z2]
        
        xi1 = [temp_barra_coords[1], temp_barra_coords[4], temp_barra_coords[7], temp_barra_coords[10]]
        yi1 = [temp_barra_coords[2], temp_barra_coords[5], temp_barra_coords[8], temp_barra_coords[11]]
        zi1 = [temp_barra_coords[3], temp_barra_coords[6], temp_barra_coords[9], temp_barra_coords[12]]
        xi2 = [temp_barra_coords[13], temp_barra_coords[16], temp_barra_coords[19], temp_barra_coords[22]]
        yi2 = [temp_barra_coords[14], temp_barra_coords[17], temp_barra_coords[20], temp_barra_coords[23]]
        zi2 = [temp_barra_coords[15], temp_barra_coords[18], temp_barra_coords[21], temp_barra_coords[24]]
        
        # Assegnazioni in-place a `ri`
        ri[1, :] = [xi1[1], yi1[1], zi1[1]]
        ri[2, :] = [xi1[2], yi1[2], zi1[2]]
        ri[3, :] = [xi1[3], yi1[3], zi1[3]]
        ri[4, :] = [xi1[4], yi1[4], zi1[4]]
        ri[5, :] = [xi2[1], yi2[1], zi2[1]]
        ri[6, :] = [xi2[2], yi2[2], zi2[2]]
        ri[7, :] = [xi2[3], yi2[3], zi2[3]]
        ri[8, :] = [xi2[4], yi2[4], zi2[4]]
        
        # Calcolo dei coefficienti in-place
        @inbounds for i in 1:3
            rmi[i] = 0.125 * (ri[1,i] + ri[2,i] + ri[3,i] + ri[4,i] + ri[5,i] + ri[6,i] + ri[7,i] + ri[8,i])
            rai[i] = 0.125 * (-ri[1,i] + ri[2,i] + ri[4,i] - ri[3,i] - ri[5,i] + ri[6,i] + ri[8,i] - ri[7,i])
            rbi[i] = 0.125 * (-ri[1,i] - ri[2,i] + ri[4,i] + ri[3,i] - ri[5,i] - ri[6,i] + ri[8,i] + ri[7,i])
            rci[i] = 0.125 * (-ri[1,i] - ri[2,i] - ri[4,i] - ri[3,i] + ri[5,i] + ri[6,i] + ri[8,i] + ri[7,i])
            rabi[i] = 0.125 * (ri[1,i] - ri[2,i] + ri[4,i] - ri[3,i] + ri[5,i] - ri[6,i] + ri[8,i] - ri[7,i])
            rbci[i] = 0.125 * (ri[1,i] + ri[2,i] - ri[4,i] - ri[3,i] - ri[5,i] - ri[6,i] + ri[8,i] + ri[7,i])
            raci[i] = 0.125 * (ri[1,i] - ri[2,i] - ri[4,i] + ri[3,i] - ri[5,i] + ri[6,i] + ri[8,i] - ri[7,i])
            rabci[i] = 0.125 * (-ri[1,i] + ri[2,i] - ri[4,i] + ri[3,i] + ri[5,i] - ri[6,i] + ri[8,i] - ri[7,i])
        end
        
        for cc = 1:numCentri
            x_o = centriOss[cc, 1]
            y_o = centriOss[cc, 2]
            z_o = centriOss[cc, 3]
            
            ha[cont, cc] = compute_ha(x_o, y_o, z_o, rootkx, wekx, rootky, weky, rootkz, wekz, beta, rmi, rai, rbi, rci, rabi, rbci, raci, rabci)
        end
    end
    
    return ha
end


function compute_ha(xo::Float64, yo::Float64, zo::Float64,
                    rootkx::Vector{Float64}, wekx::Vector{Float64}, rootky::Vector{Float64}, weky::Vector{Float64}, rootkz::Vector{Float64}, wekz::Vector{Float64}, beta::ComplexF64,
                    rmi::Vector{Float64}, rai::Vector{Float64}, rbi::Vector{Float64}, rci::Vector{Float64},
                    rabi::Vector{Float64}, rbci::Vector{Float64}, raci::Vector{Float64}, rabci::Vector{Float64})
    
    nlkx = length(wekx)
    nlky = length(weky)
    nlkz = length(wekz)
    
    sum_a1 = 0.0 + im*0.0
    
    _drai_temp = Vector{Float64}(undef, 3)
    _drbi_temp = Vector{Float64}(undef, 3)
    _drci_temp = Vector{Float64}(undef, 3)
    _r1_temp = Vector{Float64}(undef, 3)
    
    @inbounds for a1 = 1:nlkx
        sum_b1 = 0.0 + im*0.0
        @inbounds for b1 = 1:nlky
            sum_c1 = 0.0 + im*0.0
             @inbounds for c1 = 1:nlkz
                dkx_val = rootkx[a1]
                dky_val = rootky[b1]
                dkz_val = rootkz[c1]

                # *** MODIFICA IMPORTANTE: Sostituzione di .= con assegnazioni elemento per elemento ***
                _drai_temp[1] = rai[1] + rabi[1] * dky_val + raci[1] * dkz_val + rabci[1] * dky_val * dkz_val
                _drai_temp[2] = rai[2] + rabi[2] * dky_val + raci[2] * dkz_val + rabci[2] * dky_val * dkz_val
                _drai_temp[3] = rai[3] + rabi[3] * dky_val + raci[3] * dkz_val + rabci[3] * dky_val * dkz_val

                _drbi_temp[1] = rbi[1] + rabi[1] * dkx_val + rbci[1] * dkz_val + rabci[1] * dkx_val * dkz_val
                _drbi_temp[2] = rbi[2] + rabi[2] * dkx_val + rbci[2] * dkz_val + rabci[2] * dkx_val * dkz_val
                _drbi_temp[3] = rbi[3] + rabi[3] * dkx_val + rbci[3] * dkz_val + rabci[3] * dkx_val * dkz_val

                _drci_temp[1] = rci[1] + raci[1] * dkx_val + rbci[1] * dky_val + rabci[1] * dkx_val * dky_val
                _drci_temp[2] = rci[2] + raci[2] * dkx_val + rbci[2] * dky_val + rabci[2] * dkx_val * dky_val
                _drci_temp[3] = rci[3] + raci[3] * dkx_val + rbci[3] * dky_val + rabci[3] * dkx_val * dky_val
                
                # Calcolo manuale della norma L2 (già suggerito e implementato, se non lo era)
                draim = sqrt(_drai_temp[1]^2 + _drai_temp[2]^2 + _drai_temp[3]^2)
                drbim = sqrt(_drbi_temp[1]^2 + _drbi_temp[2]^2 + _drbi_temp[3]^2)
                drcim = sqrt(_drci_temp[1]^2 + _drci_temp[2]^2 + _drci_temp[3]^2)
                
                # *** MODIFICA IMPORTANTE: Sostituzione di .= con assegnazioni elemento per elemento ***
                _r1_temp[1] = rmi[1] + rai[1] * dkx_val + rbi[1] * dky_val + rci[1] * dkz_val + 
                              rabi[1] * dkx_val * dky_val + raci[1] * dkx_val * dkz_val + 
                              rbci[1] * dky_val * dkz_val + rabci[1] * dkx_val * dky_val * dkz_val
                _r1_temp[2] = rmi[2] + rai[2] * dkx_val + rbi[2] * dky_val + rci[2] * dkz_val + 
                              rabi[2] * dkx_val * dky_val + raci[2] * dkx_val * dkz_val + 
                              rbci[2] * dky_val * dkz_val + rabci[2] * dkx_val * dky_val * dkz_val
                _r1_temp[3] = rmi[3] + rai[3] * dkx_val + rbi[3] * dky_val + rci[3] * dkz_val + 
                              rabi[3] * dkx_val * dky_val + raci[3] * dkx_val * dkz_val + 
                              rbci[3] * dky_val * dkz_val + rabci[3] * dkx_val * dky_val * dkz_val
                
                x = _r1_temp[1]
                y = _r1_temp[2]
                z = _r1_temp[3]
                
                delta_x = (xo - x)
                delta_y = (yo - y)
                delta_z = (zo - z)
                
                R = sqrt(delta_x^2 + delta_y^2 + delta_z^2)
                
                G = ifelse(R > eps(Float64), (1/R) * exp(-im * beta * R), ComplexF64(0.0))
                
                f = draim * drbim * drcim * G # Questa riga (o la riga 77 con l'accumulatore) potrebbe essere il prossimo hotspot.
                
                sum_c1 += wekz[c1] * f
            end
            sum_b1 += weky[b1] * sum_c1
        end
        sum_a1 += wekx[a1] * sum_b1 # Questa riga è REPL[8] line 77, con overhead 66
    end
    
    return sum_a1
end


function qrule(n::Int)
    iter = 2
    m = trunc((n + 1) / 2)
    e1 = n * (n + 1)
    mm = 4 * m - 1
    t = (pi / (4 * n + 2)) * (3:4:mm)
    nn = (1 - (1 - 1 / n) / (8 * n * n))
    xo = nn * cos.(t)
    den = []
    d1 = []
    dpn = []
    d2pn = []
    d3pn = []
    d4pn = []
    u = []
    v = []
    h = []
    p = []
    dp = []
    pk = []

    for kk ∈ 1:iter
        pkm1 = zeros(size(xo))
        pkm1[1:size(xo, 1)] .= 1
        pk = xo
        for k ∈ 2:n
            t1 = xo .* pk
            pkp1 = t1 - pkm1 - (t1 - pkm1) / k + t1
            pkm1 = pk
            pk = pkp1
        end
        den = 1 .- xo .^ 2
        d1 = n * (pkm1 - xo .* pk)
        dpn = d1 ./ den
        d2pn = (2 .* xo .* dpn - e1 .* pk) ./ den
        d3pn = (4 * xo .* d2pn + (2 - e1) .* dpn) ./ den
        d4pn = (6 * xo .* d3pn + (6 - e1) .* d2pn) ./ den
        u = pk ./ dpn
        v = d2pn ./ dpn
        h = -u .* (1 .+ (0.5 * u) .* (v + u .* (v .* v - u .* d3pn ./ (3 * dpn))))
        p = pk + h .* (dpn + (0.5 * h) .* (d2pn + (h / 3) .* (d3pn + 0.25 * h .* d4pn)))
        dp = dpn + h .* (d2pn + (0.5 * h) .* (d3pn + h .* d4pn / 3))
        h = h - p ./ dp
        xo = xo + h
    end
    bp = zeros(1, n)
    wf = zeros(1, n)
    bp[1:size(xo, 1)] .= -xo .- h
    fx = d1 - h .* e1 .* (pk + (h / 2) .* (dpn + (h / 3) .* (
        d2pn + (h / 4) .* (d3pn + (0.2 * h) .* d4pn))))
    wf[1:size(xo, 1)] .= 2 * (1 .- bp[1:size(xo, 1)] .^ 2) ./ (fx .* fx)
    if (m + m) > n
        bp[Int64(m)] = 0
    end
    if !((m + m) == n)
        m = m - 1
    end
    jj = 1:m
    n1j = (n + 1) .- jj
    bp[Int64.(n1j)] .= -bp[Int64.(jj)]
    wf[Int64.(n1j)] .= wf[Int64.(jj)]

    return vec(bp), vec(wf)
end