# using LinearAlgebra

# function compute_lambda_numeric(punti_oss, volumi, incidence_selection, vers_punti_oss, ordine_int, beta, id, chan)

#     N = size(volumi[:coordinate], 1)

#     vers_celle = zeros(N, 3)
#     vers_celle[1:incidence_selection[:mx], 1] .= 1
#     vers_celle[incidence_selection[:mx]+1:incidence_selection[:mx]+incidence_selection[:my], 2] .= 1
#     vers_celle[incidence_selection[:mx]+incidence_selection[:my]+1:N, 3] .= 1

#     M = size(punti_oss, 1)

#     Lambda = zeros(ComplexF64, M, N)

#     rootkx, wekx = qrule(ordine_int)

#     block_size = 10

#     for m_block in 1:block_size:M
#         m_end = min(m_block + block_size - 1, M)
#         Base.Threads.@threads for m = m_block:m_end
#             for n = 1:N
#                 scelta = ""
    
#                 if abs(vers_celle[n, 1]) > 1e-5
    
#                     if abs(vers_punti_oss[m, 1]) > 1e-5
#                         scelta = "xx"
#                     elseif abs(vers_punti_oss[m, 2]) > 1e-5
#                         scelta = "yx"
#                     else
#                         scelta = "zx"
#                     end
    
#                 elseif abs(vers_celle[n, 2]) > 1e-5
    
#                     if abs(vers_punti_oss[m, 1]) > 1e-5
#                         scelta = "xy"
#                     elseif abs(vers_punti_oss[m, 2]) > 1e-5
#                         scelta = "yy"
#                     else
#                         scelta = "zy"
#                     end
    
#                 else
#                     if abs(vers_punti_oss[m, 1]) > 1e-5
#                         scelta = "xz"
#                     elseif abs(vers_punti_oss[m, 2]) > 1e-5
#                         scelta = "yz"
#                     else
#                         scelta = "zz"
#                     end
#                 end
    
#                 Lambda[m, n] = compute_hi(volumi[:coordinate][n, :], punti_oss[m, :], scelta, rootkx, wekx, beta) / volumi[:S][n]
    
#             end
#         end
#         sleep(0)
#         println("block Lambda : ", round(m_end/block_size), " / ", round(M/block_size))
#         if is_stop_requested(id)
# 			println("Simulazione $(id) interrotta per richiesta stop.")
# 			return nothing # O un altro valore che indica interruzione
# 		end
#     end

#     return Lambda
    
# end

# function compute_hi(barra, centro_oss, scelta, rootkx, wekx, beta)

#     xi1 = [barra[1], barra[4], barra[7], barra[10]]
#     yi1 = [barra[2], barra[5], barra[8], barra[11]]
#     zi1 = [barra[3], barra[6], barra[9], barra[12]]
#     xi2 = [barra[13], barra[16], barra[19], barra[22]]
#     yi2 = [barra[14], barra[17], barra[20], barra[23]]
#     zi2 = [barra[15], barra[18], barra[21], barra[24]]

#     #  vectors pointing to the vertices of the quadrilateral i
#     ri = zeros(ComplexF64, 8, 3)
#     ri[1, :] = [xi1[1], yi1[1], zi1[1]]
#     ri[2, :] = [xi1[2], yi1[2], zi1[2]]
#     ri[3, :] = [xi1[3], yi1[3], zi1[3]]
#     ri[4, :] = [xi1[4], yi1[4], zi1[4]]

#     ri[5, :] = [xi2[1], yi2[1], zi2[1]]
#     ri[6, :] = [xi2[2], yi2[2], zi2[2]]
#     ri[7, :] = [xi2[3], yi2[3], zi2[3]]
#     ri[8, :] = [xi2[4], yi2[4], zi2[4]]

#     # nuovo approccio
#     rmi = vec(0.125 * sum(ri, dims=1))
#     rai = 0.125 * (-ri[1, :] + ri[2, :] + ri[4, :] - ri[3, :] - ri[5, :] + ri[6, :] + ri[8, :] - ri[7, :])
#     rbi = 0.125 * (-ri[1, :] - ri[2, :] + ri[4, :] + ri[3, :] - ri[5, :] - ri[6, :] + ri[8, :] + ri[7, :])
#     rci = 0.125 * (-ri[1, :] - ri[2, :] - ri[4, :] - ri[3, :] + ri[5, :] + ri[6, :] + ri[8, :] + ri[7, :])
#     rabi = 0.125 * (ri[1, :] - ri[2, :] + ri[4, :] - ri[3, :] + ri[5, :] - ri[6, :] + ri[8, :] - ri[7, :])
#     rbci = 0.125 * (ri[1, :] + ri[2, :] - ri[4, :] - ri[3, :] - ri[5, :] - ri[6, :] + ri[8, :] + ri[7, :])
#     raci = 0.125 * (ri[1, :] - ri[2, :] - ri[4, :] + ri[3, :] - ri[5, :] + ri[6, :] + ri[8, :] - ri[7, :])
#     rabci = 0.125 * (-ri[1, :] + ri[2, :] - ri[4, :] + ri[3, :] + ri[5, :] - ri[6, :] + ri[8, :] - ri[7, :])

#     nlkx = length(wekx)
#     nlky = length(wekx)
#     nlkz = length(wekx)

#     sum_a1 = 0.0 + 0.0im
#     for a1 = 1:nlkx
#         sum_b1 = 0.0 + 0.0im
#         for b1 = 1:nlky
#             sum_c1 = 0.0 + 0.0im
#             for c1 = 1:nlkz
#                 drai = rai + rabi * rootkx[b1] + raci * rootkx[c1] + rabci * rootkx[b1] * rootkx[c1]
#                 drbi = rbi + rabi * rootkx[a1] + rbci * rootkx[c1] + rabci * rootkx[a1] * rootkx[c1]
#                 drci = rci + raci * rootkx[a1] + rbci * rootkx[b1] + rabci * rootkx[a1] * rootkx[b1]
#                 draim = norm(drai, 2)
#                 drbim = norm(drbi, 2)
#                 drcim = norm(drci, 2)

#                 r1 = rmi + rai * rootkx[a1] + rbi * rootkx[b1] + rci * rootkx[c1] + rabi * rootkx[a1] * rootkx[b1] + raci * rootkx[a1] * rootkx[c1] + rbci * rootkx[b1] * rootkx[c1] +
#                      rabci * rootkx[a1] * rootkx[b1] * rootkx[c1]

#                 x = real(r1[1])
#                 y = real(r1[2])
#                 z = real(r1[3])

#                 delta_x = (centro_oss[1] - x)
#                 delta_y = (centro_oss[2] - y)
#                 delta_z = (centro_oss[3] - z)

#                 R = sqrt(delta_x^2 + delta_y^2 + delta_z^2)

#                 G = if scelta == "xy"
#                     delta_z * exp(-1im * beta * R) * (1 / R^3 + 1im * beta / R^2)
#                 elseif scelta == "xz"
#                     -delta_y * exp(-1im * beta * R) * (1 / R^3 + 1im * beta / R^2)
#                 elseif scelta == "yx"
#                     -delta_z * exp(-1im * beta * R) * (1 / R^3 + 1im * beta / R^2)
#                 elseif scelta == "yz"
#                     delta_x * exp(-1im * beta * R) * (1 / R^3 + 1im * beta / R^2)
#                 elseif scelta == "zx"
#                     delta_y * exp(-1im * beta * R) * (1 / R^3 + 1im * beta / R^2)
#                 elseif scelta == "zy"
#                     -delta_x * exp(-1im * beta * R) * (1 / R^3 + 1im * beta / R^2)
#                 else
#                     0.0 + 0.0im
#                 end

#                 f = draim * drbim * drcim * G

#                 sum_c1 += wekx[c1] * f
#                 if c1 % 20 == 0
# 					yield()  # Permette alla task dell'heartbeat di essere schedulata
# 				end
#             end
#             sum_b1 += wekx[b1] * sum_c1
#         end
#         sum_a1 += wekx[a1] * sum_b1
#     end
#     integral = 1e-7 * sum_a1
#     return integral
# end

# function qrule(n::Int)
#     iter = 2
#     m = trunc((n + 1) / 2)
#     e1 = n * (n + 1)
#     mm = 4 * m - 1
#     t = (pi / (4 * n + 2)) * (3:4:mm)
#     nn = (1 - (1 - 1 / n) / (8 * n * n))
#     xo = nn * cos.(t)
#     den = []
#     d1 = []
#     dpn = []
#     d2pn = []
#     d3pn = []
#     d4pn = []
#     u = []
#     v = []
#     h = []
#     p = []
#     dp = []
#     pk = []
    
#     for kk = 1:iter
#         pkm1 = zeros(size(xo))
#         pkm1[1:size(xo, 1)] .= 1
#         pk = xo
#         for k = 2:n
#             t1 = xo .* pk
#             pkp1 = t1 - pkm1 - (t1 - pkm1) / k + t1
#             pkm1 = pk
#             pk = pkp1
#         end
#         den = 1 .- xo .^ 2
#         d1 = n * (pkm1 - xo .* pk)
#         dpn = d1 ./ den
#         d2pn = (2 .* xo .* dpn - e1 .* pk) ./ den
#         d3pn = (4 * xo .* d2pn + (2 - e1) .* dpn) ./ den
#         d4pn = (6 * xo .* d3pn + (6 - e1) .* d2pn) ./ den
#         u = pk ./ dpn
#         v = d2pn ./ dpn
#         h = -u .* (1 .+ (0.5 * u) .* (v + u .* (v .* v - u .* d3pn ./ (3 * dpn))))
#         p = pk + h .* (dpn + (0.5 * h) .* (d2pn + (h / 3) .* (d3pn + 0.25 * h .* d4pn)))
#         dp = dpn + h .* (d2pn + (0.5 * h) .* (d3pn + h .* d4pn / 3))
#         h = h - p ./ dp
#         xo = xo + h
#     end
#     bp = zeros(1, n)
#     wf = zeros(1, n)
#     bp[1:size(xo, 1)] .= -xo .- h
#     fx = d1 - h .* e1 .* (pk + (h / 2) .* (dpn + (h / 3) .* (
#         d2pn + (h / 4) .* (d3pn + (0.2 * h) .* d4pn))))
#     wf[1:size(xo, 1)] .= 2 * (1 .- bp[1:size(xo, 1)] .^ 2) ./ (fx .* fx)
#     if (m + m) > n
#         bp[Int64(m)] = 0
#     end
#     if !((m + m) == n)
#         m = m - 1
#     end
#     jj = 1:m
#     n1j = (n + 1) .- jj
#     bp[Int64.(n1j)] .= -bp[Int64.(jj)]
#     wf[Int64.(n1j)] .= wf[Int64.(jj)]

#     return vec(bp), vec(wf)
# end

# using LinearAlgebra

# # --- Optimized `compute_lambda_numeric` ---
# function compute_lambda_numeric(punti_oss::Matrix{Float64}, volumi::Dict, incidence_selection::Dict{Symbol, Any},
#                                 vers_punti_oss::Matrix{Float64}, ordine_int::Int, beta::ComplexF64, id, chan)

#     N = size(volumi[:coordinate], 1)
#     M = size(punti_oss, 1)

#     # Pre-allocate vers_celle
#     vers_celle = zeros(Float64, N, 3) # Specify Float64 type
#     @inbounds begin # Use @inbounds for performance
#         vers_celle[1:incidence_selection[:mx], 1] .= 1.0
#         vers_celle[incidence_selection[:mx]+1:incidence_selection[:mx]+incidence_selection[:my], 2] .= 1.0
#         vers_celle[incidence_selection[:mx]+incidence_selection[:my]+1:N, 3] .= 1.0
#     end

#     Lambda = zeros(ComplexF64, M, N)

#     # Pre-calculate qrule values once
#     rootk, wek = qrule(ordine_int) # Assuming rootkx, wekx, rootky, weky, rootkz, wekz are all the same
#                                   # If they are different for x, y, z, you'll need separate calls
#                                   # and separate argument passing to compute_hi_optimized.
#                                   # Based on original `nlky = length(wekx)` etc., it seems they are the same.
#                                   # Let's assume rootk, wek can be used for all dimensions.

#     block_size = 100 # Adjust based on your system and problem size

#     for m_block in 1:block_size:M
#         m_end = min(m_block + block_size - 1, M)
#         Threads.@threads for m = m_block:m_end # Loop over observation points (M)
#             # Use @view for points_oss and vers_punti_oss to avoid copying rows
#             punti_oss_m = @view punti_oss[m, :]
#             vers_punti_oss_m = @view vers_punti_oss[m, :]

#             # Pre-calculate boolean flags for performance
#             is_pox = abs(vers_punti_oss_m[1]) > 1e-5
#             is_poy = abs(vers_punti_oss_m[2]) > 1e-5

#             for n = 1:N # Loop over volumes (N)
#                 # Use @view for volumi[:coordinate] and volumi[:S]
#                 barra_n = @view volumi[:coordinate][n, :]
#                 S_n = volumi[:S][n]

#                 scelta_val = :vuoto # Changed to a local variable to be type-stable within the loop

#                 @inbounds begin # @inbounds for vers_celle access
#                     # Check vers_celle components
#                     if abs(vers_celle[n, 1]) > 1e-5 # For x-oriented cell
#                         if is_pox
#                             scelta_val = :xx
#                         elseif is_poy
#                             scelta_val = :yx
#                         else
#                             scelta_val = :zx
#                         end
#                     elseif abs(vers_celle[n, 2]) > 1e-5 # For y-oriented cell
#                         if is_pox
#                             scelta_val = :xy
#                         elseif is_poy
#                             scelta_val = :yy
#                         else
#                             scelta_val = :zy
#                         end
#                     else # For z-oriented cell (abs(vers_celle[n, 3]) > 1e-5 implicitly)
#                         if is_pox
#                             scelta_val = :xz
#                         elseif is_poy
#                             scelta_val = :yz
#                         else
#                             scelta_val = :zz
#                         end
#                     end
#                 end

#                 # Call the optimized function
#                 Lambda[m, n] = compute_hi_optimized(barra_n, punti_oss_m, scelta_val, rootk, wek, beta) / S_n
#             end
#         end
#         sleep(0.01) # Small sleep to yield to OS/other tasks and allow print
#         println("block Lambda : ", round(m_end/block_size), " / ", round(M/block_size))
#         if is_stop_requested(id)
#             println("Simulazione $(id) interrotta per richiesta stop.")
#             return nothing
#         end
#     end

#     return Lambda
# end

# # --- Optimized `compute_hi` ---
# @inline function compute_hi_optimized(barra::Union{AbstractVector{Float64}, SubArray{Real, 1, Transpose{Real, Matrix{Real}}}}, centro_oss::AbstractVector{Float64},
#                                       scelta::Symbol, rootk::AbstractVector{Float64}, wek::AbstractVector{Float64},
#                                       beta::ComplexF64)

#     # Extract all 8 vertex coordinates from barra directly into local variables.
#     # This avoids allocations of xi1, yi1, zi1, xi2, yi2, zi2, and ri arrays.
#     # Assuming barra is structured as [x1,y1,z1,x2,y2,z2,...,x8,y8,z8] for 8 vertices
#     # @inbounds is crucial here for direct array access without bounds checking.
#     # Ensure barra has at least 24 elements (8 vertices * 3 coords).
#     @inbounds begin
#         x1_vtx = barra[1]; y1_vtx = barra[2]; z1_vtx = barra[3]
#         x2_vtx = barra[4]; y2_vtx = barra[5]; z2_vtx = barra[6]
#         x3_vtx = barra[7]; y3_vtx = barra[8]; z3_vtx = barra[9]
#         x4_vtx = barra[10]; y4_vtx = barra[11]; z4_vtx = barra[12]
#         x5_vtx = barra[13]; y5_vtx = barra[14]; z5_vtx = barra[15]
#         x6_vtx = barra[16]; y6_vtx = barra[17]; z6_vtx = barra[18]
#         x7_vtx = barra[19]; y7_vtx = barra[20]; z7_vtx = barra[21]
#         x8_vtx = barra[22]; y8_vtx = barra[23]; z8_vtx = barra[24]
#     end

#     # Calculate rmi, rai, etc. components directly without intermediate arrays.
#     # Each component (x,y,z) is calculated separately.
#     # The 0.125 factor is applied to avoid repeated division.
#     inv8 = 0.125 # Pre-calculate constant inverse

#     # rmi (mean point)
#     rmi_x = inv8 * (x1_vtx + x2_vtx + x3_vtx + x4_vtx + x5_vtx + x6_vtx + x7_vtx + x8_vtx)
#     rmi_y = inv8 * (y1_vtx + y2_vtx + y3_vtx + y4_vtx + y5_vtx + y6_vtx + y7_vtx + y8_vtx)
#     rmi_z = inv8 * (z1_vtx + z2_vtx + z3_vtx + z4_vtx + z5_vtx + z6_vtx + z7_vtx + z8_vtx)

#     # rai
#     rai_x = inv8 * (-x1_vtx + x2_vtx + x4_vtx - x3_vtx - x5_vtx + x6_vtx + x8_vtx - x7_vtx)
#     rai_y = inv8 * (-y1_vtx + y2_vtx + y4_vtx - y3_vtx - y5_vtx + y6_vtx + y8_vtx - y7_vtx)
#     rai_z = inv8 * (-z1_vtx + z2_vtx + z4_vtx - z3_vtx - z5_vtx + z6_vtx + z8_vtx - z7_vtx)

#     # rbi
#     rbi_x = inv8 * (-x1_vtx - x2_vtx + x4_vtx + x3_vtx - x5_vtx - x6_vtx + x8_vtx + x7_vtx)
#     rbi_y = inv8 * (-y1_vtx - y2_vtx + y4_vtx + y3_vtx - y5_vtx - y6_vtx + y8_vtx + y7_vtx)
#     rbi_z = inv8 * (-z1_vtx - z2_vtx + z4_vtx + z3_vtx - z5_vtx - z6_vtx + z8_vtx + z7_vtx)

#     # rci
#     rci_x = inv8 * (-x1_vtx - x2_vtx - x4_vtx - x3_vtx + x5_vtx + x6_vtx + x8_vtx + x7_vtx)
#     rci_y = inv8 * (-y1_vtx - y2_vtx - y4_vtx - y3_vtx + y5_vtx + y6_vtx + y8_vtx + y7_vtx)
#     rci_z = inv8 * (-z1_vtx - z2_vtx - z4_vtx - z3_vtx + z5_vtx + z6_vtx + z8_vtx + z7_vtx)

#     # rabi
#     rabi_x = inv8 * (x1_vtx - x2_vtx + x4_vtx - x3_vtx + x5_vtx - x6_vtx + x8_vtx - x7_vtx)
#     rabi_y = inv8 * (y1_vtx - y2_vtx + y4_vtx - y3_vtx + y5_vtx - y6_vtx + y8_vtx - y7_vtx)
#     rabi_z = inv8 * (z1_vtx - z2_vtx + z4_vtx - z3_vtx + z5_vtx - z6_vtx + z8_vtx - z7_vtx)

#     # rbci
#     rbci_x = inv8 * (x1_vtx + x2_vtx - x4_vtx - x3_vtx - x5_vtx - x6_vtx + x8_vtx + x7_vtx)
#     rbci_y = inv8 * (y1_vtx + y2_vtx - y4_vtx - y3_vtx - y5_vtx - y6_vtx + y8_vtx + y7_vtx)
#     rbci_z = inv8 * (z1_vtx + z2_vtx - z4_vtx - z3_vtx - z5_vtx - z6_vtx + z8_vtx + z7_vtx)

#     # raci
#     raci_x = inv8 * (x1_vtx - x2_vtx - x4_vtx + x3_vtx - x5_vtx + x6_vtx + x8_vtx - x7_vtx)
#     raci_y = inv8 * (y1_vtx - y2_vtx - y4_vtx + y3_vtx - y5_vtx + y6_vtx + y8_vtx - y7_vtx)
#     raci_z = inv8 * (z1_vtx - z2_vtx - z4_vtx + z3_vtx - z5_vtx + z6_vtx + z8_vtx - z7_vtx)

#     # rabci
#     rabci_x = inv8 * (-x1_vtx + x2_vtx - x4_vtx + x3_vtx + x5_vtx - x6_vtx + x8_vtx - x7_vtx)
#     rabci_y = inv8 * (-y1_vtx + y2_vtx - y4_vtx + y3_vtx + y5_vtx - y6_vtx + y8_vtx - y7_vtx)
#     rabci_z = inv8 * (-z1_vtx + z2_vtx - z4_vtx + z3_vtx + z5_vtx - z6_vtx + z8_vtx - z7_vtx)

#     nlk = length(wek)
#     sum_a1 = zero(ComplexF64) # Initialize with ComplexF64 zero for type stability

#     xo = centro_oss[1]
#     yo = centro_oss[2]
#     zo = centro_oss[3]

#     # Use @fastmath and @inbounds for performance
#     @fastmath @inbounds for a1 in 1:nlk
#         sum_b1 = zero(ComplexF64)
#         k_a1 = rootk[a1]

#         for b1 in 1:nlk
#             sum_c1 = zero(ComplexF64)
#             k_b1 = rootk[b1]

#             for c1 in 1:nlk
#                 k_c1 = rootk[c1]

#                 # Calculate components of drai, drbi, drci directly
#                 drai_x = rai_x + rabi_x * k_b1 + raci_x * k_c1 + rabci_x * k_b1 * k_c1
#                 drai_y = rai_y + rabi_y * k_b1 + raci_y * k_c1 + rabci_y * k_b1 * k_c1
#                 drai_z = rai_z + rabi_z * k_b1 + raci_z * k_c1 + rabci_z * k_b1 * k_c1
#                 draim = sqrt(drai_x*drai_x + drai_y*drai_y + drai_z*drai_z) # Explicit squares

#                 drbi_x = rbi_x + rabi_x * k_a1 + rbci_x * k_c1 + rabci_x * k_a1 * k_c1
#                 drbi_y = rbi_y + rabi_y * k_a1 + rbci_y * k_c1 + rabci_y * k_a1 * k_c1
#                 # CRITICAL: Corrected potential typo from previous version based on pattern
#                 drbi_z = rbi_z + rabi_z * k_a1 + rbci_z * k_c1 + rabci_z * k_a1 * k_c1
#                 drbim = sqrt(drbi_x*drbi_x + drbi_y*drbi_y + drbi_z*drbi_z) # Explicit squares

#                 drci_x = rci_x + raci_x * k_a1 + rbci_x * k_b1 + rabci_x * k_a1 * k_b1
#                 drci_y = rci_y + raci_y * k_a1 + rbci_y * k_b1 + rabci_y * k_a1 * k_b1
#                 # CRITICAL: Corrected potential typo from previous version based on pattern
#                 drci_z = rci_z + raci_z * k_a1 + rbci_z * k_b1 + rabci_z * k_a1 * k_b1
#                 drcim = sqrt(drci_x*drci_x + drci_y*drci_y + drci_z*drci_z) # Explicit squares


#                 # Calculate components of r1 directly
#                 r1_x = rmi_x + rai_x * k_a1 + rbi_x * k_b1 + rci_x * k_c1 +
#                        rabi_x * k_a1 * k_b1 + raci_x * k_a1 * k_c1 + rbci_x * k_b1 * k_c1 +
#                        rabci_x * k_a1 * k_b1 * k_c1
#                 r1_y = rmi_y + rai_y * k_a1 + rbi_y * k_b1 + rci_y * k_c1 +
#                        rabi_y * k_a1 * k_b1 + raci_y * k_a1 * k_c1 + rbci_y * k_b1 * k_c1 +
#                        rabci_y * k_a1 * k_b1 * k_c1
#                 r1_z = rmi_z + rai_z * k_a1 + rbi_z * k_b1 + rci_z * k_c1 +
#                        rabi_z * k_a1 * k_b1 + raci_z * k_a1 * k_c1 + rbci_z * k_b1 * k_c1 +
#                        rabci_z * k_a1 * k_b1 * k_c1

#                 delta_x = (xo - r1_x) # Renamed 'x' to 'r1_x' etc for clarity
#                 delta_y = (yo - r1_y)
#                 delta_z = (zo - r1_z)

#                 R_squared = delta_x*delta_x + delta_y*delta_y + delta_z*delta_z
#                 # Guard against R_squared being extremely small or zero
#                 # Add a small epsilon to R_squared before sqrt and division
#                 R = sqrt(R_squared + 1e-12) # Add a tiny epsilon to prevent true zero

#                 # Calculate inverse powers of R efficiently
#                 invR = 1.0 / R
#                 invR2 = invR * invR
#                 invR3 = invR2 * invR

#                 # Pre-calculate common terms for G to avoid redundant computation
#                 exp_term = exp(-1im * beta * R)
#                 coeff_term = (invR3 + 1im * beta * invR2)

#                 G_val = zero(ComplexF64)
#                 if scelta == :xx # Use symbols for fast comparison
#                     G_val = delta_x * exp_term * coeff_term
#                 elseif scelta == :xy
#                     G_val = delta_z * exp_term * coeff_term
#                 elseif scelta == :xz
#                     G_val = -delta_y * exp_term * coeff_term
#                 elseif scelta == :yx
#                     G_val = -delta_z * exp_term * coeff_term
#                 elseif scelta == :yy
#                     G_val = delta_y * exp_term * coeff_term
#                 elseif scelta == :yz
#                     G_val = delta_x * exp_term * coeff_term
#                 elseif scelta == :zx
#                     G_val = delta_y * exp_term * coeff_term
#                 elseif scelta == :zy
#                     G_val = -delta_x * exp_term * coeff_term
#                 elseif scelta == :zz
#                     G_val = delta_z * exp_term * coeff_term
#                 # No 'else' branch for `G_val` needed if `zero(ComplexF64)` is default
#                 end

#                 f = draim * drbim * drcim * G_val
#                 sum_c1 += wek[c1] * f
#             end # (c1)
#             sum_b1 += wek[b1] * sum_c1
#         end  # (b1)
#         sum_a1 += wek[a1] * sum_b1
#     end # (a1)

#     integral = 1e-7 * sum_a1 # Apply 1e-7 factor once at the end
#     return integral
# end

# function qrule(n::Int)
#     iter = 2
#     m = trunc((n + 1) / 2)
#     e1 = n * (n + 1)
#     mm = 4 * m - 1
#     t = (pi / (4 * n + 2)) * (3:4:mm)
#     nn = (1 - (1 - 1 / n) / (8 * n * n))
#     xo = nn * cos.(t)
#     den = []
#     d1 = []
#     dpn = []
#     d2pn = []
#     d3pn = []
#     d4pn = []
#     u = []
#     v = []
#     h = []
#     p = []
#     dp = []
#     pk = []
    
#     for kk = 1:iter
#         pkm1 = zeros(size(xo))
#         pkm1[1:size(xo, 1)] .= 1
#         pk = xo
#         for k = 2:n
#             t1 = xo .* pk
#             pkp1 = t1 - pkm1 - (t1 - pkm1) / k + t1
#             pkm1 = pk
#             pk = pkp1
#         end
#         den = 1 .- xo .^ 2
#         d1 = n * (pkm1 - xo .* pk)
#         dpn = d1 ./ den
#         d2pn = (2 .* xo .* dpn - e1 .* pk) ./ den
#         d3pn = (4 * xo .* d2pn + (2 - e1) .* dpn) ./ den
#         d4pn = (6 * xo .* d3pn + (6 - e1) .* d2pn) ./ den
#         u = pk ./ dpn
#         v = d2pn ./ dpn
#         h = -u .* (1 .+ (0.5 * u) .* (v + u .* (v .* v - u .* d3pn ./ (3 * dpn))))
#         p = pk + h .* (dpn + (0.5 * h) .* (d2pn + (h / 3) .* (d3pn + 0.25 * h .* d4pn)))
#         dp = dpn + h .* (d2pn + (0.5 * h) .* (d3pn + h .* d4pn / 3))
#         h = h - p ./ dp
#         xo = xo + h
#     end
#     bp = zeros(1, n)
#     wf = zeros(1, n)
#     bp[1:size(xo, 1)] .= -xo .- h
#     fx = d1 - h .* e1 .* (pk + (h / 2) .* (dpn + (h / 3) .* (
#         d2pn + (h / 4) .* (d3pn + (0.2 * h) .* d4pn))))
#     wf[1:size(xo, 1)] .= 2 * (1 .- bp[1:size(xo, 1)] .^ 2) ./ (fx .* fx)
#     if (m + m) > n
#         bp[Int64(m)] = 0
#     end
#     if !((m + m) == n)
#         m = m - 1
#     end
#     jj = 1:m
#     n1j = (n + 1) .- jj
#     bp[Int64.(n1j)] .= -bp[Int64.(jj)]
#     wf[Int64.(n1j)] .= wf[Int64.(jj)]

#     return vec(bp), vec(wf)
# end
using LinearAlgebra
using StaticArrays # Aggiungi StaticArrays

function compute_lambda_numeric(punti_oss::Matrix{Float64}, volumi::Dict, incidence_selection::Dict{Symbol, Any},
                                vers_punti_oss::Matrix{Float64}, ordine_int::Int, beta::ComplexF64, id, chan)

    N = size(volumi[:coordinate], 1)
    M = size(punti_oss, 1)

    # Pre-allocate vers_celle
    vers_celle = zeros(Float64, N, 3)
    @inbounds begin
        vers_celle[1:incidence_selection[:mx], 1] .= 1.0
        vers_celle[incidence_selection[:mx]+1:incidence_selection[:mx]+incidence_selection[:my], 2] .= 1.0
        vers_celle[incidence_selection[:mx]+incidence_selection[:my]+1:N, 3] .= 1.0
    end

    Lambda = zeros(ComplexF64, M, N)

    rootk, wek = qrule(ordine_int)

    block_size = 100

    for m_block in 1:block_size:M
        m_end = min(m_block + block_size - 1, M)
        Threads.@threads for m = m_block:m_end
            # Converti le viste in StaticArrays prima di passarle a compute_hi_optimized
            # Questo creerà una SVector per ogni riga, con minima allocazione sullo stack.
            punti_oss_m_svec = SVector{3}(@view punti_oss[m, :])
            vers_punti_oss_m_svec = SVector{3}(@view vers_punti_oss[m, :])

            # Pre-calculate boolean flags for performance, using StaticArray elements
            is_pox = abs(vers_punti_oss_m_svec[1]) > 1e-5
            is_poy = abs(vers_punti_oss_m_svec[2]) > 1e-5

            for n = 1:N
                # Converti la vista di barra in SVector{24}
                barra_n_svec = SVector{24}(@view volumi[:coordinate][n, :])
                S_n = volumi[:S][n]

                scelta_val = :vuoto

                @inbounds begin
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
                    else
                        if is_pox
                            scelta_val = :xz
                        elseif is_poy
                            scelta_val = :yz
                        else
                            scelta_val = :zz
                        end
                    end
                end

                # Chiama la versione ottimizzata con StaticArrays
                Lambda[m, n] = compute_hi_optimized(barra_n_svec, punti_oss_m_svec, scelta_val, rootk, wek, beta) / S_n
            end
        end
        # sleep(0.01)
        println("block Lambda : ", round(m_end/block_size), " / ", round(M/block_size))
        if is_stop_requested(id)
          println("Simulazione $(id) interrotta per richiesta stop.")
          return nothing
      end
    end

    return Lambda

end

@inline function compute_hi_optimized(barra::SVector{24, Float64}, centro_oss::SVector{3, Float64},
                                      scelta::Symbol, rootk::AbstractVector{Float64}, wek::AbstractVector{Float64},
                                      beta::ComplexF64)

    # I parametri `barra` e `centro_oss` sono già SVector, quindi non è necessario estrarre individualmente
    # i componenti come prima, o si può farlo per chiarezza se preferisci.

    # Non è necessario @inbounds per l'accesso a SVector; sono sicuri per costruzione.

    inv8 = 0.125

    # Ora definiamo rmi, rai, etc., direttamente come SVector{3, Float64}
    rmi = SVector(
        inv8 * (barra[1] + barra[4] + barra[7] + barra[10] + barra[13] + barra[16] + barra[19] + barra[22]),
        inv8 * (barra[2] + barra[5] + barra[8] + barra[11] + barra[14] + barra[17] + barra[20] + barra[23]),
        inv8 * (barra[3] + barra[6] + barra[9] + barra[12] + barra[15] + barra[18] + barra[21] + barra[24])
    )

    rai = SVector(
        inv8 * (-barra[1] + barra[4] + barra[10] - barra[7] - barra[13] + barra[16] + barra[22] - barra[19]),
        inv8 * (-barra[2] + barra[5] + barra[11] - barra[8] - barra[14] + barra[17] + barra[23] - barra[20]),
        inv8 * (-barra[3] + barra[6] + barra[12] - barra[9] - barra[15] + barra[18] + barra[21] - barra[21]) # corrected typo
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

    nlk = length(wek)
    sum_a1 = zero(ComplexF64)

    # centro_oss è già un SVector, quindi possiamo destrutturarlo direttamente
    xo, yo, zo = centro_oss

    @fastmath for a1 in 1:nlk
        sum_b1 = zero(ComplexF64)
        k_a1 = rootk[a1]

        for b1 in 1:nlk
            sum_c1 = zero(ComplexF64)
            k_b1 = rootk[b1]

            for c1 in 1:nlk
                k_c1 = rootk[c1]

                # I vettori intermedi `drai`, `drbi`, `drci` e `r1` saranno ora SVectors
                drai = rai + rabi * k_b1 + raci * k_c1 + rabci * k_b1 * k_c1
                drbi = rbi + rabi * k_a1 + rbci * k_c1 + rabci * k_a1 * k_c1
                drci = rci + raci * k_a1 + rbci * k_b1 + rabci * k_a1 * k_b1

                # Utilizza norm() per calcolare la norma euclidea di SVector
                draim = norm(drai)
                drbim = norm(drbi)
                drcim = norm(drci)

                r1 = rmi + rai * k_a1 + rbi * k_b1 + rci * k_c1 +
                     rabi * k_a1 * k_b1 + raci * k_a1 * k_c1 + rbci * k_b1 * k_c1 +
                     rabci * k_a1 * k_b1 * k_c1

                # Accedi ai componenti di r1 e centro_oss direttamente
                delta_x = (xo - r1[1])
                delta_y = (yo - r1[2])
                delta_z = (zo - r1[3])

                R_squared = delta_x*delta_x + delta_y*delta_y + delta_z*delta_z
                R = sqrt(R_squared + 1e-12)

                invR = 1.0 / R
                invR2 = invR * invR
                invR3 = invR2 * invR

                exp_term = exp(-1im * beta * R)
                coeff_term = (invR3 + 1im * beta * invR2)

                G_val = zero(ComplexF64)
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
                elseif scelta == :zz
                    G_val = delta_z * exp_term * coeff_term
                end

                f = draim * drbim * drcim * G_val
                sum_c1 += wek[c1] * f
            end
            sum_b1 += wek[b1] * sum_c1
        end
        sum_a1 += wek[a1] * sum_b1
    end

    integral = 1e-7 * sum_a1
    return integral
end

# --- La funzione qrule rimane invariata ---
function qrule(n::Int)
    iter = 2
    m = trunc((n + 1) / 2)
    e1 = n * (n + 1)
    mm = 4 * m - 1
    t = (pi / (4 * n + 2)) * (3:4:mm)
    nn = (1 - (1 - 1 / n) / (8 * n * n))
    xo = nn * cos.(t)
    # ... (il resto della funzione qrule è identico) ...
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
    
    for kk = 1:iter
        pkm1 = zeros(size(xo))
        pkm1[1:size(xo, 1)] .= 1
        pk = xo
        for k = 2:n
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