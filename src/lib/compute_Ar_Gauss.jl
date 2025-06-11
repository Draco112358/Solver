# using LinearAlgebra

# function compute_Ar_Gauss(barre, centriOss, ordine, beta, simulation_id, chan)
# 	numCentri = size(centriOss, 1)
# 	numBarre = size(barre, 1)

# 	ha = zeros(ComplexF64, numBarre, numCentri)
# 	rootkx, wekx = qrule(ordine)
# 	rootky, weky = qrule(ordine)
# 	rootkz, wekz = qrule(ordine)

# 	block_size = 100
# 	for m_block in 1:block_size:numBarre
# 		m_end = min(m_block + block_size - 1, numBarre)
# 		Base.Threads.@threads for cont in m_block:m_end
# 			barra = barre[cont, :]

# 			xb = barra[[1, 4, 7, 10, 13, 16, 19, 22]]
# 			yb = barra[[2, 5, 8, 11, 14, 17, 20, 23]]
# 			zb = barra[[3, 6, 9, 12, 15, 18, 21, 24]]

# 			x_bar = [minimum(xb), maximum(xb)]
# 			y_bar = [minimum(yb), maximum(yb)]
# 			z_bar = [minimum(zb), maximum(zb)]

# 			for cc in 1:numCentri
# 				x_o = centriOss[cc, 1]
# 				y_o = centriOss[cc, 2]
# 				z_o = centriOss[cc, 3]

# 				ha[cont, cc] = compute_ha(x_o, x_bar, y_o, y_bar, z_o, z_bar, ordine, beta, rootkx, wekx, rootky, weky, rootkz, wekz)
# 				if cc % 100 == 0
# 					yield()  # Permette alla task dell'heartbeat di essere schedulata
# 				end
# 			end
# 		end
# 		sleep(0)
# 		println("block Ar : ", round(m_end / block_size), " / ", round(numBarre / block_size))
# 		if is_stop_requested(simulation_id)
# 			println("Simulazione $(simulation_id) interrotta per richiesta stop.")
# 			return nothing # O un altro valore che indica interruzione
# 		end
# 	end


# 	if is_stop_requested(simulation_id)
# 		println("Simulazione $(simulation_id) interrotta per richiesta stop.")
# 		return nothing # O un altro valore che indica interruzione
# 	else
# 		return ha
# 	end
# end

# function compute_ha(xo::Float64, x_vect_bar::Vector{Float64}, yo::Float64, y_vect_bar::Vector{Float64}, zo::Float64, z_vect_bar::Vector{Float64}, ordine::Int, beta::Float64, rootkx::Vector{Float64}, wekx::Vector{Float64}, rootky::Vector{Float64}, weky::Vector{Float64}, rootkz::Vector{Float64}, wekz::Vector{Float64})
# 	x1 = x_vect_bar[1]
# 	x2 = x_vect_bar[2]
# 	y1 = y_vect_bar[1]
# 	y2 = y_vect_bar[2]
# 	z1 = z_vect_bar[1]
# 	z2 = z_vect_bar[2]
# 	barra = [x1 y1 z1; x2 y1 z1; x1 y2 z1; x2 y2 z1; x1 y1 z2; x2 y1 z2; x1 y2 z2; x2 y2 z2]

# 	xi1 = barra[[1, 2, 3, 4], 1]
# 	yi1 = barra[[1, 2, 3, 4], 2]
# 	zi1 = barra[[1, 2, 3, 4], 3]
# 	xi2 = barra[[5, 6, 7, 8], 1]
# 	yi2 = barra[[5, 6, 7, 8], 2]
# 	zi2 = barra[[5, 6, 7, 8], 3]

# 	# vectors pointing to the vertices of the quadrilateral i
# 	ri = zeros(ComplexF64, 8, 3)
# 	ri[1, :] = [xi1[1], yi1[1], zi1[1]]
# 	ri[2, :] = [xi1[2], yi1[2], zi1[2]]
# 	ri[3, :] = [xi1[3], yi1[3], zi1[3]]
# 	ri[4, :] = [xi1[4], yi1[4], zi1[4]]

# 	ri[5, :] = [xi2[1], yi2[1], zi2[1]]
# 	ri[6, :] = [xi2[2], yi2[2], zi2[2]]
# 	ri[7, :] = [xi2[3], yi2[3], zi2[3]]
# 	ri[8, :] = [xi2[4], yi2[4], zi2[4]]

# 	# nuovo approccio
# 	rmi = vec(0.125 * sum(ri, dims = 1))
# 	rai = 0.125 * (-ri[1, :] + ri[2, :] + ri[4, :] - ri[3, :] - ri[5, :] + ri[6, :] + ri[8, :] - ri[7, :])
# 	rbi = 0.125 * (-ri[1, :] - ri[2, :] + ri[4, :] + ri[3, :] - ri[5, :] - ri[6, :] + ri[8, :] + ri[7, :])
# 	rci = 0.125 * (-ri[1, :] - ri[2, :] - ri[4, :] - ri[3, :] + ri[5, :] + ri[6, :] + ri[8, :] + ri[7, :])
# 	rabi = 0.125 * (ri[1, :] - ri[2, :] + ri[4, :] - ri[3, :] + ri[5, :] - ri[6, :] + ri[8, :] - ri[7, :])
# 	rbci = 0.125 * (ri[1, :] + ri[2, :] - ri[4, :] - ri[3, :] - ri[5, :] - ri[6, :] + ri[8, :] + ri[7, :])
# 	raci = 0.125 * (ri[1, :] - ri[2, :] - ri[4, :] + ri[3, :] - ri[5, :] + ri[6, :] + ri[8, :] - ri[7, :])
# 	rabci = 0.125 * (-ri[1, :] + ri[2, :] - ri[4, :] + ri[3, :] + ri[5, :] - ri[6, :] + ri[8, :] - ri[7, :])

# 	nlkx = length(wekx)
# 	nlky = length(weky)
# 	nlkz = length(wekz)

# 	sum_a1 = 0.0 + 0.0im
# 	for a1 in 1:nlkx
# 		sum_b1 = 0.0 + 0.0im
# 		for b1 in 1:nlky
# 			sum_c1 = 0.0 + 0.0im
# 			for c1 in 1:nlkz
# 				drai = rai + rabi * rootky[b1] + raci * rootkz[c1] + rabci * rootky[b1] * rootkz[c1]
# 				drbi = rbi + rabi * rootkx[a1] + rbci * rootkz[c1] + rabci * rootkx[a1] * rootkz[c1]
# 				drci = rci + raci * rootkx[a1] + rbci * rootky[b1] + rabci * rootkx[a1] * rootky[b1]
# 				draim = norm(drai)
# 				drbim = norm(drbi)
# 				drcim = norm(drci)

# 				r1 = rmi + rai * rootkx[a1] + rbi * rootky[b1] + rci * rootkz[c1] + rabi * rootkx[a1] * rootky[b1] + raci * rootkx[a1] * rootkz[c1] + rbci * rootky[b1] * rootkz[c1] +
# 					 rabci * rootkx[a1] * rootky[b1] * rootkz[c1]

# 				x = real(r1[1])
# 				y = real(r1[2])
# 				z = real(r1[3])

# 				delta_x = (xo - x)
# 				delta_y = (yo - y)
# 				delta_z = (zo - z)

# 				R = sqrt(delta_x^2 + delta_y^2 + delta_z^2)

# 				G = (1 / R) * exp(-1im * beta * R)

# 				f = draim * drbim * drcim * G

# 				sum_c1 += wekz[c1] * f
# 			end   # (c1)
# 			sum_b1 += weky[b1] * sum_c1
# 		end  # (b1)
# 		sum_a1 += wekx[a1] * sum_b1
# 	end # (a1)

# 	return sum_a1
# end

# function qrule(n::Int)
# 	iter = 2
# 	m = trunc((n + 1) / 2)
# 	e1 = n * (n + 1)
# 	mm = 4 * m - 1
# 	t = (pi / (4 * n + 2)) * (3:4:mm)
# 	nn = (1 - (1 - 1 / n) / (8 * n * n))
# 	xo = nn * cos.(t)
# 	den = []
# 	d1 = []
# 	dpn = []
# 	d2pn = []
# 	d3pn = []
# 	d4pn = []
# 	u = []
# 	v = []
# 	h = []
# 	p = []
# 	dp = []
# 	pk = []

# 	for kk ∈ 1:iter
# 		pkm1 = zeros(size(xo))
# 		pkm1[1:size(xo, 1)] .= 1
# 		pk = xo
# 		for k ∈ 2:n
# 			t1 = xo .* pk
# 			pkp1 = t1 - pkm1 - (t1 - pkm1) / k + t1
# 			pkm1 = pk
# 			pk = pkp1
# 		end
# 		den = 1 .- xo .^ 2
# 		d1 = n * (pkm1 - xo .* pk)
# 		dpn = d1 ./ den
# 		d2pn = (2 .* xo .* dpn - e1 .* pk) ./ den
# 		d3pn = (4 * xo .* d2pn + (2 - e1) .* dpn) ./ den
# 		d4pn = (6 * xo .* d3pn + (6 - e1) .* d2pn) ./ den
# 		u = pk ./ dpn
# 		v = d2pn ./ dpn
# 		h = -u .* (1 .+ (0.5 * u) .* (v + u .* (v .* v - u .* d3pn ./ (3 * dpn))))
# 		p = pk + h .* (dpn + (0.5 * h) .* (d2pn + (h / 3) .* (d3pn + 0.25 * h .* d4pn)))
# 		dp = dpn + h .* (d2pn + (0.5 * h) .* (d3pn + h .* d4pn / 3))
# 		h = h - p ./ dp
# 		xo = xo + h
# 	end
# 	bp = zeros(1, n)
# 	wf = zeros(1, n)
# 	bp[1:size(xo, 1)] .= -xo .- h
# 	fx = d1 - h .* e1 .* (pk + (h / 2) .* (dpn + (h / 3) .* (
# 		d2pn + (h / 4) .* (d3pn + (0.2 * h) .* d4pn))))
# 	wf[1:size(xo, 1)] .= 2 * (1 .- bp[1:size(xo, 1)] .^ 2) ./ (fx .* fx)
# 	if (m + m) > n
# 		bp[Int64(m)] = 0
# 	end
# 	if !((m + m) == n)
# 		m = m - 1
# 	end
# 	jj = 1:m
# 	n1j = (n + 1) .- jj
# 	bp[Int64.(n1j)] .= -bp[Int64.(jj)]
# 	wf[Int64.(n1j)] .= wf[Int64.(jj)]

# 	return vec(bp), vec(wf)
# end
using LinearAlgebra
using StaticArrays # Necessario per SVector, che useremo per alcune pre-allocazioni
using Base.Threads
using Printf # Per un output formattato

# Assicurati che `qrule` e `compute_ha_optimized` siano definite e disponibili.
# --- ADATTAMENTO NECESSARIO PER compute_ha_optimized ---
# La firma di compute_ha_optimized dovrà cambiare per accettare 6 Float64 per i min/max.
# Esempio:
# function compute_ha_optimized(xo::Float64, x_min_bar::Float64, x_max_bar::Float64,
#                               yo::Float64, y_min_bar::Float64, y_max_bar::Float64,
#                               zo::Float64, z_min_bar::Float64, z_max_bar::Float64,
#                               beta::ComplexF64, rootkx, wekx, rootky, weky, rootkz, wekz)
#     # ... il corpo della funzione rimane simile, usando direttamente x_min_bar, x_max_bar, ecc.
# end
# --------------------------------------------------------

# function compute_Ar_Gauss(barre::Union{Transpose{Float64, Matrix{Float64}}, Transpose{Real, Matrix{Real}}}, # Accetta Matrix{Float64} direttamente
#                                     centriOss::Matrix{Float64}, ordine::Int, beta::ComplexF64,
#                                     simulation_id=nothing, chan=nothing)

#     # Converti barre a Matrix{Float64} se è Transpose.
#     # Questo assicura che size(barre, 1) sia il numero di barre
#     # e che l'accesso alle righe sia efficiente.
#     barre = barre isa Transpose ? collect(barre) : barre # collect() crea una copia, ma evita problemi di Transpose
#     # Oppure, se vuoi evitare la copia, potresti fare `size(barre, 2)` se è Transpose e `barre[1,cont]` ecc.
#     # ma `collect` è spesso più semplice per la compatibilità con il resto del codice.

#     numCentri = size(centriOss, 1)
#     numBarre = size(barre, 1) # Usa barre per il numero di barre

#     # --- Pre-calcoli globali (una volta per la funzione) ---

#     # 1. qrule values (già fuori, ottimo)
#     rootkx, wekx = qrule(ordine)
#     rootky, weky = qrule(ordine)
#     rootkz, wekz = qrule(ordine)

#     # 2. Pre-calcola i dati min/max per tutte le barre
#     # Questi array conterranno i singoli valori min/max per ogni barra.
#     x_min_bar_arr = Vector{Float64}(undef, numBarre)
#     x_max_bar_arr = Vector{Float64}(undef, numBarre)
#     y_min_bar_arr = Vector{Float64}(undef, numBarre)
#     y_max_bar_arr = Vector{Float64}(undef, numBarre)
#     z_min_bar_arr = Vector{Float64}(undef, numBarre)
#     z_max_bar_arr = Vector{Float64}(undef, numBarre)

#     @inbounds Threads.@threads for cont = 1:numBarre
#         # Pre-estrai la riga come SVector per calcoli efficienti su piccola scala
#         barra_svec_24 = SVector{24, Float64}(@view barre[cont, :]) # Assuming 24 elements per bar

#         # Estrai tutti i valori x, y, z e trova min/max
#         # Vettori temporanei per min/max (o usare iteratori direttamente)
#         x_coords_all_svec = SVector{8, Float64}(
#             barra_svec_24[1], barra_svec_24[4], barra_svec_24[7], barra_svec_24[10],
#             barra_svec_24[13], barra_svec_24[16], barra_svec_24[19], barra_svec_24[22]
#         )
#         y_coords_all_svec = SVector{8, Float64}(
#             barra_svec_24[2], barra_svec_24[5], barra_svec_24[8], barra_svec_24[11],
#             barra_svec_24[14], barra_svec_24[17], barra_svec_24[20], barra_svec_24[23]
#         )
#         z_coords_all_svec = SVector{8, Float64}(
#             barra_svec_24[3], barra_svec_24[6], barra_svec_24[9], barra_svec_24[12],
#             barra_svec_24[15], barra_svec_24[18], barra_svec_24[21], barra_svec_24[24]
#         )

#         x_min_bar_arr[cont] = minimum(x_coords_all_svec)
#         x_max_bar_arr[cont] = maximum(x_coords_all_svec)
#         y_min_bar_arr[cont] = minimum(y_coords_all_svec)
#         y_max_bar_arr[cont] = maximum(y_coords_all_svec)
#         z_min_bar_arr[cont] = minimum(z_coords_all_svec)
#         z_max_bar_arr[cont] = maximum(z_coords_all_svec)
#     end

#     # 3. Pre-calcola le coordinate x_o, y_o, z_o per tutti i centri di osservazione
#     xo_arr = Vector{Float64}(undef, numCentri)
#     yo_arr = Vector{Float64}(undef, numCentri)
#     zo_arr = Vector{Float64}(undef, numCentri)

#     @inbounds Threads.@threads for cc = 1:numCentri
#         xo_arr[cc] = centriOss[cc, 1]
#         yo_arr[cc] = centriOss[cc, 2]
#         zo_arr[cc] = centriOss[cc, 3]
#     end

#     # Pre-allocazione della matrice finale `ha`
#     ha = zeros(ComplexF64, numBarre, numCentri)

#     # --- Loop di elaborazione parallelizzata ---
#     block_size = 20 # Adjust this based on your hardware and problem size

#     @inbounds for m_block_start in 1:block_size:numBarre
#         m_end = min(m_block_start + block_size - 1, numBarre)
#         # Use Threads.@threads for parallelization
#         Threads.@threads for cont in m_block_start:m_end
#             # Accedi ai valori min/max pre-calcolati per la barra corrente
#             current_x_min_bar = x_min_bar_arr[cont]
#             current_x_max_bar = x_max_bar_arr[cont]
#             current_y_min_bar = y_min_bar_arr[cont]
#             current_y_max_bar = y_max_bar_arr[cont]
#             current_z_min_bar = z_min_bar_arr[cont]
#             current_z_max_bar = z_max_bar_arr[cont]

#             for cc in 1:numCentri
#                 # Accedi alle coordinate pre-calcolate del centro di osservazione
#                 x_o = xo_arr[cc]
#                 y_o = yo_arr[cc]
#                 z_o = zo_arr[cc]

#                 ha[cont, cc] = compute_ha_optimized(
#                     x_o, current_x_min_bar, current_x_max_bar,
#                     y_o, current_y_min_bar, current_y_max_bar,
#                     z_o, current_z_min_bar, current_z_max_bar,
#                     beta, rootkx, wekx, rootky, weky, rootkz, wekz
#                 )
#             end
#         end
#         sleep(0) # Keep for progress reporting and yielding to other Julia tasks
#         @printf "Block Ar : %.0f / %.0f\n" round(m_end / block_size) round(numBarre / block_size)
#         if is_stop_requested(simulation_id)
#             println("Simulazione $(simulation_id) interrotta per richiesta stop.")
#             return nothing
#         end
#     end

#     if is_stop_requested(simulation_id)
#         println("Simulazione $(simulation_id) interrotta per richiesta stop.")
#         return nothing
#     else
#         return ha
#     end
# end

# # --- Optimized `compute_ha` ---
# # Make it `@inline` to encourage Julia to insert its code directly at the call site,
# # potentially eliminating function call overhead.
# @inline function compute_ha_optimized(
#     xo::Float64, x1::Float64, x2::Float64,
#     yo::Float64, y1::Float64, y2::Float64,
#     zo::Float64, z1::Float64, z2::Float64,
#     beta::ComplexF64, # Beta must be ComplexF64 for `exp(-1im * beta * R)`
#     rootkx::AbstractVector{Float64}, wekx::AbstractVector{Float64},
#     rootky::AbstractVector{Float64}, weky::AbstractVector{Float64},
#     rootkz::AbstractVector{Float64}, wekz::AbstractVector{Float64}
# )
#     # No need to recalculate x1, x2, etc. from vectors; they are passed directly.
#     # No need to create `barra` matrix or its sub-arrays (xi1, yi1, etc.) here.
#     # The 8 vectors ri[i, :] are just fixed vertices of a rectangular prism.
#     # We can define them based on x1, x2, y1, y2, z1, z2.

#     # These are fixed for a rectangular prism defined by min/max coordinates
#     # and can be hardcoded or calculated once per call to compute_ha.
#     # Using StaticArrays for these small vectors is highly recommended if `StaticArrays.jl` is added.
#     # Otherwise, regular Vectors are fine, but be mindful of allocations if they change per call.
#     # Since these are constant patterns based on x1,x2,y1,y2,z1,z2, they don't allocate in the inner loops.

#     # Pre-calculate the 'r' vectors (rmi, rai, etc.) using the x1,x2,y1,y2,z1,z2 directly.
#     # This avoids creating `barra` and `ri` arrays, which cause significant allocations.
#     # Each of these `r` vectors is a 3-element vector.
#     # For maximum performance, if `StaticArrays.jl` is in use:
#     # using StaticArrays
#     # rmi = SVector{3,ComplexF64}(...)
#     # ... and so on.
#     # If not using StaticArrays, normal Vectors are fine, but note the initial allocation cost.
#     # For clarity, I'll write them out as direct calculations.

#     # Original 'ri' logic based on x1,y1,z1 and x2,y2,z2.
#     # Assuming the original barra definition for ri:
#     # ri[1, :] = [x1, y1, z1]
#     # ri[2, :] = [x2, y1, z1]
#     # ri[3, :] = [x1, y2, z1]
#     # ri[4, :] = [x2, y2, z1]
#     # ri[5, :] = [x1, y1, z2]
#     # ri[6, :] = [x2, y1, z2]
#     # ri[7, :] = [x1, y2, z2]
#     # ri[8, :] = [x2, y2, z2]

#     # Calculate rmi, rai, etc. directly without intermediate `ri` array.
#     # Each component (x,y,z) is calculated separately for performance.
#     # The `0.125` can be pre-multiplied for efficiency.
#     # Also, ensuring results are ComplexF64 from the start for type stability.

#     # rmi (mean point)
#     # Equivalent to 0.125 * sum(ri, dims=1)
#     rmi_x = 0.125 * (x1 + x2 + x1 + x2 + x1 + x2 + x1 + x2)
#     rmi_y = 0.125 * (y1 + y1 + y2 + y2 + y1 + y1 + y2 + y2)
#     rmi_z = 0.125 * (z1 + z1 + z1 + z1 + z2 + z2 + z2 + z2)

#     # rai (partial derivatives)
#     # 0.125 * (-ri[1, :] + ri[2, :] + ri[4, :] - ri[3, :] - ri[5, :] + ri[6, :] + ri[8, :] - ri[7, :])
#     rai_x = 0.125 * (-x1 + x2 + x2 - x1 - x1 + x2 + x2 - x1)
#     rai_y = 0.125 * (-y1 + y1 + y2 - y2 - y1 + y1 + y2 - y2)
#     rai_z = 0.125 * (-z1 + z1 + z1 - z1 - z2 + z2 + z2 - z2)

#     # rbi
#     # 0.125 * (-ri[1, :] - ri[2, :] + ri[4, :] + ri[3, :] - ri[5, :] - ri[6, :] + ri[8, :] + ri[7, :])
#     rbi_x = 0.125 * (-x1 - x2 + x2 + x1 - x1 - x2 + x2 + x1)
#     rbi_y = 0.125 * (-y1 - y1 + y2 + y2 - y1 - y1 + y2 + y2)
#     rbi_z = 0.125 * (-z1 - z1 + z1 + z1 - z2 - z2 + z2 + z2)

#     # rci
#     # 0.125 * (-ri[1, :] - ri[2, :] - ri[4, :] - ri[3, :] + ri[5, :] + ri[6, :] + ri[8, :] + ri[7, :])
#     rci_x = 0.125 * (-x1 - x2 - x2 - x1 + x1 + x2 + x2 + x1)
#     rci_y = 0.125 * (-y1 - y1 - y2 - y2 + y1 + y1 + y2 + y2)
#     rci_z = 0.125 * (-z1 - z1 - z1 - z1 + z2 + z2 + z2 + z2)

#     # rabi
#     # 0.125 * (ri[1, :] - ri[2, :] + ri[4, :] - ri[3, :] + ri[5, :] - ri[6, :] + ri[8, :] - ri[7, :])
#     rabi_x = 0.125 * (x1 - x2 + x2 - x1 + x1 - x2 + x2 - x1)
#     rabi_y = 0.125 * (y1 - y1 + y2 - y2 + y1 - y1 + y2 - y2)
#     rabi_z = 0.125 * (z1 - z1 + z1 - z1 + z2 - z2 + z2 - z2)

#     # rbci
#     # 0.125 * (ri[1, :] + ri[2, :] - ri[4, :] - ri[3, :] - ri[5, :] - ri[6, :] + ri[8, :] + ri[7, :])
#     rbci_x = 0.125 * (x1 + x2 - x2 - x1 - x1 - x2 + x2 + x1)
#     rbci_y = 0.125 * (y1 + y1 - y2 - y2 - y1 - y1 + y2 + y2)
#     rbci_z = 0.125 * (z1 + z1 - z1 - z1 - z2 - z2 + z2 + z2)

#     # raci
#     # 0.125 * (ri[1, :] - ri[2, :] - ri[4, :] + ri[3, :] - ri[5, :] + ri[6, :] + ri[8, :] - ri[7, :])
#     raci_x = 0.125 * (x1 - x2 - x2 + x1 - x1 + x2 + x2 - x1)
#     raci_y = 0.125 * (y1 - y1 - y2 + y2 - y1 + y1 + y2 - y2)
#     raci_z = 0.125 * (z1 - z1 - z1 + z1 - z2 + z2 + z2 - z2)

#     # rabci
#     # 0.125 * (-ri[1, :] + ri[2, :] - ri[4, :] + ri[3, :] + ri[5, :] - ri[6, :] + ri[8, :] - ri[7, :])
#     rabci_x = 0.125 * (-x1 + x2 - x2 + x1 + x1 - x2 + x2 - x1)
#     rabci_y = 0.125 * (-y1 + y1 - y2 + y2 + y1 - y1 + y2 - y2)
#     rabci_z = 0.125 * (-z1 + z1 - z1 + z1 + z2 - z2 + z2 - z2)

#     nlkx = length(wekx)
#     nlky = length(weky)
#     nlkz = length(wekz)

#     sum_a1 = zero(ComplexF64) # Initialize with ComplexF64 zero for type stability

#     # Use @fastmath and @inbounds for performance in the nested loops
#     @fastmath @inbounds for a1 in 1:nlkx
#         sum_b1 = zero(ComplexF64)
#         kx_a1 = rootkx[a1] # Store for reuse

#         for b1 in 1:nlky
#             sum_c1 = zero(ComplexF64)
#             ky_b1 = rootky[b1] # Store for reuse

#             for c1 in 1:nlkz
#                 kz_c1 = rootkz[c1] # Store for reuse

#                 # Calculate components of drai, drbi, drci directly
#                 drai_x = rai_x + rabi_x * ky_b1 + raci_x * kz_c1 + rabci_x * ky_b1 * kz_c1
#                 drai_y = rai_y + rabi_y * ky_b1 + raci_y * kz_c1 + rabci_y * ky_b1 * kz_c1
#                 drai_z = rai_z + rabi_z * ky_b1 + raci_z * kz_c1 + rabci_z * ky_b1 * kz_c1
#                 draim = sqrt(drai_x^2 + drai_y^2 + drai_z^2)

#                 drbi_x = rbi_x + rabi_x * kx_a1 + rbci_x * kz_c1 + rabci_x * kx_a1 * kz_c1
#                 drbi_y = rbi_y + rabi_y * kx_a1 + rbci_y * kz_c1 + rabci_y * kx_a1 * kz_c1
#                 drbi_z = rbi_z + raci_z * kx_a1 + rbci_z * ky_b1 + rabci_z * kx_a1 * ky_b1 # Check original formula for drbi_z
#                 drbim = sqrt(drbi_x^2 + drbi_y^2 + drbi_z^2)

#                 drci_x = rci_x + raci_x * kx_a1 + rbci_x * ky_b1 + rabci_x * kx_a1 * ky_b1
#                 drci_y = rci_y + raci_y * kx_a1 + rbci_y * ky_b1 + rabci_y * kx_a1 * ky_b1
#                 drci_z = rci_z + raci_z * kx_a1 + rbci_z * ky_b1 + rabci_z * kx_a1 * ky_b1 # Check original formula for drci_z
#                 drcim = sqrt(drci_x^2 + drci_y^2 + drci_z^2)

#                 # Calculate components of r1 directly
#                 r1_x = rmi_x + rai_x * kx_a1 + rbi_x * ky_b1 + rci_x * kz_c1 +
#                        rabi_x * kx_a1 * ky_b1 + raci_x * kx_a1 * kz_c1 + rbci_x * ky_b1 * kz_c1 +
#                        rabci_x * kx_a1 * ky_b1 * kz_c1
#                 r1_y = rmi_y + rai_y * kx_a1 + rbi_y * ky_b1 + rci_y * kz_c1 +
#                        rabi_y * kx_a1 * ky_b1 + raci_y * kx_a1 * kz_c1 + rbci_y * ky_b1 * kz_c1 +
#                        rabci_y * kx_a1 * ky_b1 * kz_c1
#                 r1_z = rmi_z + rai_z * kx_a1 + rbi_z * ky_b1 + rci_z * kz_c1 +
#                        rabi_z * kx_a1 * ky_b1 + raci_z * kx_a1 * kz_c1 + rbci_z * ky_b1 * kz_c1 +
#                        rabci_z * kx_a1 * ky_b1 * kz_c1

#                 # Original code used `real(r1[1])` - this is unnecessary if r1_x, r1_y, r1_z are Float64.
#                 # All components are real (x,y,z coordinates).
#                 x = r1_x
#                 y = r1_y
#                 z = r1_z

#                 delta_x = (xo - x)
#                 delta_y = (yo - y)
#                 delta_z = (zo - z)

#                 R = sqrt(delta_x^2 + delta_y^2 + delta_z^2)

#                 # Handle potential R=0: add a small epsilon to R if it could be zero.
#                 # If R=0, (1/R) is Inf, exp(-1im*beta*0) = 1, so Inf*1.
#                 # A robust solution needs to handle singularities properly.
#                 # For now, assuming R > 0.
#                 invR = 1.0 / R
#                 G = invR * exp(-1im * beta * R) # `beta` must be ComplexF64 for this to work correctly

#                 # `f` should be ComplexF64 since `G` is ComplexF64.
#                 f = draim * drbim * drcim * G

#                 sum_c1 += wekz[c1] * f
#             end
#             sum_b1 += weky[b1] * sum_c1
#         end
#         sum_a1 += wekx[a1] * sum_b1
#     end

#     return sum_a1
# end

using LinearAlgebra
using StaticArrays # Essenziale per prestazioni elevate con array di piccole dimensioni
using Base.Threads
using Printf

# --- ADATTAMENTO NECESSARIO PER compute_ha_optimized ---
# La firma di compute_ha_optimized deve accettare i singoli valori min/max e
# i componenti dei vettori rmi, rai, rbi, rci, rabi, rbci, raci, rabci.
# NON deve semplificare questi termini a zero internamente.

@inline function compute_ha_optimized_corrected(
    xo::Float64, yo::Float64, zo::Float64,
    rmi_x::Float64, rmi_y::Float64, rmi_z::Float64,
    rai_x::Float64, rai_y::Float64, rai_z::Float64,
    rbi_x::Float64, rbi_y::Float64, rbi_z::Float64,
    rci_x::Float64, rci_y::Float64, rci_z::Float64,
    rabi_x::Float64, rabi_y::Float64, rabi_z::Float64,
    rbci_x::Float64, rbci_y::Float64, rbci_z::Float64,
    raci_x::Float64, raci_y::Float64, raci_z::Float64,
    rabci_x::Float64, rabci_y::Float64, rabci_z::Float64,
    beta::ComplexF64,
    rootkx::Vector{Float64}, wekx::Vector{Float64},
    rootky::Vector{Float64}, weky::Vector{Float64},
    rootkz::Vector{Float64}, wekz::Vector{Float64}
)
    nlkx = length(wekx)
    nlky = length(weky)
    nlkz = length(wekz)

    sum_a1 = zero(ComplexF64)

    @fastmath @inbounds for a1 in 1:nlkx
        sum_b1 = zero(ComplexF64)
        kx_a1 = rootkx[a1]

        for b1 in 1:nlky
            sum_c1 = zero(ComplexF64)
            ky_b1 = rootky[b1]

            for c1 in 1:nlkz
                kz_c1 = rootkz[c1]

                # Calculate components of drai, drbi, drci
                drai_x = rai_x + rabi_x * ky_b1 + raci_x * kz_c1 + rabci_x * ky_b1 * kz_c1
                drai_y = rai_y + rabi_y * ky_b1 + raci_y * kz_c1 + rabci_y * ky_b1 * kz_c1
                drai_z = rai_z + rabi_z * ky_b1 + raci_z * kz_c1 + rabci_z * ky_b1 * kz_c1
                draim = sqrt(drai_x^2 + drai_y^2 + drai_z^2)

                drbi_x = rbi_x + rabi_x * kx_a1 + rbci_x * kz_c1 + rabci_x * kx_a1 * kz_c1
                drbi_y = rbi_y + rabi_y * kx_a1 + rbci_y * kz_c1 + rabci_y * kx_a1 * kz_c1
                drbi_z = rbi_z + raci_z * kx_a1 + rbci_z * ky_b1 + rabci_z * kx_a1 * ky_b1
                drbim = sqrt(drbi_x^2 + drbi_y^2 + drbi_z^2)

                drci_x = rci_x + raci_x * kx_a1 + rbci_x * ky_b1 + rabci_x * kx_a1 * ky_b1
                drci_y = rci_y + raci_y * kx_a1 + rbci_y * ky_b1 + rabci_y * kx_a1 * ky_b1
                drci_z = rci_z + raci_z * kx_a1 + rbci_z * ky_b1 + rabci_z * kx_a1 * ky_b1
                drcim = sqrt(drci_x^2 + drci_y^2 + drci_z^2)

                # Calculate components of r1
                r1_x = rmi_x + rai_x * kx_a1 + rbi_x * ky_b1 + rci_x * kz_c1 +
                       rabi_x * kx_a1 * ky_b1 + raci_x * kx_a1 * kz_c1 + rbci_x * ky_b1 * kz_c1 +
                       rabci_x * kx_a1 * ky_b1 * kz_c1
                r1_y = rmi_y + rai_y * kx_a1 + rbi_y * ky_b1 + rci_y * kz_c1 +
                       rabi_y * kx_a1 * ky_b1 + raci_y * kx_a1 * kz_c1 + rbci_y * ky_b1 * kz_c1 +
                       rabci_y * kx_a1 * ky_b1 * kz_c1
                r1_z = rmi_z + rai_z * kx_a1 + rbi_z * ky_b1 + rci_z * kz_c1 +
                       rabi_z * kx_a1 * ky_b1 + raci_z * kx_a1 * kz_c1 + rbci_z * ky_b1 * kz_c1 +
                       rabci_z * kx_a1 * ky_b1 * kz_c1

                delta_x = xo - r1_x
                delta_y = yo - r1_y
                delta_z = zo - r1_z

                R = sqrt(delta_x^2 + delta_y^2 + delta_z^2)

                invR = 1.0 / R
                G = invR * exp(-1im * beta * R)

                f = draim * drbim * drcim * G
                sum_c1 += wekz[c1] * f
            end
            sum_b1 += weky[b1] * sum_c1
        end
        sum_a1 += wekx[a1] * sum_b1
    end

    return sum_a1
end


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

    # 2. Pre-calcola le coordinate x_o, y_o, z_o per tutti i centri di osservazione
    xo_arr = Vector{Float64}(undef, numCentri)
    yo_arr = Vector{Float64}(undef, numCentri)
    zo_arr = Vector{Float64}(undef, numCentri)
    @inbounds Threads.@threads for cc = 1:numCentri
        xo_arr[cc] = centriOss[cc, 1]
        yo_arr[cc] = centriOss[cc, 2]
        zo_arr[cc] = centriOss[cc, 3]
    end


    # Questi array conterranno i 3 componenti (x,y,z) di ogni 'r' vettore per ciascuna barra
    # li memorizziamo come tuple di 3 Float64 per efficienza, o SVector{3,Float64}
    # Useremo array di SVector{3,Float64} per massimizzare la performance.
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
    block_size = 200

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
                # Accedi alle coordinate pre-calcolate del centro di osservazione
                x_o = xo_arr[cc]
                y_o = yo_arr[cc]
                z_o = zo_arr[cc]

                ha[cont, cc] = compute_ha_optimized_corrected(
                    x_o, y_o, z_o,
                    current_rmi[1], current_rmi[2], current_rmi[3], # Pass as individual components
                    current_rai[1], current_rai[2], current_rai[3],
                    current_rbi[1], current_rbi[2], current_rbi[3],
                    current_rci[1], current_rci[2], current_rci[3],
                    current_rabi[1], current_rabi[2], current_rabi[3],
                    current_rbci[1], current_rbci[2], current_rbci[3],
                    current_raci[1], current_raci[2], current_raci[3],
                    current_rabci[1], current_rabci[2], current_rabci[3],
                    beta, rootkx, wekx, rootky, weky, rootkz, wekz
                )
            end
        end
        sleep(0)
        @printf "Block Ar : %.0f / %.0f\n" round(m_end / block_size) round(numBarre / block_size)
        if is_stop_requested(simulation_id)
            println("Simulazione $(simulation_id) interrotta per richiesta stop.")
            return nothing
        end
    end

    if is_stop_requested(simulation_id)
        println("Simulazione $(simulation_id) interrotta per richiesta stop.")
        return nothing
    else
        return ha
    end
end

# --- qrule function (remains unchanged from previous optimization) ---
# function qrule(n::Int)
#     iter = 2
#     m = trunc((n + 1) / 2)
#     e1 = n * (n + 1)
#     mm = 4 * m - 1
#     t = (pi / (4 * n + 2)) * (3:4:mm)
#     nn = (1 - (1 - 1 / n) / (8 * n * n))
#     xo = nn * cos.(t)

#     for kk ∈ 1:iter
#         pkm1 = ones(size(xo))
#         pk = copy(xo)

#         @inbounds for k ∈ 2:n
#             pkp1 = @. (2k - 1) / k * xo * pk - (k - 1) / k * pkm1
#             pkm1 = pk
#             pk = pkp1
#         end

#         den = 1 .- xo .^ 2
#         d1 = n .* (pkm1 .- xo .* pk)

#         dpn = d1 ./ den
#         d2pn = (2 .* xo .* dpn - e1 .* pk) ./ den
#         d3pn = (4 * xo .* d2pn + (2 - e1) .* dpn) ./ den
#         d4pn = (6 * xo .* d3pn + (6 - e1) .* d2pn) ./ den
#         u = pk ./ dpn
#         v = d2pn ./ dpn
#         h = -u .* (1 .+ (0.5 * u) .* (v .+ u .* (v .* v .- u .* d3pn ./ (3 * dpn))))
#         p = pk .+ h .* (dpn .+ (0.5 * h) .* (d2pn .+ (h ./ 3) .* (d3pn .+ 0.25 .* h .* d4pn)))
#         dp = dpn .+ h .* (d2pn .+ (0.5 .* h) .* (d3pn .+ h .* d4pn ./ 3))
#         h = h .- p ./ dp
#         xo = xo .+ h
#     end

#     bp = Vector{Float64}(undef, n)
#     wf = Vector{Float64}(undef, n)

#     bp .= -xo .- h
#     fx = d1 .- h .* e1 .* (pk .+ (h ./ 2) .* (dpn .+ (h ./ 3) .* (
#         d2pn .+ (h ./ 4) .* (d3pn .+ (0.2 .* h) .* d4pn))))
#     wf .= 2 .* (1 .- bp .^ 2) ./ (fx .* fx)

#     if (m + m) > n
#         bp[m] = 0.0
#     end
#     if !((m + m) == n)
#         m = m - 1
#     end
#     jj = 1:m
#     n1j = (n + 1) .- jj
#     bp[n1j] .= -bp[jj]
#     wf[n1j] .= wf[jj]

#     return bp, wf
# end