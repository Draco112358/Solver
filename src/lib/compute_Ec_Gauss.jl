using LinearAlgebra
using StaticArrays
using Base.Threads
using Printf

function compute_Ec_Gauss(barre::Matrix{Float64}, normale::Vector{Vector{Float64}}, centriOss::Matrix{Float64}, ordine::Int64, beta::ComplexF64, simulation_id::String, chan=nothing)
	num_barre = size(barre, 1)
	num_centri_oss = size(centriOss, 1)
	hc = zeros(ComplexF64, num_barre, 3, num_centri_oss)

	for i in eachindex(normale)
		normale[i] = convert(Vector{Float64}, normale[i])
	end

	block_size = 100
	normale = hcat(normale...)
	rx, wx = qrule(ordine)
	ry, wy = qrule(ordine)
	perm1 = [3, 2, 1, 6, 5, 4, 9, 8, 7, 12, 11, 10]
	perm1_2 = [2, 3, 1, 5, 6, 4, 8, 9, 7, 11, 12, 10]
	perm2 = [1, 3, 2, 4, 6, 5, 7, 9, 8, 10, 12, 11]
	perm2_2 = [3, 1, 2, 6, 4, 5, 9, 7, 8, 12, 10, 11]
	perm3 = [2, 1, 3, 5, 4, 6, 8, 7, 9, 11, 10, 12]
	# 4. Pre-calcola le versioni SVectore di `centriOss` per ogni permutazione
    # Questo è un array di array di SVectore, uno per ogni riga di centriOss, e per ogni permutazione
    centriOss_precomputed = Dict{Symbol, Vector{SVector{3, Float64}}}()
    centriOss_precomputed[:orig] = [SVector{3}(@view centriOss[cc, :]) for cc in 1:num_centri_oss]
    centriOss_precomputed[:perm1] = [SVector{3}(centriOss[cc, 3], centriOss[cc, 2], centriOss[cc, 1]) for cc in 1:num_centri_oss]
    centriOss_precomputed[:perm1_2] = [SVector{3}(centriOss[cc, 2], centriOss[cc, 3], centriOss[cc, 1]) for cc in 1:num_centri_oss]
    centriOss_precomputed[:perm2] = [SVector{3}(centriOss[cc, 1], centriOss[cc, 3], centriOss[cc, 2]) for cc in 1:num_centri_oss]
    centriOss_precomputed[:perm2_2] = [SVector{3}(centriOss[cc, 3], centriOss[cc, 1], centriOss[cc, 2]) for cc in 1:num_centri_oss]
    centriOss_precomputed[:perm3] = [SVector{3}(centriOss[cc, 2], centriOss[cc, 1], centriOss[cc, 3]) for cc in 1:num_centri_oss]

    # 5. Pre-calcola le versioni SVectore di `barre` per ogni permutazione
    # Questo è un array di array di SVectore, uno per ogni riga di barre, e per ogni permutazione
    # NOTA: Assumiamo che le `barre` in `compute_hcz_xy` e `compute_hcx_xy` siano di dimensione 12.
    # Se `barre` ha più colonne (es. 24 come in compute_lambda_numeric), dovrai estrarre solo le prime 12 o quelle rilevanti.
    barre_precomputed = Dict{Symbol, Vector{SVector{12, Float64}}}()
    barre_precomputed[:orig] = [SVector{12}(@view barre[bb, 1:12]) for bb in 1:num_barre] # Assuming only first 12 cols are used by the perms
    barre_precomputed[:perm1] = [SVector{12}(barre[bb, perm1]) for bb in 1:num_barre]
    barre_precomputed[:perm1_2] = [SVector{12}(barre[bb, perm1_2]) for bb in 1:num_barre]
    barre_precomputed[:perm2] = [SVector{12}(barre[bb, perm2]) for bb in 1:num_barre]
    barre_precomputed[:perm2_2] = [SVector{12}(barre[bb, perm2_2]) for bb in 1:num_barre]
    barre_precomputed[:perm3] = [SVector{12}(barre[bb, perm3]) for bb in 1:num_barre]

    # Pre-allocazione della matrice finale `hc`
    hc = zeros(ComplexF64, num_barre, 3, num_centri_oss)
	@inbounds for m_block_start in 1:block_size:num_barre
        m_end = min(m_block_start + block_size - 1, num_barre)

        # Loop interno sulle `barre` all'interno del blocco assegnato al thread
        Base.Threads.@threads for cont in m_block_start:m_end
            # Accediamo ai dati pre-calcolati per `normale`
            norm_x = normale[cont, 1] # normale[1,cont] in 3xN, x-component
            norm_y = normale[cont, 2] # y-component

            if abs(norm_x) > 1e-10
                for cc in 1:num_centri_oss
                    hc[cont, 1, cc] = compute_hcz_xy(barre_precomputed[:perm1][cont], centriOss_precomputed[:perm1][cc], beta, rx, wx, ry, wy)
                    hc[cont, 2, cc] = compute_hcx_xy(barre_precomputed[:perm1_2][cont], centriOss_precomputed[:perm1_2][cc], beta, rx, wx, ry, wy)
                    hc[cont, 3, cc] = compute_hcx_xy(barre_precomputed[:perm1][cont], centriOss_precomputed[:perm1][cc], beta, rx, wx, ry, wy)
                end
            elseif abs(norm_y) > 1e-10
                for cc in 1:num_centri_oss
                    hc[cont, 1, cc] = compute_hcx_xy(barre_precomputed[:perm2][cont], centriOss_precomputed[:perm2][cc], beta, rx, wx, ry, wy)
                    hc[cont, 2, cc] = compute_hcz_xy(barre_precomputed[:perm2][cont], centriOss_precomputed[:perm2][cc], beta, rx, wx, ry, wy)
                    hc[cont, 3, cc] = compute_hcx_xy(barre_precomputed[:perm2_2][cont], centriOss_precomputed[:perm2_2][cc], beta, rx, wx, ry, wy)
                end
            else # Presumiamo che sia la componente Z ad essere non-nulla
                for cc in 1:num_centri_oss
                    hc[cont, 1, cc] = compute_hcx_xy(barre_precomputed[:orig][cont], centriOss_precomputed[:orig][cc], beta, rx, wx, ry, wy)
                    hc[cont, 2, cc] = compute_hcx_xy(barre_precomputed[:perm3][cont], centriOss_precomputed[:perm3][cc], beta, rx, wx, ry, wy)
                    hc[cont, 3, cc] = compute_hcz_xy(barre_precomputed[:orig][cont], centriOss_precomputed[:orig][cc], beta, rx, wx, ry, wy)
                end
            end
        end
        # Output di progresso (potrebbe essere necessario sincronizzare l'output per i thread)
        @printf "Block Ec: %.0f / %.0f\n" round(m_end / block_size) round(num_barre / block_size)
        #Check per la richiesta di stop (se is_stop_requested è thread-safe)
        if is_stop_requested(simulation_id)
            println("Simulazione $(simulation_id) interrotta per richiesta stop.")
            return nothing
        end
    end

	return hc
end

function compute_hcz_xy(barra::SVector{12, Float64}, centriOss::SVector{3, Float64}, beta::ComplexF64, rx::Vector{Float64}, wx::Vector{Float64}, ry::Vector{Float64}, wy::Vector{Float64})

    # Non serve più `numCentri`, dato che il risultato finale è un singolo scalare.
    # numCentri = size(centriOss, 1)

    # Pre-calcolo delle coordinate min/max e degli h
    # (Queste righe non allocano e sono efficienti con StaticArrays)
    x1 = min(barra[1], barra[4], barra[7], barra[10])
    x2 = max(barra[1], barra[4], barra[7], barra[10])
    y1 = min(barra[2], barra[5], barra[8], barra[11])
    y2 = max(barra[2], barra[5], barra[8], barra[11])
    h1 = (x2 - x1) / 2
    h2 = (x2 + x1) / 2
    h3 = (y2 - y1) / 2
    h4 = (y2 + y1) / 2

    nx = length(wx)
    my = length(wy)

    # Inizializza `res2_1` come SCALARE ComplexF64, non un array.
    res2_1_scalar = zero(ComplexF64)

    # Estrai i componenti di centriOss una volta
    x = centriOss[1]
    y = centriOss[2]
    z = centriOss[3]
    zp = barra[3]

    @inbounds for a1 in 1:nx
        # Inizializza `res1_1` come SCALARE ComplexF64 ad ogni iterazione del loop esterno.
        res1_1_scalar = zero(ComplexF64)
        xp = h1 * rx[a1] + h2

        for a2 in 1:my
            yp = h3 * ry[a2] + h4
            # delta_x = (x - xp)
            # delta_y = (y - yp)
            # delta_z = (z - zp)
            # R = sqrt(delta_x^2 + delta_y^2 + delta_z^2)
			delta_vec = SVector(x - xp, y - yp, z - zp)
			R = norm(delta_vec)

            # Modifica per evitare allocazioni in `mul_fast` con molti argomenti:
            # abbiamo aggiunto parentesi per forzare la moltiplicazione in coppie.
            #G_1 = delta_z * ((1 / R^3 + 1im * beta / R^2) * exp(-1im * beta * R))
            G_1 = delta_vec[3] * (1 / R^3 + 1im * beta / R^2) * exp(-1im * beta * R)

            # Somma direttamente allo scalare, non un'operazione su array
            res1_1_scalar += wy[a2] * G_1
        end
        # Somma direttamente allo scalare, non un'operazione su array
        res2_1_scalar += wx[a1] * h3 * res1_1_scalar
    end

    res2_1_scalar = h1 * res2_1_scalar

    # Ritorna il singolo valore scalare
    return res2_1_scalar
end

function compute_hcx_xy(barra::SVector{12, Float64}, centriOss::SVector{3, Float64}, beta::ComplexF64, rx::Vector{Float64}, wx::Vector{Float64}, ry::Vector{Float64}, wy::Vector{Float64})
    # numCentri non è più necessario, poiché il risultato è un singolo scalare.
    # numCentri = size(centriOss, 1) # Rimosso

    # Pre-calcolo delle coordinate min/max e degli h.
    # Queste operazioni sono efficienti con StaticArrays e non allocano.
    x1 = min(barra[1], barra[4], barra[7], barra[10])
    x2 = max(barra[1], barra[4], barra[7], barra[10])
    y1 = min(barra[2], barra[5], barra[8], barra[11])
    y2 = max(barra[2], barra[5], barra[8], barra[11])
    h1 = (x2 - x1) / 2
    h2 = (x2 + x1) / 2
    h3 = (y2 - y1) / 2
    h4 = (y2 + y1) / 2

    nx = length(wx)
    my = length(wy)

    # Inizializziamo `res2_1` come SCALARE ComplexF64, non un array.
    res2_1_scalar = zero(ComplexF64)

    # Estraiamo i componenti di centriOss e barra una volta.
    x = centriOss[1]
    y = centriOss[2]
    z = centriOss[3]
    zp = barra[3]

    @inbounds for a1 in 1:nx
        # Inizializziamo `res1_1` come SCALARE ComplexF64 ad ogni iterazione del loop esterno.
        res1_1_scalar = zero(ComplexF64)
        xp = h1 * rx[a1] + h2

        for a2 in 1:my
            yp = h3 * ry[a2] + h4
            # delta_x = (x - xp)
            # delta_y = (y - yp)
            # delta_z = (z - zp)
            # R = sqrt(delta_x^2 + delta_y^2 + delta_z^2)
			delta_vec = SVector(x - xp, y - yp, z - zp)
			R = norm(delta_vec)

            # Modifica per evitare allocazioni in `mul_fast` con molti argomenti:
            # abbiamo aggiunto parentesi per forzare la moltiplicazione in coppie.
            # La differenza rispetto a compute_hcz_xy è che qui usiamo `delta_x` invece di `delta_z`.
            #G_1 = delta_x * ((1 / R^3 + 1im * beta / R^2) * exp(-1im * beta * R))
            G_1 = delta_vec[1] * (1 / R^3 + 1im * beta / R^2) * exp(-1im * beta * R)

            # Sommiamo direttamente allo scalare, non un'operazione su array.
            res1_1_scalar += wy[a2] * G_1
        end
        # Sommiamo direttamente allo scalare, non un'operazione su array.
        res2_1_scalar += wx[a1] * h3 * res1_1_scalar
    end

    res2_1_scalar = h1 * res2_1_scalar

    # La funzione ora ritorna direttamente il singolo valore scalare.
    return res2_1_scalar
end


# using LinearAlgebra

# function compute_Ec_Gauss(barre, normale, centriOss, ordine, beta, simulation_id, chan)
# 	num_barre = size(barre, 1)
# 	num_centri_oss = size(centriOss, 1)
# 	hc = zeros(ComplexF64, num_barre, 3, num_centri_oss)

# 	for i in eachindex(normale)
# 		normale[i] = convert(Vector{Float64}, normale[i])
# 	end

# 	block_size = 100
# 	normale = hcat(normale...)
#     rx, wx = qrule(ordine)
# 	ry, wy = qrule(ordine)
# 	for m_block in 1:block_size:num_barre
# 		m_end = min(m_block + block_size - 1, num_barre)
# 		Base.Threads.@threads for cont in m_block:m_end
# 			if abs(normale[cont, 1]) > 1e-10
# 				perm = [3, 2, 1, 6, 5, 4, 9, 8, 7, 12, 11, 10]
# 				perm2 = [2, 3, 1, 5, 6, 4, 8, 9, 7, 11, 12, 10]
# 				for cc in 1:num_centri_oss
# 					hc[cont, 1, cc] = compute_hcz_xy(barre[cont, perm], centriOss[cc, [3, 2, 1]], ordine, beta, rx, wx, ry, wy)
# 					hc[cont, 2, cc] = compute_hcx_xy(barre[cont, perm2], centriOss[cc, [2, 3, 1]], ordine, beta, rx, wx, ry, wy)
# 					hc[cont, 3, cc] = compute_hcx_xy(barre[cont, perm], centriOss[cc, [3, 2, 1]], ordine, beta, rx, wx, ry, wy)
# 					if cc % 100 == 0
# 						yield()  # Permette alla task dell'heartbeat di essere schedulata
# 					end
# 				end
# 			elseif abs(normale[cont, 2]) > 1e-10
# 				perm = [1, 3, 2, 4, 6, 5, 7, 9, 8, 10, 12, 11]
# 				perm2 = [3, 1, 2, 6, 4, 5, 9, 7, 8, 12, 10, 11]
# 				for cc in 1:num_centri_oss
# 					hc[cont, 1, cc] = compute_hcx_xy(barre[cont, perm], centriOss[cc, [1, 3, 2]], ordine, beta, rx, wx, ry, wy)
# 					hc[cont, 2, cc] = compute_hcz_xy(barre[cont, perm], centriOss[cc, [1, 3, 2]], ordine, beta, rx, wx, ry, wy)
# 					hc[cont, 3, cc] = compute_hcx_xy(barre[cont, perm2], centriOss[cc, [3, 1, 2]], ordine, beta, rx, wx, ry, wy)
# 					if cc % 100 == 0
# 						yield()  # Permette alla task dell'heartbeat di essere schedulata
# 					end
# 				end
# 			else
# 				perm = [2, 1, 3, 5, 4, 6, 8, 7, 9, 11, 10, 12]
# 				for cc in 1:num_centri_oss
# 					hc[cont, 1, cc] = compute_hcx_xy(barre[cont, :], centriOss[cc, :], ordine, beta, rx, wx, ry, wy)
# 					hc[cont, 2, cc] = compute_hcx_xy(barre[cont, perm], centriOss[cc, [2, 1, 3]], ordine, beta, rx, wx, ry, wy)
# 					hc[cont, 3, cc] = compute_hcz_xy(barre[cont, :], centriOss[cc, :], ordine, beta, rx, wx, ry, wy)
# 					if cc % 100 == 0
# 						yield()  # Permette alla task dell'heartbeat di essere schedulata
# 					end
# 				end
# 			end
# 		end
# 		sleep(0)
# 		println("block Ec : ", round(m_end / block_size), " / ", round(num_barre / block_size))
# 		if is_stop_requested(simulation_id)
# 			println("Simulazione $(simulation_id) interrotta per richiesta stop.")
# 			return nothing # O un altro valore che indica interruzione
# 		end
# 	end

# 	# Base.Threads.@threads for cont in 1:num_barre

# 	# end
# 	if is_stop_requested(simulation_id)
# 		println("Simulazione $(simulation_id) interrotta per richiesta stop.")
# 		return nothing # O un altro valore che indica interruzione
# 	else
# 		return hc
# 	end
# end

# function compute_hcz_xy(barra, centriOss, ordine, beta, rx, wx, ry, wy)
# 	numCentri = size(centriOss, 1) # Julia's size on a vector gives a tuple, so we don't need the second dimension
# 	x1 = minimum(barra[[1, 4, 7, 10]])
# 	x2 = maximum(barra[[1, 4, 7, 10]])
# 	y1 = minimum(barra[[2, 5, 8, 11]])
# 	y2 = maximum(barra[[2, 5, 8, 11]])
# 	h1 = (x2 - x1) / 2
# 	h2 = (x2 + x1) / 2
# 	h3 = (y2 - y1) / 2
# 	h4 = (y2 + y1) / 2
# 	nx = length(wx)
# 	my = length(wy)
# 	res2_1 = zeros(ComplexF64, numCentri)
# 	x = centriOss[1] # Accessing single elements directly
# 	y = centriOss[2]
# 	z = centriOss[3]
# 	zp = barra[3]
# 	for a1 in 1:nx
# 		res1_1 = zeros(ComplexF64, numCentri)
# 		xp = h1 * rx[a1] + h2
# 		for a2 in 1:my
# 			yp = h3 * ry[a2] + h4
# 			delta_x = (x - xp)
# 			delta_y = (y - yp)
# 			delta_z = (z - zp)
# 			R = sqrt(delta_x^2 + delta_y^2 + delta_z^2)
# 			G_1 = delta_z * (1 / R^3 + 1im * beta / R^2) * exp(-1im * beta * R)
# 			res1_1 .+= wy[a2] * G_1
# 		end
# 		res2_1 .+= wx[a1] * h3 * res1_1
# 	end
# 	res2_1 = h1 * res2_1
# 	return res2_1[1]
# end

# function compute_hcx_xy(barra, centriOss, ordine, beta, rx, wx, ry, wy)
# 	numCentri = size(centriOss, 1)
# 	x1 = minimum(barra[[1, 4, 7, 10]])
# 	x2 = maximum(barra[[1, 4, 7, 10]])
# 	y1 = minimum(barra[[2, 5, 8, 11]])
# 	y2 = maximum(barra[[2, 5, 8, 11]])
# 	h1 = (x2 - x1) / 2
# 	h2 = (x2 + x1) / 2
# 	h3 = (y2 - y1) / 2
# 	h4 = (y2 + y1) / 2
# 	nx = length(wx)
# 	my = length(wy)
# 	res2_1 = zeros(ComplexF64, numCentri)
# 	x = centriOss[1]
# 	y = centriOss[2]
# 	z = centriOss[3]
# 	zp = barra[3]
# 	for a1 in 1:nx
# 		res1_1 = zeros(ComplexF64, numCentri)
# 		xp = h1 * rx[a1] + h2
# 		for a2 in 1:my
# 			yp = h3 * ry[a2] + h4
# 			delta_x = (x - xp)
# 			delta_y = (y - yp)
# 			delta_z = (z - zp)
# 			R = sqrt(delta_x^2 + delta_y^2 + delta_z^2)
# 			G_1 = delta_x * (1 / R^3 + 1im * beta / R^2) * exp(-1im * beta * R)
# 			res1_1 .+= wy[a2] * G_1
# 		end
# 		res2_1 .+= wx[a1] * h3 * res1_1
# 	end
# 	res2_1 = h1 * res2_1
# 	return res2_1[1]
# end

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