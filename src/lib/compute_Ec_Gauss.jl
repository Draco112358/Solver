using LinearAlgebra

function compute_Ec_Gauss(barre, normale, centriOss, ordine, beta, simulation_id, chan)
	num_barre = size(barre, 1)
	num_centri_oss = size(centriOss, 1)
	hc = zeros(ComplexF64, num_barre, 3, num_centri_oss)

	for i in eachindex(normale)
		normale[i] = convert(Vector{Float64}, normale[i])
	end

	block_size = 100
	normale = hcat(normale...)
	for m_block in 1:block_size:num_barre
		m_end = min(m_block + block_size - 1, num_barre)
		Base.Threads.@threads for cont in m_block:m_end
			if abs(normale[cont, 1]) > 1e-10
				perm = [3, 2, 1, 6, 5, 4, 9, 8, 7, 12, 11, 10]
				perm2 = [2, 3, 1, 5, 6, 4, 8, 9, 7, 11, 12, 10]
				for cc in 1:num_centri_oss
					hc[cont, 1, cc] = compute_hcz_xy(barre[cont, perm], centriOss[cc, [3, 2, 1]], ordine, beta)
					hc[cont, 2, cc] = compute_hcx_xy(barre[cont, perm2], centriOss[cc, [2, 3, 1]], ordine, beta)
					hc[cont, 3, cc] = compute_hcx_xy(barre[cont, perm], centriOss[cc, [3, 2, 1]], ordine, beta)
					if cc % 100 == 0
						yield()  # Permette alla task dell'heartbeat di essere schedulata
					end
				end
			elseif abs(normale[cont, 2]) > 1e-10
				perm = [1, 3, 2, 4, 6, 5, 7, 9, 8, 10, 12, 11]
				perm2 = [3, 1, 2, 6, 4, 5, 9, 7, 8, 12, 10, 11]
				for cc in 1:num_centri_oss
					hc[cont, 1, cc] = compute_hcx_xy(barre[cont, perm], centriOss[cc, [1, 3, 2]], ordine, beta)
					hc[cont, 2, cc] = compute_hcz_xy(barre[cont, perm], centriOss[cc, [1, 3, 2]], ordine, beta)
					hc[cont, 3, cc] = compute_hcx_xy(barre[cont, perm2], centriOss[cc, [3, 1, 2]], ordine, beta)
					if cc % 100 == 0
						yield()  # Permette alla task dell'heartbeat di essere schedulata
					end
				end
			else
				perm = [2, 1, 3, 5, 4, 6, 8, 7, 9, 11, 10, 12]
				for cc in 1:num_centri_oss
					hc[cont, 1, cc] = compute_hcx_xy(barre[cont, :], centriOss[cc, :], ordine, beta)
					hc[cont, 2, cc] = compute_hcx_xy(barre[cont, perm], centriOss[cc, [2, 1, 3]], ordine, beta)
					hc[cont, 3, cc] = compute_hcz_xy(barre[cont, :], centriOss[cc, :], ordine, beta)
					if cc % 100 == 0
						yield()  # Permette alla task dell'heartbeat di essere schedulata
					end
				end
			end
		end
		sleep(0)
		println("block Ec : ", round(m_end / block_size), " / ", round(num_barre / block_size))
	end

	# Base.Threads.@threads for cont in 1:num_barre

	# end
	if is_stop_requested(simulation_id)
		println("Simulazione $(simulation_id) interrotta per richiesta stop.")
		return nothing # O un altro valore che indica interruzione
	else
		return hc
	end
end

function compute_hcz_xy(barra, centriOss, ordine, beta)
	numCentri = size(centriOss, 1) # Julia's size on a vector gives a tuple, so we don't need the second dimension
	x1 = minimum(barra[[1, 4, 7, 10]])
	x2 = maximum(barra[[1, 4, 7, 10]])
	y1 = minimum(barra[[2, 5, 8, 11]])
	y2 = maximum(barra[[2, 5, 8, 11]])
	rx, wx = qrule(ordine)
	ry, wy = qrule(ordine)
	h1 = (x2 - x1) / 2
	h2 = (x2 + x1) / 2
	h3 = (y2 - y1) / 2
	h4 = (y2 + y1) / 2
	nx = length(wx)
	my = length(wy)
	res2_1 = zeros(ComplexF64, numCentri)
	x = centriOss[1] # Accessing single elements directly
	y = centriOss[2]
	z = centriOss[3]
	zp = barra[3]
	for a1 in 1:nx
		res1_1 = zeros(ComplexF64, numCentri)
		xp = h1 * rx[a1] + h2
		for a2 in 1:my
			yp = h3 * ry[a2] + h4
			delta_x = (x - xp)
			delta_y = (y - yp)
			delta_z = (z - zp)
			R = sqrt(delta_x^2 + delta_y^2 + delta_z^2)
			G_1 = delta_z * (1 / R^3 + 1im * beta / R^2) * exp(-1im * beta * R)
			res1_1 .+= wy[a2] * G_1
		end
		res2_1 .+= wx[a1] * h3 * res1_1
	end
	res2_1 = h1 * res2_1
	return res2_1[1]
end

function compute_hcx_xy(barra, centriOss, ordine, beta)
	numCentri = size(centriOss, 1)
	x1 = minimum(barra[[1, 4, 7, 10]])
	x2 = maximum(barra[[1, 4, 7, 10]])
	y1 = minimum(barra[[2, 5, 8, 11]])
	y2 = maximum(barra[[2, 5, 8, 11]])
	rx, wx = qrule(ordine)
	ry, wy = qrule(ordine)
	h1 = (x2 - x1) / 2
	h2 = (x2 + x1) / 2
	h3 = (y2 - y1) / 2
	h4 = (y2 + y1) / 2
	nx = length(wx)
	my = length(wy)
	res2_1 = zeros(ComplexF64, numCentri)
	x = centriOss[1]
	y = centriOss[2]
	z = centriOss[3]
	zp = barra[3]
	for a1 in 1:nx
		res1_1 = zeros(ComplexF64, numCentri)
		xp = h1 * rx[a1] + h2
		for a2 in 1:my
			yp = h3 * ry[a2] + h4
			delta_x = (x - xp)
			delta_y = (y - yp)
			delta_z = (z - zp)
			R = sqrt(delta_x^2 + delta_y^2 + delta_z^2)
			G_1 = delta_x * (1 / R^3 + 1im * beta / R^2) * exp(-1im * beta * R)
			res1_1 .+= wy[a2] * G_1
		end
		res2_1 .+= wx[a1] * h3 * res1_1
	end
	res2_1 = h1 * res2_1
	return res2_1[1]
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
