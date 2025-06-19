using LinearAlgebra
using StaticArrays # Essenziale per prestazioni elevate con array di piccole dimensioni
using Base.Threads
using Printf

# function compute_Ar_Gauss(barre, centriOss, ordine, beta, simulation_id, chan)
# 	numCentri = size(centriOss, 1)
# 	numBarre = size(barre, 1)

# 	ha = zeros(ComplexF64, numBarre, numCentri)
#     rootkx, wekx = qrule(ordine)
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

# function compute_ha(
#     xo::Float64, x_vect_bar::Vector{Float64}, yo::Float64, y_vect_bar::Vector{Float64}, 
#     zo::Float64, z_vect_bar::Vector{Float64}, 
#     ordine::Int, beta,
#     rootkx::Vector{Float64}, wekx::Vector{Float64},
#     rootky::Vector{Float64}, weky::Vector{Float64},
#     rootkz::Vector{Float64}, wekz::Vector{Float64})
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
# 	@inbounds for a1 in 1:nlkx
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
function compute_Ar_Gauss(barre::Transpose{Float64, Matrix{Float64}},
                                    centriOss::Matrix{Float64}, ordine::Int, beta::ComplexF64,
                                    simulation_id=nothing, chan=nothing)
    numCentri = size(centriOss, 1)
    N = size(barre, 1)

    # Determina il numero di segmenti per la parallelizzazione
    # La logica è simile a quella di MATLAB, ma Julia usa `Int` per gli indici
    num_punti = ceil(Int, N / 100)
    vect_ind = ceil.(Int, LinRange(1, N, num_punti))

    # Assicura che vect_ind abbia almeno due elementi e che l'ultimo indice sia N
    if length(vect_ind) < 2
        vect_ind = [1, N]
    end
    if vect_ind[end] > N
        vect_ind[end] = N
    end

    Nq = length(vect_ind) - 1
    barrec = Vector{Matrix{Float64}}(undef, Nq) # Pre-allocazione per efficienza

    # Suddivisione delle 'barre' in segmenti
    for pos = 1:Nq
        # Aggiustamento per l'ultimo segmento
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
        println("block Ar : ", c, " / ", Nq)
    end

    # Inizializzazione della matrice dei risultati con numeri complessi
    ha = zeros(ComplexF64, N, numCentri)

    # Ri-assemblaggio dei risultati
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
    
    # Initialize 'ha' as a matrix of complex numbers
    ha = zeros(ComplexF64, size(barre, 1), numCentri)
    
    # Loop through each 'barra'
    for cont = 1:size(barre, 1)
        barra = barre[cont, :]
        
        # Extracting coordinates. Julia's array indexing makes this concise.
        # Note: MATLAB allows direct indexing like barra(1), barra(4), etc.
        # In Julia, you explicitly select elements: barra[1], barra[4].
        xb = [barra[1], barra[4], barra[7], barra[10], barra[13], barra[16], barra[19], barra[22]]
        yb = [barra[2], barra[5], barra[8], barra[11], barra[14], barra[17], barra[20], barra[23]]
        zb = [barra[3], barra[6], barra[9], barra[12], barra[15], barra[18], barra[21], barra[24]]
        
        # Finding min and max for each coordinate
        x_bar = [minimum(xb), maximum(xb)]
        y_bar = [minimum(yb), maximum(yb)]
        z_bar = [minimum(zb), maximum(zb)]

        x1, x2 = x_bar[1], x_bar[2]
        y1, y2 = y_bar[1], y_bar[2]
        z1, z2 = z_bar[1], z_bar[2]
        
        # Define 'barra' (vertices of the 3D box)
        barra = [x1, y1, z1, x2, y1, z1, x1, y2, z1, x2, y2, z1, x1, y1, z2, x2, y1, z2, x1, y2, z2, x2, y2, z2]
        
        # Extracting groups of coordinates for two quadrilaterals
        xi1 = [barra[1], barra[4], barra[7], barra[10]]
        yi1 = [barra[2], barra[5], barra[8], barra[11]]
        zi1 = [barra[3], barra[6], barra[9], barra[12]]
        xi2 = [barra[13], barra[16], barra[19], barra[22]]
        yi2 = [barra[14], barra[17], barra[20], barra[23]]
        zi2 = [barra[15], barra[18], barra[21], barra[24]]
        
        # Initialize 'ri' as a matrix of complex numbers (vectors pointing to vertices)
        # Although MATLAB initialized it with `1j*zeros`, the assignments immediately overwrite this
        # with real values, then it's used in complex calculations later.
        # In Julia, it's often clearer to start with the correct type.
        # However, since `ri` values are directly assigned from `xi1`, `yi1`, etc. which are real,
        # we can just make `ri` a Matrix{Float64} and let Julia handle type promotion in subsequent calculations.
        ri = Matrix{Float64}(undef, 8, 3) 
        ri[1, :] = [xi1[1], yi1[1], zi1[1]]
        ri[2, :] = [xi1[2], yi1[2], zi1[2]]
        ri[3, :] = [xi1[3], yi1[3], zi1[3]]
        ri[4, :] = [xi1[4], yi1[4], zi1[4]]
        ri[5, :] = [xi2[1], yi2[1], zi2[1]]
        ri[6, :] = [xi2[2], yi2[2], zi2[2]]
        ri[7, :] = [xi2[3], yi2[3], zi2[3]]
        ri[8, :] = [xi2[4], yi2[4], zi2[4]]
        
        # New approach: Calculate coefficients for the tensor product basis functions
        # Julia's broadcasting and element-wise operations are powerful here.
        # Note: `sum(ri, 1)` in MATLAB sums along dimension 1 (columns).
        # In Julia, `sum(ri, dims=1)` sums along the first dimension, returning a 1x3 matrix.
        # To get a 1D vector like MATLAB's `sum(ri,1)` often behaves, you might do `vec(sum(ri, dims=1))`.
        # However, for these vector math operations, keeping them as row vectors (1x3 matrices) is fine.
        # Or, if these are meant to be 3-element vectors, you can explicitly reshape.
        # For now, I'll assume they should be 1x3 matrices to match MATLAB's `sum(ri,1)` behavior.
        
        # Update: Based on subsequent usage (`rmi + rai*rootkx(a1)`), these are indeed 3-element vectors.
        # So `vec(sum(ri, dims=1))` or simply ensuring the sum results in a 1D array.
        # A cleaner way in Julia for summing rows to get a single vector: `sum(eachrow(ri))`.
        rmi = 0.125 * sum(eachrow(ri)) 
        rai = 0.125 * (-ri[1,:] .+ ri[2,:] .+ ri[4,:] .- ri[3,:] .- ri[5,:] .+ ri[6,:] .+ ri[8,:] .- ri[7,:])
        rbi = 0.125 * (-ri[1,:] .- ri[2,:] .+ ri[4,:] .+ ri[3,:] .- ri[5,:] .- ri[6,:] .+ ri[8,:] .+ ri[7,:])
        rci = 0.125 * (-ri[1,:] .- ri[2,:] .- ri[4,:] .- ri[3,:] .+ ri[5,:] .+ ri[6,:] .+ ri[8,:] .+ ri[7,:])
        rabi = 0.125 * (ri[1,:] .- ri[2,:] .+ ri[4,:] .- ri[3,:] .+ ri[5,:] .- ri[6,:] .+ ri[8,:] .- ri[7,:])
        rbci = 0.125 * (ri[1,:] .+ ri[2,:] .- ri[4,:] .- ri[3,:] .- ri[5,:] .- ri[6,:] .+ ri[8,:] .+ ri[7,:])
        raci = 0.125 * (ri[1,:] .- ri[2,:] .- ri[4,:] .+ ri[3,:] .- ri[5,:] .+ ri[6,:] .+ ri[8,:] .- ri[7,:])
        rabci = 0.125 * (-ri[1,:] .+ ri[2,:] .- ri[4,:] .+ ri[3,:] .+ ri[5,:] .- ri[6,:] .+ ri[8,:] .- ri[7,:])
        
        # Loop through each observation center
        for cc = 1:numCentri
            x_o = centriOss[cc, 1]
            y_o = centriOss[cc, 2]
            z_o = centriOss[cc, 3]
            
            # Call the helper function compute_ha
            rai_s  = SVector{3}(rai)
            rbi_s  = SVector{3}(rbi)
            rci_s  = SVector{3}(rci)
            rmi_s  = SVector{3}(rmi)
            rabi_s = SVector{3}(rabi)
            rbci_s = SVector{3}(rbci)
            raci_s = SVector{3}(raci)
            rabci_s= SVector{3}(rabci)
            ha[cont, cc] = compute_ha_fast(x_o, y_o, z_o, rootkx, wekx, rootky, weky, rootkz, wekz, beta, rmi_s, rai_s, rbi_s, rci_s, rabi_s, rbci_s, raci_s, rabci_s)
        end
    end
    
    return ha
end



@inline norm3(v::SVector{3,Float64}) = sqrt(muladd(v[1],v[1],
                                          muladd(v[2],v[2],
                                                 v[3]*v[3])))

function compute_ha_fast(xo::Float64, yo::Float64, zo::Float64,
                         rootkx::Vector{Float64}, wekx::Vector{Float64},
                         rootky::Vector{Float64}, weky::Vector{Float64},
                         rootkz::Vector{Float64}, wekz::Vector{Float64},
                         beta::ComplexF64,
                         rmi::SVector{3,Float64}, rai::SVector{3,Float64},
                         rbi::SVector{3,Float64}, rci::SVector{3,Float64},
                         rabi::SVector{3,Float64}, rbci::SVector{3,Float64},
                         raci::SVector{3,Float64}, rabci::SVector{3,Float64})

    nlkx, nlky, nlkz = length(wekx), length(weky), length(wekz)
    sum_a1 = 0.0 + 0.0im

    @inbounds for a1 in 1:nlkx
        dkx = rootkx[a1]
        wx  = wekx[a1]
        sum_b1 = 0.0 + 0.0im

        @inbounds for b1 in 1:nlky
            dky = rootky[b1]
            wy  = weky[b1]
            sum_c1 = 0.0 + 0.0im

            # Questi pezzi dipendono solo da dkx/dky → si possono precalcolare
            drai_xy = rai           .+ rabi .* dky
            drbi_x  = rbi           .+ rabi .* dkx
            drci_xy = rci           .+ rbci .* dky    # parte comune usata sotto
            # r1_xy sono i termini fino a dkx,dky (manca dkz)
            r1_xy   = rmi .+ rai .* dkx .+ rbi .* dky .+ rabi .* dkx .* dky

            @inbounds for c1 in 1:nlkz
                dkz = rootkz[c1]
                wz  = wekz[c1]

                # --- vettori 3-componenti (tutto su registri) ---
                dra  = drai_xy         .+ raci .* dkz .+ rabci .* dky .* dkz
                drb  = drbi_x          .+ rbci .* dkz .+ rabci .* dkx .* dkz
                drc  = drci_xy         .+ raci .* dkx .* dkz .+ rbci .* dkz

                dra_m = norm3(dra)
                drb_m = norm3(drb)
                drc_m = norm3(drc)

                r1 = r1_xy .+ rci .* dkz .+ raci .* dkx .* dkz .+
                     rbci .* dky .* dkz .+ rabci .* dkx .* dky .* dkz

                Δx = xo - r1[1]
                Δy = yo - r1[2]
                Δz = zo - r1[3]
                R  = sqrt(muladd(Δx,Δx, muladd(Δy,Δy, Δz*Δz)))

                G = R > eps(Float64) ? (inv(R) * exp(-im*beta*R)) : 0.0 + 0.0im

                sum_c1 += wz * (dra_m * drb_m * drc_m) * G
            end
            sum_b1 += wy * sum_c1
        end
        sum_a1 += wx * sum_b1
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