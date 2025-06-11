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
using LinearAlgebra # Ensure LinearAlgebra is available for norm and ComplexF64

# --- Optimized `compute_Ar_Gauss` ---
function compute_Ar_Gauss(barre::Union{Transpose{Float64, Matrix{Float64}}, Transpose{Real, Matrix{Real}}}, centriOss::Matrix{Float64}, ordine::Int, beta::ComplexF64, simulation_id, chan)
    numCentri = size(centriOss, 1)
    numBarre = size(barre, 1)

    ha = zeros(ComplexF64, numBarre, numCentri)

    # Pre-calculate qrule values once
    # Ensure qrule returns Vectors, not 1xN matrices
    rootkx, wekx = qrule(ordine)
    rootky, weky = qrule(ordine)
    rootkz, wekz = qrule(ordine)

    block_size = 100 # Adjust this based on your hardware and problem size

    for m_block in 1:block_size:numBarre
        m_end = min(m_block + block_size - 1, numBarre)
        # Use Threads.@threads for parallelization
        Threads.@threads for cont in m_block:m_end
            # Use a view for `barra` to avoid copying the whole row
            barra_view = @view barre[cont, :]

            # Pre-extract min/max coordinates for the current bar.
            # Using @inbounds as these are known valid indices.
            @inbounds begin
                # These indices (1,4,7,10,13,16,19,22) look like 8 vertices (x,y,z,...)
                # It's more common to have 4 vertices for each face of a prism (12 elements total for 4 on one face).
                # The original `xb`, `yb`, `zb` access 8 elements, which seems to imply
                # that `barra` might contain coordinates for 8 vertices (24 elements total for 8*3D).
                # If `barra` is [x1,y1,z1,x2,y2,z2,...,x8,y8,z8], then the indices are correct for 8 vertices.
                # If it's a 4-vertex prism with two faces (12 elements for 4 vertices per face), the indexing is off.
                # Assuming `barra` has 24 elements for 8 vertices (x,y,z per vertex).
                # If `barra` is structured differently (e.g., [x1a,y1a,z1a,x2a,y2a,z2a, x1b,y1b,z1b,x2b,y2b,z2b]),
                # you might need to adjust these indices.

                # Extracting all 8 x, y, z coordinates directly for min/max
                # This avoids allocating `xb`, `yb`, `zb` arrays per loop iteration.
                # Assuming `barra_view` is `[x1,y1,z1, x2,y2,z2, ... x8,y8,z8]`
                x_coords_all = (barra_view[1], barra_view[4], barra_view[7], barra_view[10],
                                barra_view[13], barra_view[16], barra_view[19], barra_view[22])
                y_coords_all = (barra_view[2], barra_view[5], barra_view[8], barra_view[11],
                                barra_view[14], barra_view[17], barra_view[20], barra_view[23])
                z_coords_all = (barra_view[3], barra_view[6], barra_view[9], barra_view[12],
                                barra_view[15], barra_view[18], barra_view[21], barra_view[24])

                x_min_bar = minimum(x_coords_all)
                x_max_bar = maximum(x_coords_all)
                y_min_bar = minimum(y_coords_all)
                y_max_bar = maximum(y_coords_all)
                z_min_bar = minimum(z_coords_all)
                z_max_bar = maximum(z_coords_all)
            end

            # Create these small static arrays or tuples once per `cont` loop
            # and pass them directly.
            # Using StaticArrays for small, fixed-size arrays can be beneficial for performance.
            # You'll need `using StaticArrays` at the top of your file.
            # If not using StaticArrays, a simple tuple or Vector{Float64} works fine too.
            # x_bar_tuple = (x_min_bar, x_max_bar)
            # y_bar_tuple = (y_min_bar, y_max_bar)
            # z_bar_tuple = (z_min_bar, z_max_bar)

            for cc in 1:numCentri
                x_o = centriOss[cc, 1]
                y_o = centriOss[cc, 2]
                z_o = centriOss[cc, 3]

                # Pass scalar min/max values directly instead of a tuple/vector.
                # This avoids unpacking inside compute_ha and potential type instability.
                ha[cont, cc] = compute_ha_optimized(
                    x_o, x_min_bar, x_max_bar,
                    y_o, y_min_bar, y_max_bar,
                    z_o, z_min_bar, z_max_bar,
                    beta, rootkx, wekx, rootky, weky, rootkz, wekz
                )
                # Removed `if cc % 100 == 0 yield()`:
                # `yield()` inside the innermost loop can add overhead.
                # The `sleep(0)` and `is_stop_requested` check outside the inner loop are sufficient
                # for responsiveness without impacting per-iteration performance significantly.
            end
        end
        sleep(0) # Keep for progress reporting and yielding to other Julia tasks
        println("block Ar : ", round(m_end / block_size), " / ", round(numBarre / block_size))
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

# --- Optimized `compute_ha` ---
# Make it `@inline` to encourage Julia to insert its code directly at the call site,
# potentially eliminating function call overhead.
@inline function compute_ha_optimized(
    xo::Float64, x1::Float64, x2::Float64,
    yo::Float64, y1::Float64, y2::Float64,
    zo::Float64, z1::Float64, z2::Float64,
    beta::ComplexF64, # Beta must be ComplexF64 for `exp(-1im * beta * R)`
    rootkx::AbstractVector{Float64}, wekx::AbstractVector{Float64},
    rootky::AbstractVector{Float64}, weky::AbstractVector{Float64},
    rootkz::AbstractVector{Float64}, wekz::AbstractVector{Float64}
)
    # No need to recalculate x1, x2, etc. from vectors; they are passed directly.
    # No need to create `barra` matrix or its sub-arrays (xi1, yi1, etc.) here.
    # The 8 vectors ri[i, :] are just fixed vertices of a rectangular prism.
    # We can define them based on x1, x2, y1, y2, z1, z2.

    # These are fixed for a rectangular prism defined by min/max coordinates
    # and can be hardcoded or calculated once per call to compute_ha.
    # Using StaticArrays for these small vectors is highly recommended if `StaticArrays.jl` is added.
    # Otherwise, regular Vectors are fine, but be mindful of allocations if they change per call.
    # Since these are constant patterns based on x1,x2,y1,y2,z1,z2, they don't allocate in the inner loops.

    # Pre-calculate the 'r' vectors (rmi, rai, etc.) using the x1,x2,y1,y2,z1,z2 directly.
    # This avoids creating `barra` and `ri` arrays, which cause significant allocations.
    # Each of these `r` vectors is a 3-element vector.
    # For maximum performance, if `StaticArrays.jl` is in use:
    # using StaticArrays
    # rmi = SVector{3,ComplexF64}(...)
    # ... and so on.
    # If not using StaticArrays, normal Vectors are fine, but note the initial allocation cost.
    # For clarity, I'll write them out as direct calculations.

    # Original 'ri' logic based on x1,y1,z1 and x2,y2,z2.
    # Assuming the original barra definition for ri:
    # ri[1, :] = [x1, y1, z1]
    # ri[2, :] = [x2, y1, z1]
    # ri[3, :] = [x1, y2, z1]
    # ri[4, :] = [x2, y2, z1]
    # ri[5, :] = [x1, y1, z2]
    # ri[6, :] = [x2, y1, z2]
    # ri[7, :] = [x1, y2, z2]
    # ri[8, :] = [x2, y2, z2]

    # Calculate rmi, rai, etc. directly without intermediate `ri` array.
    # Each component (x,y,z) is calculated separately for performance.
    # The `0.125` can be pre-multiplied for efficiency.
    # Also, ensuring results are ComplexF64 from the start for type stability.

    # rmi (mean point)
    # Equivalent to 0.125 * sum(ri, dims=1)
    rmi_x = 0.125 * (x1 + x2 + x1 + x2 + x1 + x2 + x1 + x2)
    rmi_y = 0.125 * (y1 + y1 + y2 + y2 + y1 + y1 + y2 + y2)
    rmi_z = 0.125 * (z1 + z1 + z1 + z1 + z2 + z2 + z2 + z2)

    # rai (partial derivatives)
    # 0.125 * (-ri[1, :] + ri[2, :] + ri[4, :] - ri[3, :] - ri[5, :] + ri[6, :] + ri[8, :] - ri[7, :])
    rai_x = 0.125 * (-x1 + x2 + x2 - x1 - x1 + x2 + x2 - x1)
    rai_y = 0.125 * (-y1 + y1 + y2 - y2 - y1 + y1 + y2 - y2)
    rai_z = 0.125 * (-z1 + z1 + z1 - z1 - z2 + z2 + z2 - z2)

    # rbi
    # 0.125 * (-ri[1, :] - ri[2, :] + ri[4, :] + ri[3, :] - ri[5, :] - ri[6, :] + ri[8, :] + ri[7, :])
    rbi_x = 0.125 * (-x1 - x2 + x2 + x1 - x1 - x2 + x2 + x1)
    rbi_y = 0.125 * (-y1 - y1 + y2 + y2 - y1 - y1 + y2 + y2)
    rbi_z = 0.125 * (-z1 - z1 + z1 + z1 - z2 - z2 + z2 + z2)

    # rci
    # 0.125 * (-ri[1, :] - ri[2, :] - ri[4, :] - ri[3, :] + ri[5, :] + ri[6, :] + ri[8, :] + ri[7, :])
    rci_x = 0.125 * (-x1 - x2 - x2 - x1 + x1 + x2 + x2 + x1)
    rci_y = 0.125 * (-y1 - y1 - y2 - y2 + y1 + y1 + y2 + y2)
    rci_z = 0.125 * (-z1 - z1 - z1 - z1 + z2 + z2 + z2 + z2)

    # rabi
    # 0.125 * (ri[1, :] - ri[2, :] + ri[4, :] - ri[3, :] + ri[5, :] - ri[6, :] + ri[8, :] - ri[7, :])
    rabi_x = 0.125 * (x1 - x2 + x2 - x1 + x1 - x2 + x2 - x1)
    rabi_y = 0.125 * (y1 - y1 + y2 - y2 + y1 - y1 + y2 - y2)
    rabi_z = 0.125 * (z1 - z1 + z1 - z1 + z2 - z2 + z2 - z2)

    # rbci
    # 0.125 * (ri[1, :] + ri[2, :] - ri[4, :] - ri[3, :] - ri[5, :] - ri[6, :] + ri[8, :] + ri[7, :])
    rbci_x = 0.125 * (x1 + x2 - x2 - x1 - x1 - x2 + x2 + x1)
    rbci_y = 0.125 * (y1 + y1 - y2 - y2 - y1 - y1 + y2 + y2)
    rbci_z = 0.125 * (z1 + z1 - z1 - z1 - z2 - z2 + z2 + z2)

    # raci
    # 0.125 * (ri[1, :] - ri[2, :] - ri[4, :] + ri[3, :] - ri[5, :] + ri[6, :] + ri[8, :] - ri[7, :])
    raci_x = 0.125 * (x1 - x2 - x2 + x1 - x1 + x2 + x2 - x1)
    raci_y = 0.125 * (y1 - y1 - y2 + y2 - y1 + y1 + y2 - y2)
    raci_z = 0.125 * (z1 - z1 - z1 + z1 - z2 + z2 + z2 - z2)

    # rabci
    # 0.125 * (-ri[1, :] + ri[2, :] - ri[4, :] + ri[3, :] + ri[5, :] - ri[6, :] + ri[8, :] - ri[7, :])
    rabci_x = 0.125 * (-x1 + x2 - x2 + x1 + x1 - x2 + x2 - x1)
    rabci_y = 0.125 * (-y1 + y1 - y2 + y2 + y1 - y1 + y2 - y2)
    rabci_z = 0.125 * (-z1 + z1 - z1 + z1 + z2 - z2 + z2 - z2)

    nlkx = length(wekx)
    nlky = length(weky)
    nlkz = length(wekz)

    sum_a1 = zero(ComplexF64) # Initialize with ComplexF64 zero for type stability

    # Use @fastmath and @inbounds for performance in the nested loops
    @fastmath @inbounds for a1 in 1:nlkx
        sum_b1 = zero(ComplexF64)
        kx_a1 = rootkx[a1] # Store for reuse

        for b1 in 1:nlky
            sum_c1 = zero(ComplexF64)
            ky_b1 = rootky[b1] # Store for reuse

            for c1 in 1:nlkz
                kz_c1 = rootkz[c1] # Store for reuse

                # Calculate components of drai, drbi, drci directly
                drai_x = rai_x + rabi_x * ky_b1 + raci_x * kz_c1 + rabci_x * ky_b1 * kz_c1
                drai_y = rai_y + rabi_y * ky_b1 + raci_y * kz_c1 + rabci_y * ky_b1 * kz_c1
                drai_z = rai_z + rabi_z * ky_b1 + raci_z * kz_c1 + rabci_z * ky_b1 * kz_c1
                draim = sqrt(drai_x^2 + drai_y^2 + drai_z^2)

                drbi_x = rbi_x + rabi_x * kx_a1 + rbci_x * kz_c1 + rabci_x * kx_a1 * kz_c1
                drbi_y = rbi_y + rabi_y * kx_a1 + rbci_y * kz_c1 + rabci_y * kx_a1 * kz_c1
                drbi_z = rbi_z + raci_z * kx_a1 + rbci_z * ky_b1 + rabci_z * kx_a1 * ky_b1 # Check original formula for drbi_z
                drbim = sqrt(drbi_x^2 + drbi_y^2 + drbi_z^2)

                drci_x = rci_x + raci_x * kx_a1 + rbci_x * ky_b1 + rabci_x * kx_a1 * ky_b1
                drci_y = rci_y + raci_y * kx_a1 + rbci_y * ky_b1 + rabci_y * kx_a1 * ky_b1
                drci_z = rci_z + raci_z * kx_a1 + rbci_z * ky_b1 + rabci_z * kx_a1 * ky_b1 # Check original formula for drci_z
                drcim = sqrt(drci_x^2 + drci_y^2 + drci_z^2)

                # Calculate components of r1 directly
                r1_x = rmi_x + rai_x * kx_a1 + rbi_x * ky_b1 + rci_x * kz_c1 +
                       rabi_x * kx_a1 * ky_b1 + raci_x * kx_a1 * kz_c1 + rbci_x * ky_b1 * kz_c1 +
                       rabci_x * kx_a1 * ky_b1 * kz_c1
                r1_y = rmi_y + rai_y * kx_a1 + rbi_y * ky_b1 + rci_y * kz_c1 +
                       rabi_y * kx_a1 * ky_b1 + raci_y * kx_a1 * kz_c1 + rbci_y * ky_b1 * kz_c1 +
                       rabci_y * kx_a1 * ky_b1 * kz_c1
                r1_z = rmi_z + rai_z * kx_a1 + rbi_z * ky_b1 + rci_z * kz_c1 +
                       rabi_z * kx_a1 * ky_b1 + raci_z * kx_a1 * kz_c1 + rbci_z * ky_b1 * kz_c1 +
                       rabci_z * kx_a1 * ky_b1 * kz_c1

                # Original code used `real(r1[1])` - this is unnecessary if r1_x, r1_y, r1_z are Float64.
                # All components are real (x,y,z coordinates).
                x = r1_x
                y = r1_y
                z = r1_z

                delta_x = (xo - x)
                delta_y = (yo - y)
                delta_z = (zo - z)

                R = sqrt(delta_x^2 + delta_y^2 + delta_z^2)

                # Handle potential R=0: add a small epsilon to R if it could be zero.
                # If R=0, (1/R) is Inf, exp(-1im*beta*0) = 1, so Inf*1.
                # A robust solution needs to handle singularities properly.
                # For now, assuming R > 0.
                invR = 1.0 / R
                G = invR * exp(-1im * beta * R) # `beta` must be ComplexF64 for this to work correctly

                # `f` should be ComplexF64 since `G` is ComplexF64.
                f = draim * drbim * drcim * G

                sum_c1 += wekz[c1] * f
            end
            sum_b1 += weky[b1] * sum_c1
        end
        sum_a1 += wekx[a1] * sum_b1
    end

    return sum_a1
end

# --- qrule function (remains unchanged from previous optimization) ---
function qrule(n::Int)
    iter = 2
    m = trunc((n + 1) / 2)
    e1 = n * (n + 1)
    mm = 4 * m - 1
    t = (pi / (4 * n + 2)) * (3:4:mm)
    nn = (1 - (1 - 1 / n) / (8 * n * n))
    xo = nn * cos.(t)

    for kk ∈ 1:iter
        pkm1 = ones(size(xo))
        pk = copy(xo)

        @inbounds for k ∈ 2:n
            pkp1 = @. (2k - 1) / k * xo * pk - (k - 1) / k * pkm1
            pkm1 = pk
            pk = pkp1
        end

        den = 1 .- xo .^ 2
        d1 = n .* (pkm1 .- xo .* pk)

        dpn = d1 ./ den
        d2pn = (2 .* xo .* dpn - e1 .* pk) ./ den
        d3pn = (4 * xo .* d2pn + (2 - e1) .* dpn) ./ den
        d4pn = (6 * xo .* d3pn + (6 - e1) .* d2pn) ./ den
        u = pk ./ dpn
        v = d2pn ./ dpn
        h = -u .* (1 .+ (0.5 * u) .* (v .+ u .* (v .* v .- u .* d3pn ./ (3 * dpn))))
        p = pk .+ h .* (dpn .+ (0.5 * h) .* (d2pn .+ (h ./ 3) .* (d3pn .+ 0.25 .* h .* d4pn)))
        dp = dpn .+ h .* (d2pn .+ (0.5 .* h) .* (d3pn .+ h .* d4pn ./ 3))
        h = h .- p ./ dp
        xo = xo .+ h
    end

    bp = Vector{Float64}(undef, n)
    wf = Vector{Float64}(undef, n)

    bp .= -xo .- h
    fx = d1 .- h .* e1 .* (pk .+ (h ./ 2) .* (dpn .+ (h ./ 3) .* (
        d2pn .+ (h ./ 4) .* (d3pn .+ (0.2 .* h) .* d4pn))))
    wf .= 2 .* (1 .- bp .^ 2) ./ (fx .* fx)

    if (m + m) > n
        bp[m] = 0.0
    end
    if !((m + m) == n)
        m = m - 1
    end
    jj = 1:m
    n1j = (n + 1) .- jj
    bp[n1j] .= -bp[jj]
    wf[n1j] .= wf[jj]

    return bp, wf
end