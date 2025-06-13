include("Song_P_improved_Ivana_strategy2.jl")

function calcola_P(superfici, escalings, QS_Rcc_FW, id)

	eps0 = 8.854187816997944e-12

	epsilon1 = 5e-3
	epsilon2 = 1e-3
	epsilon3 = 1e-3
	epsilon4 = 3e-1

	use_suppression = 1
	superfici["estremi_celle"] = hcat(superfici["estremi_celle"]...)
	superfici["centri"] = hcat(superfici["centri"]...)
	nsup = size(superfici["estremi_celle"], 1)
	# Inside calcola_P, after nsup is determined
	# Allocate arrays to store the calculated properties for each surface
	x_vs = Vector{Vector{Float64}}(undef, nsup)
	y_vs = Vector{Vector{Float64}}(undef, nsup)
	z_vs = Vector{Vector{Float64}}(undef, nsup)

	xc_s = Vector{Float64}(undef, nsup)
	yc_s = Vector{Float64}(undef, nsup)
	zc_s = Vector{Float64}(undef, nsup)

	a_s = Vector{Float64}(undef, nsup)
	b_s = Vector{Float64}(undef, nsup)
	c_s = Vector{Float64}(undef, nsup)

	sup_s = Vector{Float64}(undef, nsup)
	sup_xz_planes = Vector{Int}(undef, nsup)
	sup_yz_planes = Vector{Int}(undef, nsup)

	round_precision = 14 # Define this once

	process_surfaces!(x_vs, y_vs,z_vs,xc_s,yc_s,zc_s,a_s,b_s,c_s,sup_s,sup_xz_planes,sup_yz_planes,superfici,round_precision)

	R_cc = zeros(nsup, nsup)
	block_size1 = 200  # ad esempio, 10 iterazioni per blocco
	println("Calcolo P initialization")
	R_cc = calculate_R_cc!(R_cc, superfici["centri"], nsup, QS_Rcc_FW, id, block_size1)
	P = zeros(nsup, nsup)
	println("Calcolo P Song_P_improved_Ivana_strategy")

	calculate_P_matrix(P, x_vs, y_vs, z_vs, xc_s, yc_s, zc_s, a_s, b_s, c_s, sup_s, sup_yz_planes, sup_xz_planes, Float64.(superfici["S"]), nsup, 0, block_size1, id, escalings[:P], epsilon1, epsilon2, epsilon3, epsilon4, use_suppression, eps0)

	return Dict(
		:P => P,
		:R_cc => R_cc,
	)
end

function calculate_R_cc!(R_cc::Matrix{Float64}, centri::Matrix{Float64}, nsup::Int,
                         QS_Rcc_FW::Int, id::String, block_size1::Int)

    if QS_Rcc_FW < 2
        return nothing # Nothing to do if condition not met
    end

    num_blocks = ceil(Int, nsup / block_size1) # Calculate total number of blocks

    for (block_idx, m_block) in enumerate(1:block_size1:nsup)
        m_end = min(m_block + block_size1 - 1, nsup)

        Threads.@threads for m in m_block:m_end
            # Use StaticArrays for small, fixed-size vectors for performance.
            # No need for @view if you convert to SVector.
            centro_m = SVector{3}(centri[m, 1], centri[m, 2], centri[m, 3])
            # Or, if centri is guaranteed to be 3 columns, more robust:
            # centro_m = SVector{size(centri, 2)}(centri[m, :])

            @inbounds for n in m:nsup # @inbounds for potential speedup if bounds checks are not needed
                centro_n = SVector{3}(centri[n, 1], centri[n, 2], centri[n, 3])
                dist = norm(centro_m - centro_n)
                R_cc[m, n] = dist
                R_cc[n, m] = dist
            end
        end

        # Check for stop request after each block is fully processed
        if is_stop_requested(id)
            println("Simulazione $(id) interrotta per richiesta stop.")
            # Consider throwing an exception or returning a status for clearer interruption handling
            return nothing
        end

        # Use an actual counter for blocks for clearer output
        println("Processed block $(block_idx) / $(num_blocks)")
    end

    return R_cc # Return the populated matrix
end

function process_surfaces!(
    x_vs::Vector{Vector{Float64}}, y_vs::Vector{Vector{Float64}}, z_vs::Vector{Vector{Float64}},
    xc_s::Vector{Float64}, yc_s::Vector{Float64}, zc_s::Vector{Float64},
    a_s::Vector{Float64}, b_s::Vector{Float64}, c_s::Vector{Float64}, sup_s::Vector{Float64},
    sup_xz_planes::Vector{Int}, sup_yz_planes::Vector{Int},
    superfici::Dict{String, Any}, round_precision::Int
)
    nsup = size(superfici["estremi_celle"], 1)

    Threads.@threads for i in 1:nsup
        # Using a direct view is good, but for small, fixed-size access,
        # direct indexing might be slightly faster due to avoiding view object creation.
        # However, @view is generally robust. Let's stick with it for clarity.
        current_surface = @view superfici["estremi_celle"][i, :]

        # Extract values directly to avoid re-indexing `current_surface` multiple times
        x1 = current_surface[1]
        x10 = current_surface[10]
        y2 = current_surface[2]
        y11 = current_surface[11]
        z3 = current_surface[3]
        z12 = current_surface[12]

        # For fixed-size vectors (length 2), it's often more efficient to create them directly
        # than rely on `unique` which is more general.
        # Assuming `x_vs[i]` etc. are already `Vector{Float64}`.
        x_vs[i] = [round(min(x1, x10), digits=round_precision), round(max(x1, x10), digits=round_precision)]
        y_vs[i] = [round(min(y2, y11), digits=round_precision), round(max(y2, y11), digits=round_precision)]
        z_vs[i] = [round(min(z3, z12), digits=round_precision), round(max(z3, z12), digits=round_precision)]

        # Calculate centroids directly from the min/max values
        xc_s[i] = 0.5 * (x_vs[i][2] + x_vs[i][1]) # Assuming min/max order
        yc_s[i] = 0.5 * (y_vs[i][2] + y_vs[i][1])
        zc_s[i] = 0.5 * (z_vs[i][2] + z_vs[i][1])

        # Calculate absolute differences for side lengths
        a_val = abs(x_vs[i][2] - x_vs[i][1])
        b_val = abs(y_vs[i][2] - y_vs[i][1])
        c_val = abs(z_vs[i][2] - z_vs[i][1])

        sup1_yz_plane = 0
        sup1_xz_plane = 0

        # Use `isapprox` for floating-point comparisons if `a_val`, `b_val`, `c_val`
        # could be very close due to `round_precision`. For exact comparisons, `<= ` is fine.
        # Assuming exact comparisons are intended here.
        if (a_val <= b_val && a_val <= c_val)
            sup1_yz_plane = 1
            a_val = 1.0 # Use 1.0 for Float64 literal
        elseif (b_val <= a_val && b_val <= c_val)
            sup1_xz_plane = 1
            b_val = 1.0
        else
            c_val = 1.0
        end

        # Store results in the pre-allocated arrays
        a_s[i] = a_val
        b_s[i] = b_val
        c_s[i] = c_val
        sup_s[i] = a_val * b_val * c_val # This might not be a "surface" area but a volume if dimensions are set to 1.0
        sup_xz_planes[i] = sup1_xz_plane
        sup_yz_planes[i] = sup1_yz_plane
    end
    return nothing
end

function calculate_P_matrix(
    P_matrix::Matrix{Float64}, # P is typically real, if complex adjust type
    x_vs::Vector{Vector{Float64}}, y_vs::Vector{Vector{Float64}}, z_vs::Vector{Vector{Float64}},
    xc_s::Vector{Float64}, yc_s::Vector{Float64}, zc_s::Vector{Float64},
    a_s::Vector{Float64}, b_s::Vector{Float64}, c_s::Vector{Float64},
    sup_s::Vector{Float64}, sup_yz_planes::Vector{Int64}, sup_xz_planes::Vector{Int64},
    S_superfici::Vector{Float64}, # Pass superfici["S"] directly
    nsup::Int,
    offset::Int, # Added for potential future use with offsets
    block_size2::Int,
    id::String,
    escalings_P::Float64,
    epsilon1::Float64, epsilon2::Float64, epsilon3::Float64, epsilon4::Float64,
    use_suppression::Int64,
    eps0_val::Float64;
)
    # Validate dimensions
    #size(P_matrix) == (nsup, nsup) || throw(DimensionMismatch("P_matrix must be a $(nsup)x$(nsup) matrix."))
    # Add more dimension checks for input vectors if necessary (e.g., length(x_vs) >= nsup + offset)

    num_blocks = ceil(Int, nsup / block_size2)

    for (block_idx, m_block) in enumerate(1:block_size2:nsup)
        m_end = min(m_block + block_size2 - 1, nsup)

        Threads.@threads for m_idx in m_block:m_end
            # In this specific snippet, the offset is 0 as m goes from 1 to nsup,
            # so actual_m_idx == m_idx. But we keep it for consistency.
            actual_m_idx = m_idx + offset

            @inbounds for n_idx in m_idx:nsup
                actual_n_idx = n_idx + offset

                integ, _ = Song_P_improved_Ivana_strategy(
                    x_vs[actual_m_idx], y_vs[actual_m_idx], z_vs[actual_m_idx],
                    xc_s[actual_m_idx], yc_s[actual_m_idx], zc_s[actual_m_idx],
                    a_s[actual_m_idx], b_s[actual_m_idx], c_s[actual_m_idx],
                    sup_s[actual_m_idx], sup_yz_planes[actual_m_idx], sup_xz_planes[actual_m_idx],
                    x_vs[actual_n_idx], y_vs[actual_n_idx], z_vs[actual_n_idx],
                    xc_s[actual_n_idx], yc_s[actual_n_idx], zc_s[actual_n_idx],
                    a_s[actual_n_idx], b_s[actual_n_idx], c_s[actual_n_idx],
                    sup_s[actual_n_idx], sup_yz_planes[actual_n_idx], sup_xz_planes[actual_n_idx],
                    epsilon1, epsilon2, epsilon3, epsilon4, use_suppression
                )
                
                # Use the directly passed S_superfici
                inv_factor = 1.0 / (4 * π * eps0_val * S_superfici[actual_m_idx] * S_superfici[actual_n_idx])
                
                value = inv_factor * integ * escalings_P
                P_matrix[m_idx, n_idx] = value
                P_matrix[n_idx, m_idx] = value
            end
        end

        # Check for stop request after each block
        if is_stop_requested(id)
            println("Simulazione $(id) interrotta per richiesta stop.")
            return nothing
        end

        println("Processed block $(block_idx) / $(num_blocks) for P calculation.")
    end

    return P_matrix
end