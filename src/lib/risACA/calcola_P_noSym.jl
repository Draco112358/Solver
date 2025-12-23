function calcola_P_noSym(m::Vector{Int}, n::Vector{Int}, superfici, escalings, id)
    # @save "inpuP.jld2" superfici escalings QS_Rcc_FW id
    eps0 = 8.854187816997944e-12

    epsilon1 = 5e-3
    epsilon2 = 1e-3
    epsilon3 = 1e-3
    epsilon4 = 3e-1

    use_suppression = 1
    # if isa(superfici["estremi_celle"], Vector)
    #     superfici["estremi_celle"] = hcat(superfici["estremi_celle"]...)'
    # end
    # if isa(superfici["centri"], Vector)
    #     superfici["centri"] = hcat(superfici["centri"]...)'
    # end
    nsup1 = length(m)
    nsup2 = length(n)
    # Inside calcola_P, after nsup is determined
    # Allocate arrays to store the calculated properties for each surface
    x_vs = Vector{Vector{Float64}}(undef, nsup1)
    y_vs = Vector{Vector{Float64}}(undef, nsup1)
    z_vs = Vector{Vector{Float64}}(undef, nsup1)

    xc_s = Vector{Float64}(undef, nsup1)
    yc_s = Vector{Float64}(undef, nsup1)
    zc_s = Vector{Float64}(undef, nsup1)

    a_s = Vector{Float64}(undef, nsup1)
    b_s = Vector{Float64}(undef, nsup1)
    c_s = Vector{Float64}(undef, nsup1)

    sup_s = Vector{Float64}(undef, nsup1)
    sup_xz_planes = Vector{Int}(undef, nsup1)
    sup_yz_planes = Vector{Int}(undef, nsup1)

    x_vs2 = Vector{Vector{Float64}}(undef, nsup2)
    y_vs2 = Vector{Vector{Float64}}(undef, nsup2)
    z_vs2 = Vector{Vector{Float64}}(undef, nsup2)

    xc_s2 = Vector{Float64}(undef, nsup2)
    yc_s2 = Vector{Float64}(undef, nsup2)
    zc_s2 = Vector{Float64}(undef, nsup2)

    a_s2 = Vector{Float64}(undef, nsup2)
    b_s2 = Vector{Float64}(undef, nsup2)
    c_s2 = Vector{Float64}(undef, nsup2)

    sup_s2 = Vector{Float64}(undef, nsup2)
    sup_xz_planes2 = Vector{Int}(undef, nsup2)
    sup_yz_planes2 = Vector{Int}(undef, nsup2)



    round_precision = 14 # Define this once

    process_surface_ACA!(x_vs, y_vs, z_vs, xc_s, yc_s, zc_s, a_s, b_s, c_s, sup_s, sup_xz_planes, sup_yz_planes, superfici, round_precision, m)
    process_surface_ACA!(x_vs2, y_vs2, z_vs2, xc_s2, yc_s2, zc_s2, a_s2, b_s2, c_s2, sup_s2, sup_xz_planes2, sup_yz_planes2, superfici, round_precision, n)
    block_size1 = 200  # ad esempio, 10 iterazioni per blocco
    P = zeros(nsup1, nsup2)

    calculate_P_matrix_noSym!(P, x_vs, y_vs, z_vs, xc_s, yc_s, zc_s, a_s, b_s, c_s, sup_s, sup_yz_planes, sup_xz_planes, Float64.(superfici["S"][m]), nsup1, x_vs2, y_vs2, z_vs2, xc_s2, yc_s2, zc_s2, a_s2, b_s2, c_s2, sup_s2, sup_yz_planes2, sup_xz_planes2, Float64.(superfici["S"][n]), nsup2, 0, block_size1, id, escalings[:P], epsilon1, epsilon2, epsilon3, epsilon4, use_suppression, eps0)
    return P
end

function calculate_P_matrix_noSym!(
    P_matrix::Union{Matrix{Float64},Matrix{ComplexF64}}, # P is typically real, if complex adjust type
    x_vs::Vector{Vector{Float64}}, y_vs::Vector{Vector{Float64}}, z_vs::Vector{Vector{Float64}},
    xc_s::Vector{Float64}, yc_s::Vector{Float64}, zc_s::Vector{Float64},
    a_s::Vector{Float64}, b_s::Vector{Float64}, c_s::Vector{Float64},
    sup_s::Vector{Float64}, sup_yz_planes::Vector{Int64}, sup_xz_planes::Vector{Int64},
    S_superfici::Vector{Float64}, # Pass subset superfici["S"][m]
    nsup::Int,
    x_vs2::Vector{Vector{Float64}}, y_vs2::Vector{Vector{Float64}}, z_vs2::Vector{Vector{Float64}},
    xc_s2::Vector{Float64}, yc_s2::Vector{Float64}, zc_s2::Vector{Float64},
    a_s2::Vector{Float64}, b_s2::Vector{Float64}, c_s2::Vector{Float64},
    sup_s2::Vector{Float64}, sup_yz_planes2::Vector{Int64}, sup_xz_planes2::Vector{Int64},
    S_superfici2::Vector{Float64}, # Pass subset superfici["S"][n]
    nsup2::Int,
    offset::Int, # Added for potential future use with offsets
    block_size2::Int,
    id::String,
    escalings_P,
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

            @inbounds for n_idx in 1:nsup2
                actual_n_idx = n_idx + offset

                integ, _ = Song_P_improved_Ivana_strategy(
                    x_vs[actual_m_idx], y_vs[actual_m_idx], z_vs[actual_m_idx],
                    xc_s[actual_m_idx], yc_s[actual_m_idx], zc_s[actual_m_idx],
                    a_s[actual_m_idx], b_s[actual_m_idx], c_s[actual_m_idx],
                    sup_s[actual_m_idx], sup_yz_planes[actual_m_idx], sup_xz_planes[actual_m_idx],
                    x_vs2[actual_n_idx], y_vs2[actual_n_idx], z_vs2[actual_n_idx],
                    xc_s2[actual_n_idx], yc_s2[actual_n_idx], zc_s2[actual_n_idx],
                    a_s2[actual_n_idx], b_s2[actual_n_idx], c_s2[actual_n_idx],
                    sup_s2[actual_n_idx], sup_yz_planes2[actual_n_idx], sup_xz_planes2[actual_n_idx],
                    epsilon1, epsilon2, epsilon3, epsilon4, use_suppression
                )

                # Use the directly passed S_superfici
                inv_factor = 1.0 / (4 * π * eps0_val * S_superfici[actual_m_idx] * S_superfici2[actual_n_idx])

                value = inv_factor * integ * escalings_P
                P_matrix[m_idx, n_idx] = value
            end
        end

        # Check for stop request after each block
        # if is_stop_requested(id)
        #     println("Simulazione $(id) interrotta per richiesta stop.")
        #     return nothing
        # end

        #println("Processed block $(block_idx) / $(num_blocks) for P calculation.")
    end

    return P_matrix
end