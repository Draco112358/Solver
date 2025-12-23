function calcola_Lp_Sym(m::Vector{Int}, volumi, incidence_selection, escalings, id, component::String="all")
    epsilon1 = 5e-3
    epsilon2 = 1e-3
    epsilon3 = 1e-3
    epsilon4 = 3e-1

    use_suppression = true

    mx = incidence_selection[:mx]
    my = incidence_selection[:my]
    mz = incidence_selection[:mz]
    # ntot = mx + my + mz # Numero totale di volumi

    ntot = length(m)

    # --- INIZIO PRE-CALCOLO DELLE PROPRIETÀ DEI VOLUMI ---
    # Questa parte viene eseguita una sola volta per tutti i volumi.
    round_precision = 14

    global_offset = 0
    if component == "y"
        global_offset = mx
    elseif component == "z"
        global_offset = mx + my
    end

    # Definisci gli indici una sola volta
    indici_x_min = [1, 7, 13, 19]
    indici_y_min = [2, 5, 14, 17]
    indici_z_min = [3, 6, 9, 12]
    indici_x_max = [4, 10, 16, 22]
    indici_y_max = [8, 11, 20, 23]
    indici_z_max = [15, 18, 21, 24]

    x_vs = Vector{Vector{Float64}}(undef, ntot)
    y_vs = Vector{Vector{Float64}}(undef, ntot)
    z_vs = Vector{Vector{Float64}}(undef, ntot)

    xc_s = Vector{Float64}(undef, ntot)
    yc_s = Vector{Float64}(undef, ntot)
    zc_s = Vector{Float64}(undef, ntot)

    a_s = Vector{Float64}(undef, ntot)
    b_s = Vector{Float64}(undef, ntot)
    c_s = Vector{Float64}(undef, ntot)
    V_s = Vector{Float64}(undef, ntot)

    # --- INIZIO PRE-CALCOLO DEI PARAMETRI PER Compute_Lp_Self ---
    l_s = Vector{Float64}(undef, ntot)
    W_s = Vector{Float64}(undef, ntot)
    T_s = Vector{Float64}(undef, ntot)
    # --- FINE INIZIO PRE-CALCOLO ---

    Base.Threads.@threads for i in 1:ntot
        idx_global = m[i] + global_offset
        current_volume_coords = @view volumi[:coordinate][idx_global, :] # Usa @view per evitare copie

        x_vs[i] = unique(round.(current_volume_coords[[indici_x_min; indici_x_max]], digits=round_precision))
        y_vs[i] = unique(round.(current_volume_coords[[indici_y_min; indici_y_max]], digits=round_precision))
        z_vs[i] = unique(round.(current_volume_coords[[indici_z_min; indici_z_max]], digits=round_precision))

        xc_s[i] = 0.5 * (x_vs[i][end] + x_vs[i][1])
        yc_s[i] = 0.5 * (y_vs[i][end] + y_vs[i][1])
        zc_s[i] = 0.5 * (z_vs[i][end] + z_vs[i][1])

        a_s[i] = abs(x_vs[i][end] - x_vs[i][1])
        b_s[i] = abs(y_vs[i][end] - y_vs[i][1])
        c_s[i] = abs(z_vs[i][end] - z_vs[i][1])

        V_s[i] = a_s[i] * b_s[i] * c_s[i]
        # Pre-calcolo per Compute_Lp_Self
        curr_dir_temp = 0 # Questo deve essere gestito correttamente in base al contesto del volume

        if component == "all"
            if m[i] <= mx
                curr_dir_temp = 1
            elseif m[i] <= mx + my
                curr_dir_temp = 2
            else
                curr_dir_temp = 3
            end
        elseif component == "x"
            curr_dir_temp = 1
        elseif component == "y"
            curr_dir_temp = 2
        elseif component == "z"
            curr_dir_temp = 3
        end

        if curr_dir_temp == 1
            l_s[i] = abs(mean_length_rev(current_volume_coords, 1))
            W_s[i] = abs(mean_length_rev(current_volume_coords, 3))
            T_s[i] = abs(mean_length_rev(current_volume_coords, 2))
        elseif curr_dir_temp == 2
            T_s[i] = abs(mean_length_rev(current_volume_coords, 1))
            W_s[i] = abs(mean_length_rev(current_volume_coords, 3))
            l_s[i] = abs(mean_length_rev(current_volume_coords, 2))
        else # curr_dir_temp == 3
            W_s[i] = abs(mean_length_rev(current_volume_coords, 1))
            l_s[i] = abs(mean_length_rev(current_volume_coords, 3))
            T_s[i] = abs(mean_length_rev(current_volume_coords, 2))
        end
    end
    # --- FINE PRE-CALCOLO ---

    block_size = 200

    if component == "x"
        Lp_x = Matrix{Float64}(undef, ntot, ntot)
        calculate_Lp_matrix_Sym(Lp_x, l_s, W_s, T_s, x_vs, y_vs, z_vs, xc_s, yc_s, zc_s, a_s, b_s, c_s, V_s, Float64.(volumi[:S][m.+global_offset]), ntot, 0, block_size, id, escalings[:Lp], epsilon1, epsilon2, epsilon3, epsilon4, use_suppression, component, mx, my)
        return Lp_x
    elseif component == "y"
        Lp_y = Matrix{Float64}(undef, ntot, ntot)
        calculate_Lp_matrix_Sym(Lp_y, l_s, W_s, T_s, x_vs, y_vs, z_vs, xc_s, yc_s, zc_s, a_s, b_s, c_s, V_s, Float64.(volumi[:S][m.+global_offset]), ntot, 0, block_size, id, escalings[:Lp], epsilon1, epsilon2, epsilon3, epsilon4, use_suppression, component, mx, my)
        return Lp_y
    elseif component == "z"
        Lp_z = Matrix{Float64}(undef, ntot, ntot)
        calculate_Lp_matrix_Sym(Lp_z, l_s, W_s, T_s, x_vs, y_vs, z_vs, xc_s, yc_s, zc_s, a_s, b_s, c_s, V_s, Float64.(volumi[:S][m.+global_offset]), ntot, 0, block_size, id, escalings[:Lp], epsilon1, epsilon2, epsilon3, epsilon4, use_suppression, component, mx, my)
        return Lp_z
    else # component == "all"
        Lp_x = Matrix{Float64}(undef, mx, mx)
        Lp_y = Matrix{Float64}(undef, my, my)
        Lp_z = Matrix{Float64}(undef, mz, mz)

        calculate_Lp_matrix_Sym(Lp_x, l_s, W_s, T_s, x_vs, y_vs, z_vs, xc_s, yc_s, zc_s, a_s, b_s, c_s, V_s, Float64.(volumi[:S][m]), mx, 0, block_size, id, escalings[:Lp], epsilon1, epsilon2, epsilon3, epsilon4, use_suppression, "x", mx, my)
        calculate_Lp_matrix_Sym(Lp_y, l_s, W_s, T_s, x_vs, y_vs, z_vs, xc_s, yc_s, zc_s, a_s, b_s, c_s, V_s, Float64.(volumi[:S][m]), my, mx, block_size, id, escalings[:Lp], epsilon1, epsilon2, epsilon3, epsilon4, use_suppression, "y", mx, my)
        calculate_Lp_matrix_Sym(Lp_z, l_s, W_s, T_s, x_vs, y_vs, z_vs, xc_s, yc_s, zc_s, a_s, b_s, c_s, V_s, Float64.(volumi[:S][m]), mz, mx + my, block_size, id, escalings[:Lp], epsilon1, epsilon2, epsilon3, epsilon4, use_suppression, "z", mx, my)

        return Dict(
            :Lp_x => Lp_x,
            :Lp_y => Lp_y,
            :Lp_z => Lp_z
        )
    end
end

function calculate_Lp_matrix_Sym(
    Lp_matrix::Union{Matrix{Float64},Matrix{ComplexF64}}, # Assuming Lp can store complex values later
    l_s::Vector{Float64}, W_s::Vector{Float64}, T_s::Vector{Float64},
    x_vs::Vector{Vector{Float64}}, y_vs::Vector{Vector{Float64}}, z_vs::Vector{Vector{Float64}},
    xc_s::Vector{Float64}, yc_s::Vector{Float64}, zc_s::Vector{Float64},
    a_s::Vector{Float64}, b_s::Vector{Float64}, c_s::Vector{Float64},
    V_s::Vector{Float64}, S_volumi::Vector{Float64},
    dim::Int,
    offset::Int,
    block_size::Int,
    id::String,
    escalings_Lp::Int64,
    epsilon1::Float64, epsilon2::Float64, epsilon3::Float64, epsilon4::Float64,
    use_suppression::Bool,
    component::String,
    mx::Int, my::Int
)
    global_offset = 0
    if component == "y"
        global_offset = mx
    elseif component == "z"
        global_offset = mx + my
    end

    # Validate dimensions
    #size(Lp_matrix) == (dim, dim) || throw(DimensionMismatch("Lp_matrix must be a $(dim)x$(dim) matrix."))
    # Add more dimension checks for input vectors if necessary (e.g., length(l_s) >= dim + offset)
    num_blocks = ceil(Int, dim / block_size)

    for (block_idx, m_block) in enumerate(1:block_size:dim)
        m_end = min(m_block + block_size - 1, dim)

        Threads.@threads for m_idx in m_block:m_end
            # Calculate actual indices in the input arrays
            actual_m_idx = m_idx + offset

            # Self-inductance (diagonal element)
            Lp_matrix[m_idx, m_idx] = Compute_Lp_Self(l_s[actual_m_idx], W_s[actual_m_idx], T_s[actual_m_idx]) * escalings_Lp
            # Mutual inductance (off-diagonal elements)
            @inbounds for n_idx in m_idx+1:dim
                actual_n_idx = n_idx + offset
                integ, _ = Song_improved_Ivana_strategy2(
                    x_vs[actual_m_idx], y_vs[actual_m_idx], z_vs[actual_m_idx],
                    xc_s[actual_m_idx], yc_s[actual_m_idx], zc_s[actual_m_idx],
                    a_s[actual_m_idx], b_s[actual_m_idx], c_s[actual_m_idx], V_s[actual_m_idx],
                    x_vs[actual_n_idx], y_vs[actual_n_idx], z_vs[actual_n_idx],
                    xc_s[actual_n_idx], yc_s[actual_n_idx], zc_s[actual_n_idx],
                    a_s[actual_n_idx], b_s[actual_n_idx], c_s[actual_n_idx], V_s[actual_n_idx],
                    epsilon1, epsilon2, epsilon3, epsilon4, use_suppression
                )
                #Lp deve essere positiva, se integ risultasse negativo mìbasta mettere abs nella formula sotto
                value = 1e-7 / (S_volumi[actual_m_idx] * S_volumi[actual_n_idx]) * integ * escalings_Lp
                Lp_matrix[m_idx, n_idx] = value
                Lp_matrix[n_idx, m_idx] = value
            end
        end

        # Check for stop request after each block is fully processed
        # if is_stop_requested(id)
        #     println("Simulazione $(id) interrotta per richiesta stop.")
        #     return nothing
        # end

        #println("Processed block $(block_idx) / $(num_blocks) for Lp$(component) calculation.")
    end

    return Lp_matrix
end