function computeVs(
    time_json, # Renamed to avoid clash with Base.time
    time_delay_vs::Float64,
    signal_type_E,
    volumi::Dict{Symbol, AbstractArray},
    nodi_coord::Matrix{Float64},
    E::Vector{Float64},
    K::Vector{Float64},
    A::SparseMatrixCSC{Float64, Int64},
    tr::Int64,
    power::Int64,
    dev_stand::Int64,
    f0::Float64,
    ind_freq_interest::Vector{Int64}
)
    # --- OPTIMIZATION 6: Eagerly convert JSON types to standard Julia vectors
    time = Vector{Float64}(time_json)

    # Basic constants and dimensions
    c0 = 3e8
    num_cells = size(A, 1) # More robust way to get row count
    num_freqs = length(ind_freq_interest)
    
    # --- OPTIMIZATION 4: Pre-allocate output and temporary arrays
    Vs = zeros(ComplexF64, num_cells, num_freqs)
    ft = zeros(Float64, length(time))
    
    # Deconstruct vector E once
    E1, E2, E3 = E[1], E[2], E[3]

    # --- OPTIMIZATION 1 & 5: Setup FFT plan and pre-calculate frequency vector
    planfft = plan_fft(ft) # Plan is now based on the pre-allocated array
    freq_vector = fft_frequencies(time) # Calculate frequencies only once

    # --- OPTIMIZATION 2: Pre-calculate all geometric properties
    midpoints = zeros(Float64, num_cells, 3)
    cell_vectors = zeros(Float64, num_cells, 3)
    Base.Threads.@threads for k in 1:num_cells
        # --- OPTIMIZATION 3: Efficiently find non-zero elements in the sparse row
        row_indices, row_values = findnz(A[k, :])
        
        n1_idx_ptr = findfirst(v -> v == -1, row_values)
        n2_idx_ptr = findfirst(v -> v == 1, row_values)
        
        # Error handling for robustness
        if isnothing(n1_idx_ptr) || isnothing(n2_idx_ptr)
            error("Cell $k does not have valid node connections in sparse matrix A.")
        end

        n1_index = row_indices[n1_idx_ptr]
        n2_index = row_indices[n2_idx_ptr]

        n1 = nodi_coord[n1_index, :]
        n2 = nodi_coord[n2_index, :]

        midpoints[k, :] = (n1 .+ n2) ./ 2
        cell_vectors[k, :] = n2 .- n1
    end

    # Get signal type once, outside the loop
    signal_type = signal_type_E["type"]

    # --- MAIN LOOP (Now much lighter) ---
    Base.Threads.@threads for k in 1:num_cells
        r = @view midpoints[k, :]
        vett_cella = @view cell_vectors[k, :]
        
        # Time delay calculation for the current cell
        delay = time_delay_vs + dot(K, r) / c0
        time_shifted = time .- delay

        # In-place calculation of the time-domain signal
        if signal_type == "exponential"
            @. ft = (time_shifted / tr)^power * exp(-power * (time_shifted / tr - 1)) * (time >= delay)
        elseif signal_type == "gaussian_modulated"
            @. ft = cos(2 * pi * f0 * time_shifted) * exp(-(time_shifted)^2 / (2 * dev_stand^2)) * (time >= delay)
        elseif signal_type == "sinusoidal"
            @. ft = cos(2 * pi * f0 * time_shifted) * (time >= delay)
        end
        
        # --- OPTIMIZATION 1: Perform only ONE FFT per cell ---
        # The fft_UAq_opt function now only returns the complex spectrum
        Trasformata_ft_freq = fft_UAq(ft, planfft)

        # Reconstruct the full E-field spectrum using linearity
        # E0_freq(f) = [E1*FFT(ft); E2*FFT(ft); E3*FFT(ft)]
        # This is done implicitly in the final dot product calculation
        
        for f_idx in 1:num_freqs
            freq_index = ind_freq_interest[f_idx]
            
            # Get the single complex value from the FFT of ft at the desired frequency
            ft_val_freq = Trasformata_ft_freq[freq_index]
            
            # E0_freq at this frequency is [E1, E2, E3] * ft_val_freq
            # Vs = dot(E0_freq, vett_cella) = dot([E1, E2, E3] * ft_val_freq, vett_cella)
            # By scalar associativity: Vs = ft_val_freq * dot([E1, E2, E3], vett_cella)
            
            E_dot_vett = dot(E, vett_cella) # Dot product of two small vectors
            Vs[k, f_idx] = ft_val_freq * E_dot_vett
        end
    end

    return Vs
end