# include("fft_UAq.jl")
# using FFTW, LinearAlgebra, JSON3

# function computeVs(time::JSON3.Array{Union{Float64, Int64}, Vector{UInt8}, SubArray{UInt64, 1, Vector{UInt64}, Tuple{UnitRange{Int64}}, true}},
#     time_delay_vs::Float64, signal_type_E::JSON3.Object{Vector{UInt8}, SubArray{UInt64, 1, Vector{UInt64}, Tuple{UnitRange{Int64}}, true}},
#     volumi::Dict{Symbol, AbstractArray}, nodi_coord::Matrix{Float64}, E::Vector{Float64}, K::Vector{Float64},
#     A::SparseMatrixCSC{Float64, Int64}, tr::Int64, power::Int64, dev_stand::Int64, f0::Float64, ind_freq_interest::Vector{Int64})

#     c0 = 3e8
#     num_cells = size(volumi[:coordinate], 1)
#     num_freqs = length(ind_freq_interest)
#     Vs = zeros(ComplexF64, num_cells, num_freqs)

#     E1 = E[1]
#     E2 = E[2]
#     E3 = E[3]

#     planfft = plan_fft(E1 .* zeros(length(time)))

#     for k in 1:num_cells
#         n1_index = findfirst(x -> x == -1, A[k,:])
#         n2_index = findfirst(x -> x == 1, A[k,:])

#         n1 = nodi_coord[n1_index, :]
#         n2 = nodi_coord[n2_index, :]

#         r = reshape((n2 + n1) / 2, 1, 3) # midpoint of the cell
#         #dump(r)
#         vett_cella = reshape(n2 - n1, 1, 3) # vector along the cell direction

#         ft = zeros(length(time))
#         time_shifted = time .- time_delay_vs .- dot(K, r) / c0
        
#         if signal_type_E["type"] == "exponential"
#             # time evolution of the electric field for exponential signal
#             ft .= (time_shifted ./ tr).^power .* exp.(-power .* (time_shifted ./ tr .- 1)) .* (time .>= (time_delay_vs .+ dot(K, r) / c0))
#         elseif signal_type_E["type"] == "gaussian_modulated"
#             # time evolution of the electric field for modulated Gaussian signal
#             ft .= cos.(2 * pi * f0 * time_shifted) .* exp.(-(time_shifted).^2 / (2 * dev_stand^2)) .* (time .>= (time_delay_vs .+ dot(K, r) / c0))
#         elseif signal_type_E["type"] == "sinusoidal"
#             # time evolution of the electric field for sinusoidal signal
#             ft .= cos.(2 * pi * f0 * time_shifted) .* (time .>= (time_delay_vs .+ dot(K, r) / c0))
#         end

#         Eix = E1 .* ft
#         Eiy = E2 .* ft
#         Eiz = E3 .* ft

#         Trasformata_Ex = fft_UAq(time, Eix, planfft)
#         Trasformata_Ey = fft_UAq(time, Eiy, planfft)
#         Trasformata_Ez = fft_UAq(time, Eiz, planfft)

#         Eix_freq = Trasformata_Ex[2, :]
#         Eiy_freq = Trasformata_Ey[2, :]
#         Eiz_freq = Trasformata_Ez[2, :]

#         E0_freq = transpose([Eix_freq Eiy_freq Eiz_freq])


#         for f in 1:num_freqs
#             Vs[k, f] = dot(E0_freq[:, ind_freq_interest[f]], vec(vett_cella)) # Use vec to ensure it's a vector
#         end
#     end

#     return Vs
# end

include("fft_UAq.jl") # Using the optimized version of the helper
using FFTW, LinearAlgebra, JSON3, SparseArrays

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