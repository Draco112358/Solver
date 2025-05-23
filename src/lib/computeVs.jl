include("fft_UAq.jl")
using FFTW, LinearAlgebra

function computeVs(time, time_delay_vs, signal_type_E, volumi, nodi_coord, E, K, A, tr, power, dev_stand, f0, ind_freq_interest)

    c0 = 3e8
    num_cells = size(volumi[:coordinate], 1)
    num_freqs = length(ind_freq_interest)
    Vs = zeros(ComplexF64, num_cells, num_freqs)

    for k in 1:num_cells
        n1_index = findfirst(x -> x == -1, A[k,:])
        n2_index = findfirst(x -> x == 1, A[k,:])

        n1 = nodi_coord[n1_index, :]
        n2 = nodi_coord[n2_index, :]

        r = reshape((n2 + n1) / 2, 1, 3) # midpoint of the cell
        #dump(r)
        vett_cella = reshape(n2 - n1, 1, 3) # vector along the cell direction

        ft = zeros(length(time))
        time_shifted = time .- time_delay_vs .- dot(K, r) / c0
        
        if signal_type_E["type"] == "exponential"
            # time evolution of the electric field for exponential signal
            ft .= (time_shifted ./ tr).^power .* exp.(-power .* (time_shifted ./ tr .- 1)) .* (time .>= (time_delay_vs .+ dot(K, r) / c0))
        elseif signal_type_E["type"] == "gaussian_modulated"
            # time evolution of the electric field for modulated Gaussian signal
            ft .= cos.(2 * pi * f0 * time_shifted) .* exp.(-(time_shifted).^2 / (2 * dev_stand^2)) .* (time .>= (time_delay_vs .+ dot(K, r) / c0))
        elseif signal_type_E["type"] == "sinusoidal"
            # time evolution of the electric field for sinusoidal signal
            ft .= cos.(2 * pi * f0 * time_shifted) .* (time .>= (time_delay_vs .+ dot(K, r) / c0))
        end

        # Ei = [E[1] .* ft; E[2] .* ft; E[3] .* ft]

        # Eix = reshape(Ei[1,:], 1, length(time))
        # Eiy = reshape(Ei[2,:], 1, length(time))
        # Eiz = reshape(Ei[3,:], 1, length(time))

        Eix = E[1] .* ft
        Eiy = E[2] .* ft
        Eiz = E[3] .* ft

        Trasformata_Ex = fft_UAq(time, Eix)
        Trasformata_Ey = fft_UAq(time, Eiy)
        Trasformata_Ez = fft_UAq(time, Eiz)

        Eix_freq = Trasformata_Ex[2, :]
        Eiy_freq = Trasformata_Ey[2, :]
        Eiz_freq = Trasformata_Ez[2, :]

        E0_freq = transpose([Eix_freq Eiy_freq Eiz_freq])


        for f in 1:num_freqs
            Vs[k, f] = dot(E0_freq[:, ind_freq_interest[f]], vec(vett_cella)) # Use vec to ensure it's a vector
        end
    end

    return Vs
end