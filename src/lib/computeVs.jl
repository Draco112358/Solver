using FFTW

function computeVs(time::Vector{Float64}, time_delay_vs::Float64, signal_type_E::Int64, volumi, nodi_coord::Matrix{Float64}, E::Vector{Float64}, K::Vector{Float64}, incidence_selection, tr::Float64, power::Float64, dev_stand::Float64, f0::Float64, ind_freq_interest::Vector{Int64})

    c0 = 3e8
    num_cells = size(volumi[:coordinate], 1)
    num_freqs = length(ind_freq_interest)
    Vs = zeros(num_cells, num_freqs)

    for k in 1:num_cells
        n1_index = findfirst(x -> x == -1, incidence_selection[:A][k,:])
        n2_index = findfirst(x -> x == 1, incidence_selection[:A][k,:])

        n1 = nodi_coord[n1_index, :]
        n2 = nodi_coord[n2_index, :]

        r = reshape((n2 + n1) / 2, 1, 3) # midpoint of the cell
        vett_cella = reshape(n2 - n1, 1, 3) # vector along the cell direction

        ft = zeros(length(time))
        time_shifted = time .- time_delay_vs .- dot(K, r) / c0

        if signal_type_E == 1
            # time evolution of the electric field for exponential signal
            ft .= (time_shifted ./ tr).^power .* exp.(-power .* (time_shifted ./ tr .- 1)) .* (time_shifted .>= 0)
        elseif signal_type_E == 2
            # time evolution of the electric field for modulated Gaussian signal
            ft .= cos.(2 * pi * f0 * time_shifted) .* exp.(-(time_shifted).^2 / (2 * dev_stand^2)) .* (time_shifted .>= 0)
        elseif signal_type_E == 3
            # time evolution of the electric field for sinusoidal signal
            ft .= cos.(2 * pi * f0 * time_shifted) .* (time_shifted .>= 0)
        end

        Ei = [E[1] .* ft; E[2] .* ft; E[3] .* ft]

        Eix = reshape(Ei[1,:], 1, length(time))
        Eiy = reshape(Ei[2,:], 1, length(time))
        Eiz = reshape(Ei[3,:], 1, length(time))

        Trasformata_Ex = fft_UAq(time, Eix)
        Trasformata_Ey = fft_UAq(time, Eiy)
        Trasformata_Ez = fft_UAq(time, Eiz)

        Eix_freq = Trasformata_Ex[2, :]
        Eiy_freq = Trasformata_Ey[2, :]
        Eiz_freq = Trasformata_Ez[2, :]

        E0_freq = [Eix_freq; Eiy_freq; Eiz_freq]

        for f in 1:num_freqs
            Vs[k, f] = dot(E0_freq[:, ind_freq_interest[f]], vec(vett_cella)) # Use vec to ensure it's a vector
        end
    end

    return Vs
end

function fft_UAq(t, x) # Assuming x can be complex
    fintem = t[end] - t[1]          # Length of the time window
    ncampt = length(t)              # Number of samples in time
    dtem = t[2] - t[1]              # Interval between two consecutive samples
  
    frecamp = 1 / dtem              # Sampling frequency
    fremax = frecamp / 2            # Nyquist frequency
    frefond = 1 / fintem            # Frequency resolution
    # Perform the FFT
    XAA = fft(x) * dtem
    XA = XAA[1:floor(Int(ncampt/2)) + 1]
  
    # The MATLAB code comments out doubling the spectrum.
    # If you need to do this, uncomment the following line:
    # XA[2:end] = 2 * XA[2:end]
  
    f = (0:floor(Int(ncampt/2))) * frefond
    Trasformata = [transpose(f); transpose(XA)] # Transpose to match MATLAB's column-wise output
    return Trasformata
  end