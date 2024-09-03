using LinearAlgebra
using FFTW

function compute_Circulant_P_sup_Rcc(rows_P, escalings, freq)
    # Extract scaling factor
    escaling = escalings["P"]
    
    # Extract dimensions
    Nx = rows_P["Nx"]
    Ny = rows_P["Ny"]
    Nz = rows_P["Nz"]
    
    # Initialize FFTCP as a 3x3 array of empty arrays
    FFTCP = Array{Array{ComplexF64}}(undef, 3, 3)
    
    # Loop over combinations of c1 and c2
    for c1 in 1:3
        for c2 in c1:3
            # Build row of voxels for the given indices and frequency
            row_P = escaling * build_row_vox_Rcc(rows_P["QS"][c1, c2], rows_P["Rcc"][c1, c2], freq)
            
            # Create the circulant matrix
            CP = store_circulant(row_P, Nx[c1, c2], Ny[c1, c2], Nz[c1, c2])
            
            # Compute the FFT of the circulant matrix
            FFTCP[c1, c2] = fft(CP)
        end
    end
    
    return FFTCP
end

function store_circulant(row_P, Nx, Ny, Nz)
    i1x = 1:Nx
    i2x = (Nx + 2):(2 * Nx)
    i3x = Nx:-1:2

    i1y = 1:Ny
    i2y = (Ny + 2):(2 * Ny)
    i3y = Ny:-1:2

    i1z = 1:Nz
    i2z = (Nz + 2):(2 * Nz)
    i3z = Nz:-1:2

    Circ = zeros(Nx, Ny, Nz)
    for cont in range(1,length(row_P))
        (m, n, k) = From_1D_to_3D(Nx, Ny, cont)
        Circ[m, n, k] = row_P[cont]
    end

    Cout = zeros(2 * Nx, 2 * Ny, 2 * Nz)
    Cout[i1x, i1y, i1z] = Circ[i1x, i1y, i1z]
    Cout[i2x, i1y, i1z] = Circ[i3x, i1y, i1z]
    Cout[i1x, i2y, i1z] = Circ[i1x, i3y, i1z]
    Cout[i1x, i1y, i2z] = Circ[i1x, i1y, i3z]
    Cout[i2x, i2y, i1z] = Circ[i3x, i3y, i1z]
    Cout[i2x, i1y, i2z] = Circ[i3x, i1y, i3z]
    Cout[i1x, i2y, i2z] = Circ[i1x, i3y, i3z]
    Cout[i2x, i2y, i2z] = Circ[i3x, i3y, i3z]

    return Cout
end

function from_1D_to_3D(M, N, pos)
    pos = pos - 1
    k = floor(pos / (M * N))
    j = floor((pos - k * M * N) / M)
    i = mod(pos - k * M * N, M)
    k = k + 1
    j = j + 1
    i = i + 1
    return convert(Int64, i), convert(Int64, j), convert(Int64, k)
end