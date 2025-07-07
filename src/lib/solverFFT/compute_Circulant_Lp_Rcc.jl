include("build_row_vox_Rcc.jl")
using MKL
using LinearAlgebra, FFTW

# Main function to compute FFT of Circulant matrices for Lp using Rcc
function compute_Circulant_Lp_Rcc(rows_Lp, escalings, freq)
    escaling = escalings["Lp"]

    Nx = rows_Lp["Nx"]
    Ny = rows_Lp["Ny"]
    Nz = rows_Lp["Nz"]

    Nz2 = Nz

    FFTCLp = Array{Array{ComplexF64}}(undef, 3, 2)

    Lpx_row = escaling * build_row_vox_Rcc(rows_Lp["Lpx_QS"], rows_Lp["Rcc_x"], freq)
    CLpx = store_circulant(Lpx_row, Nx - 1, Ny, Nz)
    FFTCLp[1, 1] = fft(CLpx)

    Lpy_row = escaling * build_row_vox_Rcc(rows_Lp["Lpy_QS"], rows_Lp["Rcc_y"], freq)
    CLpy = store_circulant(Lpy_row, Nx, Ny - 1, Nz)
    FFTCLp[2, 1] = fft(CLpy)

    Lpz_row = escaling * build_row_vox_Rcc(rows_Lp["Lpz_QS"], rows_Lp["Rcc_z"], freq)
    CLpz = store_circulant(Lpz_row, Nx, Ny, Nz2 - 1)
    FFTCLp[3, 1] = fft(CLpz)

    return FFTCLp
end

# Function to store circulant matrix
function store_circulant(row_Lp, Nx, Ny, Nz)
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
    for cont in range(1,length(row_Lp))
        (m, n, k) = From_1D_to_3D(Nx, Ny, cont)
        Circ[m, n, k] = row_Lp[cont]
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