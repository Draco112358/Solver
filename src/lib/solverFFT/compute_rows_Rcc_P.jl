include("compute_row_P_sup.jl")
include("compute_P_vox_Rcc.jl")
using MKL
using FLoops

function compute_rows_Rcc_P(circulant_centers, id)
    sx = circulant_centers["sx"]
    sy = circulant_centers["sy"]
    sz = circulant_centers["sz"]

    Nx = circulant_centers["Nx"]
    Ny = circulant_centers["Ny"]
    Nz = circulant_centers["Nz"]

    # Create dictionaries to hold Nx, Ny, Nz values
    rows_P = Dict()
    rows_P["Nx"] = zeros(Int, 3, 3)
    rows_P["Ny"] = zeros(Int, 3, 3)
    rows_P["Nz"] = zeros(Int, 3, 3)

    rows_P["Nx"][1, 1] = Nx
    rows_P["Ny"][1, 1] = Ny + 1
    rows_P["Nz"][1, 1] = Nz

    rows_P["Nx"][2, 2] = Nx + 1
    rows_P["Ny"][2, 2] = Ny
    rows_P["Nz"][2, 2] = Nz

    rows_P["Nx"][3, 3] = Nx
    rows_P["Ny"][3, 3] = Ny
    rows_P["Nz"][3, 3] = Nz + 1

    rows_P["Nx"][1, 2] = 2 * (Nx + 1) - 1
    rows_P["Ny"][1, 2] = 2 * (Ny + 1) - 1
    rows_P["Nz"][1, 2] = Nz

    rows_P["Nx"][1, 3] = Nx
    rows_P["Ny"][1, 3] = 2 * (Ny + 1) - 1
    rows_P["Nz"][1, 3] = 2 * (Nz + 1) - 1

    rows_P["Nx"][2, 3] = 2 * (Nx + 1) - 1
    rows_P["Ny"][2, 3] = Ny
    rows_P["Nz"][2, 3] = 2 * (Nz + 1) - 1
    # Initialize QS and Rcc as dictionaries to store results
    rows_P["QS"] = Array{Vector{Float64}, 2}(undef, 3, 3)
    rows_P["Rcc"] = Array{Vector{Float64}, 2}(undef, 3, 3)
    

    size_subs = 10000

    # Compute sub-block indices and QS{1,1}
    compute_qs_entry!(rows_P, circulant_centers["p12_se"], sx, sy, sz, size_subs, 1, 1, id)
    compute_qs_entry!(rows_P, circulant_centers["p34_se"], sx, sy, sz, size_subs, 2, 2, id)
    compute_qs_entry!(rows_P, circulant_centers["p56_se"], sx, sy, sz, size_subs, 3, 3, id)
    compute_qs_entry!(rows_P, circulant_centers["p1234"], sx, sy, sz, size_subs, 1, 2, id)
    compute_qs_entry!(rows_P, circulant_centers["p1256"], sx, sy, sz, size_subs, 1, 3, id)
    compute_qs_entry!(rows_P, circulant_centers["p3456"], sx, sy, sz, size_subs, 2, 3, id)

    # Compute Rcc values
    rows_P["Rcc"][1, 1] = compute_P_vox_Rcc(circulant_centers["p12_se"][1, :], circulant_centers["p12_se"], sx, sy, sz, 1, 1, id)
    rows_P["Rcc"][2, 2] = compute_P_vox_Rcc(circulant_centers["p34_se"][1, :], circulant_centers["p34_se"], sx, sy, sz, 3, 3, id)
    rows_P["Rcc"][3, 3] = compute_P_vox_Rcc(circulant_centers["p56_se"][1, :], circulant_centers["p56_se"], sx, sy, sz, 5, 5, id)
    rows_P["Rcc"][1, 2] = compute_P_vox_Rcc(circulant_centers["p1234"][1, :], circulant_centers["p1234"], sx, sy, sz, 3, 2, id)
    rows_P["Rcc"][1, 3] = compute_P_vox_Rcc(circulant_centers["p1256"][1, :], circulant_centers["p1256"], sx, sy, sz, 5, 2, id)
    rows_P["Rcc"][2, 3] = compute_P_vox_Rcc(circulant_centers["p3456"][1, :], circulant_centers["p3456"], sx, sy, sz, 5, 4, id)

    return rows_P
end

# Helper function to compute entries for QS
function compute_qs_entry!(rows_P, p_se, sx, sy, sz, size_subs, row, col, id)
    num_ele_P_par = size(p_se, 1)
    VS, VE = generate_indices_sub_blocks(num_ele_P_par, size_subs)
    row_PS = Vector{Any}(undef, length(VS))
    coss = p_se[1, :]
    #deve essere fatto in parallelo
    for cpf in range(1, length(VS))
        row_PS[cpf] = compute_row_P_sup(coss, p_se[VS[cpf]:VE[cpf], :], sx, sy, sz, row, col, id)
    end
    rows_P["QS"][row, col] = zeros(num_ele_P_par)
    for cpf in range(1, length(VS))
        rows_P["QS"][row, col][VS[cpf]:VE[cpf]] = row_PS[cpf]
    end
end

# Function to generate indices sub-blocks
function generate_indices_sub_blocks(N, Ns)
    if Ns < N
        vect_ind_start = collect(1:Ns:N)
        vect_ind_end = vect_ind_end=collect(vect_ind_start[2]-1:Ns:N)
        if vect_ind_end[end] < N
            push!(vect_ind_end, N)
        end
    else
        vect_ind_start = [1]
        vect_ind_end = [N]
    end
    return vect_ind_start, vect_ind_end
end