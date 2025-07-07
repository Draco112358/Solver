include("compute_Voxels_Rcc.jl")
include("compute_Lp_Voxels_QS.jl")
using MKL

function compute_rows_Rcc_Lp(circulant_centers, id)
    # Extract properties from circulant_centers
    sx = circulant_centers["sx"]
    sy = circulant_centers["sy"]
    sz = circulant_centers["sz"]

    rows_Lp = Dict()
    rows_Lp["Nx"] = circulant_centers["Nx"]
    rows_Lp["Ny"] = circulant_centers["Ny"]
    rows_Lp["Nz"] = circulant_centers["Nz"]

    Nlx = size(circulant_centers["Lpx"], 1)
    Nly = size(circulant_centers["Lpy"], 1)
    Nlz = size(circulant_centers["Lpz"], 1)

    size_subs = 10000

    dc = 1
    if Nlx > 0
        rows_Lp["Rcc_x"] = compute_Voxels_Rcc(circulant_centers["Lpx"][1, :], circulant_centers["Lpx"])

        num_ele_Lp_par = size(circulant_centers["Lpx"], 1)
        VS, VE = generate_indices_sub_blocks(num_ele_Lp_par, size_subs)

        Lpx_rowS = Vector{Vector{Float64}}(undef, length(VS))
        coss = circulant_centers["Lpx"][1, :]

        Threads.@threads for cpf in range(1,length(VS))
            Lpx_rowS[cpf] = compute_Lp_Voxels_QS(coss, circulant_centers["Lpx"][VS[cpf]:VE[cpf], :], sx, sy, sz, sx, sy, sz, dc,id)
            if isnothing(Lpx_rowS)
                return nothing
            end
        end

        rows_Lp["Lpx_QS"] = zeros(num_ele_Lp_par)
        for cpf in range(1,length(VS))
            rows_Lp["Lpx_QS"][VS[cpf]:VE[cpf]] .= Lpx_rowS[cpf]
        end

        rows_Lp["Lpx_QS"][1] = 0

    else
        rows_Lp["Rcc_x"] = zeros(0)
        rows_Lp["Lpx_QS"] = zeros(0)
        rows_Lp["Rcc_x_2"] = zeros(0)
        rows_Lp["Lpx_QS_2"] = zeros(0)
    end

    dc = 2
    if Nly > 0
        rows_Lp["Rcc_y"] = compute_Voxels_Rcc(circulant_centers["Lpy"][1, :], circulant_centers["Lpy"])

        num_ele_Lp_par = size(circulant_centers["Lpy"], 1)
        VS, VE = generate_indices_sub_blocks(num_ele_Lp_par, size_subs)

        Lpy_rowS = Vector{Vector{Float64}}(undef, length(VS))
        coss = circulant_centers["Lpy"][1, :]

        Threads.@threads for cpf in range(1,length(VS))
            Lpy_rowS[cpf] = compute_Lp_Voxels_QS(coss, circulant_centers["Lpy"][VS[cpf]:VE[cpf], :], sx, sy, sz, sx, sy, sz, dc, id)
            if isnothing(Lpy_rowS)
                return nothing
            end
        end

        rows_Lp["Lpy_QS"] = zeros(num_ele_Lp_par)
        for cpf in range(1, length(VS))
            rows_Lp["Lpy_QS"][VS[cpf]:VE[cpf]] .= Lpy_rowS[cpf]
        end

        rows_Lp["Lpy_QS"][1] = 0

    else
        rows_Lp["Rcc_y"] = zeros(0)
        rows_Lp["Lpy_QS"] = zeros(0)
        rows_Lp["Rcc_y_2"] = zeros(0)
        rows_Lp["Lpy_QS_2"] = zeros(0)
    end

    dc = 3
    if Nlz > 0
        rows_Lp["Rcc_z"] = compute_Voxels_Rcc(circulant_centers["Lpz"][1, :], circulant_centers["Lpz"])

        num_ele_Lp_par = size(circulant_centers["Lpz"], 1)
        VS, VE = generate_indices_sub_blocks(num_ele_Lp_par, size_subs)

        Lpz_rowS = Vector{Vector{Float64}}(undef, length(VS))
        coss = circulant_centers["Lpz"][1, :]

        Threads.@threads for cpf in range(1,length(VS))
            Lpz_rowS[cpf] = compute_Lp_Voxels_QS(coss, circulant_centers["Lpz"][VS[cpf]:VE[cpf], :], sx, sy, sz, sx, sy, sz, dc, id)
            if isnothing(Lpz_rowS)
                return nothing
            end
        end

        rows_Lp["Lpz_QS"] = zeros(num_ele_Lp_par)
        for cpf in range(1,length(VS))
            rows_Lp["Lpz_QS"][VS[cpf]:VE[cpf]] .= Lpz_rowS[cpf]
        end

        rows_Lp["Lpz_QS"][1] = 0

    else
        rows_Lp["Rcc_z"] = zeros(0)
        rows_Lp["Lpz_QS"] = zeros(0)
        rows_Lp["Rcc_z_2"] = zeros(0)
        rows_Lp["Lpz_QS_2"] = zeros(0)
    end

    return rows_Lp
end

function generate_indices_sub_blocks(N, Ns)
    if Ns < N
        vect_ind_start = 1:Ns:N
        vect_ind_end = vcat(vect_ind_start[2:end] .- 1, N)

        if vect_ind_end[end] < N
            vect_ind_end[end] = N
        end

    else
        vect_ind_start = [1]
        vect_ind_end = [N]
    end

    return vect_ind_start, vect_ind_end
end
