function calcola_Lp_ACA_cpp(volumi, incidence_selection, ACA_thres, escalings, id)

    use_suppression = 1

    # Ensure Matrix format for volumi fields
    for field in [:centri, :coordinate]
        if haskey(volumi, field)
            data_field = volumi[field]
            if isa(data_field, Matrix) && eltype(data_field) <: Vector
                volumi[field] = reduce(vcat, transpose.(data_field[:, 1]))
            elseif isa(data_field, Vector) && eltype(data_field) <: Vector
                volumi[field] = reduce(vcat, transpose.(data_field))
            end
        end
    end

    mx = incidence_selection[:mx]
    my = incidence_selection[:my]
    mz = incidence_selection[:mz]

    Lp_data = Dict{Symbol,Any}(
        :mx => mx,
        :my => my,
        :mz => mz
    )

    # --- Lpx calculation ---
    Level_Lpx = ceil(Int, log2(mx / 600))
    if Level_Lpx < 1
        Level_Lpx = 1
    end
    basis_func_boxes_Lpx, Level_Lpx = prepare_multilevel(volumi[:centri][1:mx, :], Level_Lpx)
    Lp_data[:Lpx_comp] = last_level_fill(basis_func_boxes_Lpx, Level_Lpx)

    # --- Lpy calculation ---
    Level_Lpy = ceil(Int, log2(my / 600))
    if Level_Lpy < 1
        Level_Lpy = 1
    end
    basis_func_boxes_Lpy, Level_Lpy = prepare_multilevel(volumi[:centri][mx+1:mx+my, :], Level_Lpy)
    Lp_data[:Lpy_comp] = last_level_fill(basis_func_boxes_Lpy, Level_Lpy)

    # --- Lpz calculation ---
    Level_Lpz = ceil(Int, log2(mz / 600))
    if Level_Lpz < 1
        Level_Lpz = 1
    end
    basis_func_boxes_Lpz, Level_Lpz = prepare_multilevel(volumi[:centri][mx+my+1:mx+my+mz, :], Level_Lpz)
    Lp_data[:Lpz_comp] = last_level_fill(basis_func_boxes_Lpz, Level_Lpz)


    # # --- Lpx Matrix elements ---
    # estremi_celle = volumi[:coordinate][1:mx, :]
    # versori = hcat(ones(mx, 1), zeros(mx, 2))
    # centri = volumi[:centri][1:mx, :]

    num_blocks = length(Lp_data[:Lpx_comp])
    Lp_data[:Lpx] = Dict{Symbol,Any}(
        :D => Vector{Matrix{Float64}}(undef, num_blocks),
        :U => Matrix{Union{Nothing,Matrix{Float64}}}(nothing, num_blocks, num_blocks),
        :V => Matrix{Union{Nothing,Matrix{Float64}}}(nothing, num_blocks, num_blocks)
    )

    for c1 in 1:num_blocks
        m = Lp_data[:Lpx_comp][c1]

        if isempty(m)
            Lp_data[:Lpx][:D][c1] = Array{Float64}(undef, 0, 0)
        else
            Lp_data[:Lpx][:D][c1] = calcola_Lp_Sym(m, volumi, incidence_selection, escalings, id, "x")
        end
    end

    for c1 in 1:num_blocks-1
        for c2 in c1+1:num_blocks
            if !isempty(Lp_data[:Lpx_comp][c1]) && !isempty(Lp_data[:Lpx_comp][c2])
                Lp_data[:Lpx][:U][c1, c2], Lp_data[:Lpx][:V][c1, c2] = ACA_Lp(ACA_thres, Lp_data[:Lpx_comp][c1], Lp_data[:Lpx_comp][c2], volumi, incidence_selection, escalings, id, use_suppression, "x")
            end
        end
    end

    # --- Lpy Matrix elements ---
    # estremi_celle = volumi[:coordinate][mx+1:mx+my, :]
    # versori = hcat(zeros(my, 1), ones(my, 1), zeros(my, 1))
    # centri = volumi[:centri][mx+1:mx+my, :]

    num_blocks = length(Lp_data[:Lpy_comp])
    Lp_data[:Lpy] = Dict{Symbol,Any}(
        :D => Vector{Matrix{Float64}}(undef, num_blocks),
        :U => Matrix{Union{Nothing,Matrix{Float64}}}(nothing, num_blocks, num_blocks),
        :V => Matrix{Union{Nothing,Matrix{Float64}}}(nothing, num_blocks, num_blocks)
    )

    for c1 in 1:num_blocks
        m = Lp_data[:Lpy_comp][c1]

        if isempty(m)
            Lp_data[:Lpy][:D][c1] = Array{Float64}(undef, 0, 0)
        else
            Lp_data[:Lpy][:D][c1] = calcola_Lp_Sym(m, volumi, incidence_selection, escalings, id, "y")
        end
    end

    for c1 in 1:num_blocks-1
        for c2 in c1+1:num_blocks
            if !isempty(Lp_data[:Lpy_comp][c1]) && !isempty(Lp_data[:Lpy_comp][c2])
                Lp_data[:Lpy][:U][c1, c2], Lp_data[:Lpy][:V][c1, c2] = ACA_Lp(ACA_thres, Lp_data[:Lpy_comp][c1], Lp_data[:Lpy_comp][c2], volumi, incidence_selection, escalings, id, use_suppression, "y")
            end
        end
    end

    # --- Lpz Matrix elements ---
    # estremi_celle = volumi[:coordinate][mx+my+1:mx+my+mz, :]
    # versori = hcat(zeros(mz, 2), ones(mz, 1))
    # centri = volumi[:centri][mx+my+1:mx+my+mz, :]

    num_blocks = length(Lp_data[:Lpz_comp])
    Lp_data[:Lpz] = Dict{Symbol,Any}(
        :D => Vector{Matrix{Float64}}(undef, num_blocks),
        :U => Matrix{Union{Nothing,Matrix{Float64}}}(nothing, num_blocks, num_blocks),
        :V => Matrix{Union{Nothing,Matrix{Float64}}}(nothing, num_blocks, num_blocks)
    )

    for c1 in 1:num_blocks
        m = Lp_data[:Lpz_comp][c1]

        if isempty(m)
            Lp_data[:Lpz][:D][c1] = Array{Float64}(undef, 0, 0)
        else
            Lp_data[:Lpz][:D][c1] = calcola_Lp_Sym(m, volumi, incidence_selection, escalings, id, "z")
        end
    end

    for c1 in 1:num_blocks-1
        for c2 in c1+1:num_blocks
            if !isempty(Lp_data[:Lpz_comp][c1]) && !isempty(Lp_data[:Lpz_comp][c2])
                Lp_data[:Lpz][:U][c1, c2], Lp_data[:Lpz][:V][c1, c2] = ACA_Lp(ACA_thres, Lp_data[:Lpz_comp][c1], Lp_data[:Lpz_comp][c2], volumi, incidence_selection, escalings, id, use_suppression, "z")
            end
        end
    end

    return Lp_data

end
