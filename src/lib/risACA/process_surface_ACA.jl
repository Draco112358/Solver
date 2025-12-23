function process_surface_ACA!(
    x_vs::Vector{Vector{Float64}}, y_vs::Vector{Vector{Float64}}, z_vs::Vector{Vector{Float64}},
    xc_s::Vector{Float64}, yc_s::Vector{Float64}, zc_s::Vector{Float64},
    a_s::Vector{Float64}, b_s::Vector{Float64}, c_s::Vector{Float64}, sup_s::Vector{Float64},
    sup_xz_planes::Vector{Int}, sup_yz_planes::Vector{Int},
    superfici::Dict{String,Any}, round_precision::Int, ind::Vector{Int}
)
    nsup = length(ind)

    Threads.@threads for i in 1:nsup
        # Using a direct view is good, but for small, fixed-size access,
        # direct indexing might be slightly faster due to avoiding view object creation.
        # However, @view is generally robust. Let's stick with it for clarity.
        current_surface = @view superfici["estremi_celle"][ind[i], :]

        # Extract values directly to avoid re-indexing `current_surface` multiple times
        x1 = current_surface[1]
        x10 = current_surface[10]
        y2 = current_surface[2]
        y11 = current_surface[11]
        z3 = current_surface[3]
        z12 = current_surface[12]

        # For fixed-size vectors (length 2), it's often more efficient to create them directly
        # than rely on `unique` which is more general.
        # Assuming `x_vs[i]` etc. are already `Vector{Float64}`.
        x_vs[i] = [round(min(x1, x10), digits=round_precision), round(max(x1, x10), digits=round_precision)]
        y_vs[i] = [round(min(y2, y11), digits=round_precision), round(max(y2, y11), digits=round_precision)]
        z_vs[i] = [round(min(z3, z12), digits=round_precision), round(max(z3, z12), digits=round_precision)]

        # Calculate centroids directly from the min/max values
        xc_s[i] = 0.5 * (x_vs[i][2] + x_vs[i][1]) # Assuming min/max order
        yc_s[i] = 0.5 * (y_vs[i][2] + y_vs[i][1])
        zc_s[i] = 0.5 * (z_vs[i][2] + z_vs[i][1])

        # Calculate absolute differences for side lengths
        a_val = abs(x_vs[i][2] - x_vs[i][1])
        b_val = abs(y_vs[i][2] - y_vs[i][1])
        c_val = abs(z_vs[i][2] - z_vs[i][1])

        sup1_yz_plane = 0
        sup1_xz_plane = 0

        # Use `isapprox` for floating-point comparisons if `a_val`, `b_val`, `c_val`
        # could be very close due to `round_precision`. For exact comparisons, `<= ` is fine.
        # Assuming exact comparisons are intended here.
        if (a_val <= b_val && a_val <= c_val)
            sup1_yz_plane = 1
            a_val = 1.0 # Use 1.0 for Float64 literal
        elseif (b_val <= a_val && b_val <= c_val)
            sup1_xz_plane = 1
            b_val = 1.0
        else
            c_val = 1.0
        end

        # Store results in the pre-allocated arrays
        a_s[i] = a_val
        b_s[i] = b_val
        c_s[i] = c_val
        sup_s[i] = a_val * b_val * c_val # This might not be a "surface" area but a volume if dimensions are set to 1.0
        sup_xz_planes[i] = sup1_xz_plane
        sup_yz_planes[i] = sup1_yz_plane
    end
    return nothing
end