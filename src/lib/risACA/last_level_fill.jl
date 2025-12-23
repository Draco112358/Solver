"""
    last_level_fill(basis_func_boxes, L)

Fill the last level of the multilevel structure by grouping basis functions
into boxes based on their assignment at level L.

# Arguments
- `basis_func_boxes`: Matrix containing box indices for each basis function at each level
- `L`: The level to use for grouping

# Returns
- `comp`: Vector of vectors, where each element contains the indices of basis functions
         belonging to that box at level L
"""
function last_level_fill(basis_func_boxes, L)
    num_col = 2^L
    comp = Vector{Vector{Int}}(undef, num_col)

    for c2 in 1:num_col
        # Find all basis functions that belong to box (c2-1) at level L
        comp[c2] = findall(basis_func_boxes[:, L] .== (c2 - 1))
    end

    return comp
end
