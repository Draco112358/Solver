function build_Yle_S(lumped_elements, grounding_nodes, ports, escalings, n, w, val_chiusura, 
                         type_le, R_le, L_le, C_le)
    # Extract node information
    le_nodes = vcat(lumped_elements[:le_nodes][:, 1], lumped_elements[:le_nodes][:, 2])
    port_nodes = vcat(ports[:port_nodes][:, 1], ports[:port_nodes][:, 2])
    N_ele = length(unique(vcat(le_nodes, port_nodes, grounding_nodes)))
    
    # Initialize dictionaries for indices and values
    ind_dict = Dict{Tuple{Int, Int}, Float64}()
    
    # Process lumped elements
    for c1 in 1:size(lumped_elements[:le_nodes], 1)
        n1 = lumped_elements[:le_nodes][c1, 1]
        n2 = lumped_elements[:le_nodes][c1, 2]
        
        val_le = get_lumped_elements_admittance_S(type_le[c1], R_le[c1], L_le[c1], C_le[c1], w)
        
        # Add to diagonal for n1
        ind_dict[(n1, n1)] = get(ind_dict, (n1, n1), 0.0) + val_le
        
        # Add to diagonal for n2
        ind_dict[(n2, n2)] = get(ind_dict, (n2, n2), 0.0) + val_le
        
        # Add to off-diagonal n1-n2
        ind_dict[(n1, n2)] = get(ind_dict, (n1, n2), 0.0) - val_le
        
        # Add to off-diagonal n2-n1
        ind_dict[(n2, n1)] = get(ind_dict, (n2, n1), 0.0) - val_le
    end
    
    # Process ports
    for c1 in 1:size(ports[:port_nodes], 1)
        n1 = ports[:port_nodes][c1, 1]
        n2 = ports[:port_nodes][c1, 2]
        
        # Add to diagonal for n1
        ind_dict[(n1, n1)] = get(ind_dict, (n1, n1), 0.0) + 1 / val_chiusura
        
        # Add to diagonal for n2
        ind_dict[(n2, n2)] = get(ind_dict, (n2, n2), 0.0) + 1 / val_chiusura
        
        # Add to off-diagonal n1-n2
        ind_dict[(n1, n2)] = get(ind_dict, (n1, n2), 0.0) - 1 / val_chiusura
        
        # Add to off-diagonal n2-n1
        ind_dict[(n2, n1)] = get(ind_dict, (n2, n1), 0.0) - 1 / val_chiusura
    end
    
    # Process grounding nodes
    for k in grounding_nodes
        ind_dict[(k, k)] = get(ind_dict, (k, k), 0.0) + 1 / 1e12
    end
    
    # Convert dictionaries to arrays
    ind_r = Int[]
    ind_c = Int[]
    vals = Float64[]
    
    for ((r, c), v) in ind_dict
        push!(ind_r, r)
        push!(ind_c, c)
        push!(vals, v * escalings[:Yle])
    end
    
    # Build the sparse matrix
    Yle = sparse(ind_r, ind_c, vals, n, n)
    return Yle
end

function get_lumped_elements_admittance_S(type, R, L, C, w)
    y = 0.0
    if type == 1
        z = R + im*w*L + 1/(im*w*C)
        if isnan(z) || isinf(z) || abs(z) < 1e-16
            z = R + im*w*L
        end
        y = 1/z
    elseif type == 2
        y = (1/(R + im*w*L) + im*w*C)
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = (1/R + im*w*C)
        end
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = (1/(R + im*w*L))
        end
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = 1/R
        end
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = (1/(1/(im*w*C)) + im*w*C)
        end
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = im*w*C
        end
    elseif type == 3
        y = (1/(R + 1/(im*w*C)) + 1/(im*w*L))
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = (1/(R + 1/(im*w*C)))
        end
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = (1/R + 1/(im*w*L))
        end
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = 1/R
        end
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = (1/(1/(im*w*C)) + 1/(im*w*L))
        end
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = (1/(im*w*L))
        end
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = (1/(1/(im*w*C)))
        end
    elseif type == 4
        y = 1/R + im*w*C + 1/(im*w*L)
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = 1/R + im*w*C
        end
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = im*w*C + 1/(im*w*L)
        end
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = im*w*C
        end
    elseif type == 5
        y = 1/R + (1/(im*w*L + 1/(im*w*C)))
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = 1/R + (1/(1/(im*w*C)))
        end
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = 1/R + (1/(im*w*L))
        end
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = 1/R
        end
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = (1/(im*w*L + 1/(im*w*C)))
        end
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = (1/(im*w*L))
        end
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = (1/(1/(im*w*C)))
        end
    end

    if isnan(y) || isinf(y)
        y = 0
    end
    return y
end