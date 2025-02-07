using SparseArrays

function build_Yle_S_new(lumped_elements, grounding_nodes, ports, escalings, n, w, val_chiusura, 
                         type_le, R_le, L_le, C_le)
    # Extract node information
    le_nodes = vcat(lumped_elements[:le_nodes][:, 1], lumped_elements[:le_nodes][:, 2])
    port_nodes = vcat(ports[:port_nodes][:, 1], ports[:port_nodes][:, 2])
    N_ele = length(unique(vcat(le_nodes, port_nodes, grounding_nodes)))
    
    # Initialize arrays
    NNz_max = N_ele^2
    ind_r = zeros(Int, NNz_max)
    ind_c = zeros(Int, NNz_max)
    vals = zeros(Float64, NNz_max)
    cont = 0
    
    # Process lumped elements
    for c1 in 1:size(lumped_elements[:le_nodes], 1)
        n1 = lumped_elements[:le_nodes][c1, 1]
        n2 = lumped_elements[:le_nodes][c1, 2]
        
        val_le = get_lumped_elements_admittance(type_le[c1], R_le[c1], L_le[c1], C_le[c1], w)
        
        # Add to diagonal for n1
        ind = findall(x -> x == n1, ind_r[1:cont]) ∩ findall(x -> x == n1, ind_c[1:cont])
        if isempty(ind)
            cont += 1
            ind_r[cont] = n1
            ind_c[cont] = n1
            vals[cont] = val_le
        else
            vals[ind[1]] += val_le
        end
        
        # Add to diagonal for n2
        ind = findall(x -> x == n2, ind_r[1:cont]) ∩ findall(x -> x == n2, ind_c[1:cont])
        if isempty(ind)
            cont += 1
            ind_r[cont] = n2
            ind_c[cont] = n2
            vals[cont] = val_le
        else
            vals[ind[1]] += val_le
        end
        
        # Add to off-diagonal n1-n2
        ind = findall(x -> x == n1, ind_r[1:cont]) ∩ findall(x -> x == n2, ind_c[1:cont])
        if isempty(ind)
            cont += 1
            ind_r[cont] = n1
            ind_c[cont] = n2
            vals[cont] = -val_le
        else
            vals[ind[1]] -= val_le
        end
        
        # Add to off-diagonal n2-n1
        ind = findall(x -> x == n2, ind_r[1:cont]) ∩ findall(x -> x == n1, ind_c[1:cont])
        if isempty(ind)
            cont += 1
            ind_r[cont] = n2
            ind_c[cont] = n1
            vals[cont] = -val_le
        else
            vals[ind[1]] -= val_le
        end
    end
    
    # Process ports
    for c1 in 1:size(ports[:port_nodes], 1)
        n1 = ports[:port_nodes][c1, 1]
        n2 = ports[:port_nodes][c1, 2]
        
        # Add to diagonal for n1
        ind = findall(x -> x == n1, ind_r[1:cont]) ∩ findall(x -> x == n1, ind_c[1:cont])
        if isempty(ind)
            cont += 1
            ind_r[cont] = n1
            ind_c[cont] = n1
            vals[cont] = 1 / val_chiusura
        else
            vals[ind[1]] += 1 / val_chiusura
        end
        
        # Add to diagonal for n2
        ind = findall(x -> x == n2, ind_r[1:cont]) ∩ findall(x -> x == n2, ind_c[1:cont])
        if isempty(ind)
            cont += 1
            ind_r[cont] = n2
            ind_c[cont] = n2
            vals[cont] = 1 / val_chiusura
        else
            vals[ind[1]] += 1 / val_chiusura
        end
        
        # Add to off-diagonal n1-n2
        ind = findall(x -> x == n1, ind_r[1:cont]) ∩ findall(x -> x == n2, ind_c[1:cont])
        if isempty(ind)
            cont += 1
            ind_r[cont] = n1
            ind_c[cont] = n2
            vals[cont] = -1 / val_chiusura
        else
            vals[ind[1]] -= 1 / val_chiusura
        end
        
        # Add to off-diagonal n2-n1
        ind = findall(x -> x == n2, ind_r[1:cont]) ∩ findall(x -> x == n1, ind_c[1:cont])
        if isempty(ind)
            cont += 1
            ind_r[cont] = n2
            ind_c[cont] = n1
            vals[cont] = -1 / val_chiusura
        else
            vals[ind[1]] -= 1 / val_chiusura
        end
    end
    
    # Process grounding nodes
    for k in grounding_nodes
        ind = findall(x -> x == k, ind_r[1:cont]) ∩ findall(x -> x == k, ind_c[1:cont])
        if isempty(ind)
            cont += 1
            ind_r[cont] = k
            ind_c[cont] = k
            vals[cont] = 1 / 1e12
        else
            vals[ind[1]] += 1 / 1e12
        end
    end
    
    # Build the sparse matrix
    Yle = sparse(ind_r[1:cont], ind_c[1:cont], vals[1:cont] * escalings[:Yle], n, n)
    return Yle
end



function get_lumped_elements_admittance(type, R, L, C, w)
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
