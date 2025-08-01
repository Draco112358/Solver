function build_Yle_FFT(lumped_elements, grounding_nodes, ports, escalings, n, w, val_chiusura)
    le_nodes = vcat(lumped_elements["le_nodes"][:, 1], lumped_elements["le_nodes"][:, 2])
    port_nodes = vcat(ports["port_nodes"][:, 1], ports["port_nodes"][:, 2])
    N_ele = length(unique(vcat(le_nodes, port_nodes, grounding_nodes)))
    NNz_max = N_ele^2
    ind_r = zeros(NNz_max)
    ind_c = zeros(NNz_max)
    vals = zeros(ComplexF64, NNz_max)
    nlum = size(lumped_elements["le_nodes"], 1)
    cont = 0
    for c1 in range(1, nlum)
        n1 = lumped_elements["le_nodes"][c1, 1]
        n2 = lumped_elements["le_nodes"][c1, 2]
        ind1 = findall(ind_r .== n1)
        ind2 = findall(ind_c .== n1)
        ind = intersect(ind1, ind2)
        val_le=get_lumped_elements_admittance_fft(lumped_elements["type"][c1],lumped_elements["R_value"][c1],lumped_elements["L_value"][c1],lumped_elements["C_value"][c1],w);
        if isempty(ind)
            cont += 1
            ind_r[cont] = n1
            ind_c[cont] = n1
            vals[cont] = val_le
        else
            vals[ind[1]]=vals[ind[1]]+val_le
        end
        ind1 = findall(ind_r .== n2)
        ind2 = findall(ind_c .== n2)
        ind = intersect(ind1, ind2)
        if isempty(ind)
            cont += 1
            ind_r[cont] = n2
            ind_c[cont] = n2
            vals[cont] = val_le
        else
            vals[ind[1]]=vals[ind[1]]+val_le
        end
        ind1 = findall(ind_r .== n1)
        ind2 = findall(ind_c .== n2)
        ind = intersect(ind1, ind2)
        if isempty(ind)
            cont += 1
            ind_r[cont] = n1
            ind_c[cont] = n2
            vals[cont] = -val_le
        else
            vals[ind[1]]=vals[ind[1]]-val_le
        end
        ind1 = findall(ind_r .== n2)
        ind2 = findall(ind_c .== n1)
        ind = intersect(ind1, ind2)
        if isempty(ind)
            cont += 1
            ind_r[cont] = n2
            ind_c[cont] = n1
            vals[cont] = -val_le
        else
            vals[ind[1]]=vals[ind[1]]-val_le
        end
    end
    nlum = size(ports["port_nodes"], 1)
    for c1 in 1:nlum
        n1 = ports["port_nodes"][c1, 1]
        n2 = ports["port_nodes"][c1, 2]
        ind1 = findall(ind_r .== n1)
        ind2 = findall(ind_c .== n1)
        ind = intersect(ind1, ind2)
        if isempty(ind)
            cont += 1
            ind_r[cont] = n1
            ind_c[cont] = n1
            vals[cont] = 1 / val_chiusura
        else
            vals[ind[1]] += 1 / val_chiusura
        end
        ind1 = findall(ind_r .== n2)
        ind2 = findall(ind_c .== n2)
        ind = intersect(ind1, ind2)
        if isempty(ind)
            cont += 1
            ind_r[cont] = n2
            ind_c[cont] = n2
            vals[cont] = 1 / val_chiusura
        else
            vals[ind[1]] += 1 / val_chiusura
        end
        ind1 = findall(ind_r .== n1)
        ind2 = findall(ind_c .== n2)
        ind = intersect(ind1, ind2)
        if isempty(ind)
            cont += 1
            ind_r[cont] = n1
            ind_c[cont] = n2
            vals[cont] = -1 / val_chiusura
        else
            vals[ind[1]] -= 1 / val_chiusura
        end
        ind1 = findall(ind_r .== n2)
        ind2 = findall(ind_c .== n1)
        ind = intersect(ind1, ind2)
        if isempty(ind)
            cont += 1
            ind_r[cont] = n2
            ind_c[cont] = n1
            vals[cont] = -1 / val_chiusura
        else
            vals[ind[1]] -= 1 / val_chiusura
        end
    end
    for k in range(1, length(grounding_nodes))
        n1 = grounding_nodes[k]
        ind1 = findall(ind_r .== n1)
        ind2 = findall(ind_c .== n1)
        ind = intersect(ind1, ind2)
        if isempty(ind)
            cont += 1
            vals[cont] = 1 / 1e12
            ind_r[cont] = grounding_nodes[k]
            ind_c[cont] = grounding_nodes[k]
        else
            vals[ind[1]] += 1 / 1e12
        end
    end

    Yle = sparse(ind_r[1:cont], ind_c[1:cont], vals[1:cont] * escalings["Yle"], n, n)
    return Yle
end

function get_lumped_elements_admittance_fft(type, R, L, C, w)
    y = 0.0
    if (type == 1)
        z = R + 1im * w * L + 1 / (1im * w * C)
        if isnan(z) || isinf(z) || abs(z) < 1e-16
            z = R + 1im * w * L
        end
        y = 1 / z
    elseif (type == 2)
        y = (1 / (R + 1im * w * L) + 1im * w * C)
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = (1 / R + 1im * w * C)
        end
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = (1 / (R + 1im * w * L))
        end
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = 1 / R
        end
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = (1 / (1im * w * L) + 1im * w * C)
        end
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = 1im * w * C
        end
    elseif (type == 3)
        y = (1 / (R + 1 / (1im * w * C)) + 1 / (1im * w * L))
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = (1 / (R + 1 / (1im * w * C)))
        end
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = (1 / R + 1 / (1im * w * L))
        end
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = 1 / R
        end
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = (1 / (1 / (1im * w * C)) + 1 / (1im * w * L))
        end
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = (1 / (1im * w * L))
        end
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = (1 / (1 / (1im * w * C)))
        end
    elseif (type == 4)
        y = 1 / R + 1im * w * C + 1 / (1im * w * L)
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = 1 / R + 1im * w * C
        end
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = 1im * w * C + 1 / (1im * w * L)
        end
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = 1im * w * C
        end
    elseif (type == 5)
        y = 1 / R + (1 / (1im * w * L + 1 / (1im * w * C)))
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = 1 / R + (1 / (1 / (1im * w * C)))
        end
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = 1 / R + (1 / (1im * w * L))
        end
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = 1 / R
        end
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = (1 / (1im * w * L + 1 / (1im * w * C)))
        end
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = (1 / (1im * w * L))
        end
        if isnan(y) || isinf(y) || abs(y) < 1e-16
            y = (1 / (1 / (1im * w * C)))
        end
    end

    if isnan(y) || isinf(y)
        y = 0.0
    end
    return y
end