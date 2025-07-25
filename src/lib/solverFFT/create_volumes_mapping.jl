function create_volumes_mapping(grids)
    num_grids = length(grids)
    #println(size(grids))
    #Nx, Ny, Nz = size(grids[1])
    Nx = size(grids[1], 1)
    Ny = size(grids[1][1], 1)
    Nz = size(grids[1][1][1], 1)
    mapping = zeros(Nx * Ny * Nz)
    num_ele = 0
    for cont3 = 1:Nz
        for cont2 = 1:Ny
            for cont = 1:Nx
                for k = 1:num_grids
                    if grids[k][cont][cont2][cont3] != 0
                        num_ele += 1
                        mapping[from_3D_to_1D(cont, cont2, cont3, Nx, Ny)] = num_ele
                        break
                    end
                end
            end
        end
    end
    return mapping, num_ele
end