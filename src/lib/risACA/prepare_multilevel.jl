using Statistics

"""
    prepare_multilevel(centri, L)

Prepare multilevel basis function boxes for ACA algorithm.

This function creates a hierarchical spatial decomposition of basis functions
by recursively subdividing the bounding box along the largest dimension.

# Arguments
- `centri`: Matrix of center coordinates (N×3) where each row is [x, y, z]
- `L`: Initial level parameter for multilevel decomposition

# Returns
- `basis_func_boxes`: Matrix (N×Level) containing box indices for each basis function at each level
- `Level`: Actual number of levels created (may be less than L if subdivision cannot continue)
"""
function prepare_multilevel(centri, L)
    rcx = centri[:, 1]
    rcy = centri[:, 2]
    rcz = centri[:, 3]

    xmax = maximum(rcx)
    xmin = minimum(rcx)
    ymax = maximum(rcy)
    ymin = minimum(rcy)
    zmax = maximum(rcz)
    zmin = minimum(rcz)

    box_size_x = (xmax - xmin) * (1 + 1e-7)
    box_size_y = (ymax - ymin) * (1 + 1e-7)
    box_size_z = (zmax - zmin) * (1 + 1e-7)

    buffer = max(100, 2 * (2^L))  # Increased buffer size to handle subdivisions
    VtempCut = zeros(buffer)
    vCutx = zeros(buffer)
    vCuty = zeros(buffer)
    vCutz = zeros(buffer)

    vCutx[1] = xmin
    vCutx[2] = xmax
    vCuty[1] = ymin
    vCuty[2] = ymax
    vCutz[1] = zmin
    vCutz[2] = zmax

    lenx = 2
    leny = 2
    lenz = 2

    control = 0
    Level = L

    basis_func_boxes = zeros(size(rcx, 1), Level)

    for l in 1:L
        if control == 0
            # Subdivision process: divide box_size along the coordinate with maximum extent
            # based on the mean value of basis functions in that coordinate
            if max(box_size_x, box_size_y, box_size_z) == box_size_x
                VtempCut[1] = 0
                cont = 1

                for i in 1:(lenx-1)
                    sup_x = vCutx[i+1]
                    inf_x = vCutx[i]
                    res = (rcx .<= sup_x) .& (rcx .>= inf_x)
                    xavg = mean(rcx[res])

                    if xavg != sup_x && xavg != inf_x
                        if i == 1
                            VtempCut[cont] = inf_x
                            cont += 1
                            box_size_x = xavg - inf_x
                        end
                        VtempCut[cont] = xavg
                        if box_size_x < (xavg - inf_x)
                            box_size_x = xavg - inf_x
                        end
                        cont += 1
                        VtempCut[cont] = sup_x
                        cont += 1
                        if box_size_x < (sup_x - xavg)
                            box_size_x = sup_x - xavg
                        end
                    else
                        control = 1
                        Level = l - 1
                    end
                end
                lenx = cont - 1
                vCutx = zeros(length(vCutx))
                vCutx[1:lenx] = VtempCut[1:lenx]

            elseif max(box_size_x, box_size_y, box_size_z) == box_size_y
                VtempCut[1] = 0
                cont = 1

                for i in 1:(leny-1)
                    sup_y = vCuty[i+1]
                    inf_y = vCuty[i]
                    res = (rcy .<= sup_y) .& (rcy .>= inf_y)
                    yavg = mean(rcy[res])

                    if yavg != sup_y && yavg != inf_y
                        if i == 1
                            VtempCut[cont] = inf_y
                            cont += 1
                            box_size_y = yavg - inf_y
                        end
                        VtempCut[cont] = yavg

                        if box_size_y < (yavg - inf_y)
                            box_size_y = yavg - inf_y
                        end
                        cont += 1
                        VtempCut[cont] = sup_y
                        cont += 1
                        if box_size_y < (sup_y - yavg)
                            box_size_y = sup_y - yavg
                        end
                    else
                        control = 1
                        Level = l - 1
                    end
                end

                leny = cont - 1
                vCuty = zeros(length(vCuty))
                vCuty[1:leny] = VtempCut[1:leny]

            else  # max dimension is z
                VtempCut[1] = 0
                cont = 1

                for i in 1:(lenz-1)
                    sup_z = vCutz[i+1]
                    inf_z = vCutz[i]
                    res = (rcz .<= sup_z) .& (rcz .>= inf_z)
                    zavg = mean(rcz[res])

                    if zavg != sup_z && zavg != inf_z
                        if i == 1
                            VtempCut[cont] = inf_z
                            cont += 1
                            box_size_z = zavg - inf_z
                        end
                        VtempCut[cont] = zavg

                        if box_size_z < (zavg - inf_z)
                            box_size_z = zavg - inf_z
                        end

                        cont += 1
                        VtempCut[cont] = sup_z
                        cont += 1
                        if box_size_z < (sup_z - zavg)
                            box_size_z = sup_z - zavg
                        end
                    else
                        control = 1
                        Level = l - 1
                    end
                end

                lenz = cont - 1
                vCutz = zeros(length(vCutz))
                vCutz[1:lenz] = VtempCut[1:lenz]
            end

            # Assign box indices to basis functions
            if control == 0
                cl = 0
                for i in 1:(lenx-1)
                    sup_x = vCutx[i+1]
                    inf_x = vCutx[i]
                    for j in 1:(leny-1)
                        sup_y = vCuty[j+1]
                        inf_y = vCuty[j]
                        for k in 1:(lenz-1)
                            sup_z = vCutz[k+1]
                            inf_z = vCutz[k]
                            for h in range(1, length(rcx))
                                if (rcx[h] >= inf_x && rcx[h] <= sup_x) &&
                                   (rcy[h] >= inf_y && rcy[h] <= sup_y) &&
                                   (rcz[h] >= inf_z && rcz[h] <= sup_z)
                                    basis_func_boxes[h, l] = cl
                                end
                            end
                            cl += 1
                        end
                    end
                end
            end
        end
        VtempCut = zeros(length(VtempCut))
    end

    basis_func_boxes = basis_func_boxes[:, 1:Level]

    return basis_func_boxes, Level
end
