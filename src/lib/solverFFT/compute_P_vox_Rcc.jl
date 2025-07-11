function compute_P_vox_Rcc(centers_o, centers, sx, sy, sz, caser, caser2, id)
    # Define xc1, yc1, zc1 based on the value of `caser`
    if caser == 1
        xc1 = centers_o[1]
        yc1 = centers_o[2] - sy / 2
        zc1 = centers_o[3]
    elseif caser == 2
        xc1 = centers_o[1]
        yc1 = centers_o[2] + sy / 2
        zc1 = centers_o[3]
    elseif caser == 3
        xc1 = centers_o[1] - sx / 2
        yc1 = centers_o[2]
        zc1 = centers_o[3]
    elseif caser == 4
        xc1 = centers_o[1] + sx / 2
        yc1 = centers_o[2]
        zc1 = centers_o[3]
    elseif caser == 5
        xc1 = centers_o[1]
        yc1 = centers_o[2]
        zc1 = centers_o[3] - sz / 2
    else
        xc1 = centers_o[1]
        yc1 = centers_o[2]
        zc1 = centers_o[3] + sz / 2
    end

    N = size(centers, 1)
    Rcc = zeros(N)

    # Compute Rcc based on the value of `caser2`
    if caser2 == 1
        for k in 1:N
            Rcc[k] = norm([xc1, yc1, zc1] .- [centers[k, 1], centers[k, 2] - sy / 2, centers[k, 3]])
            if is_stop_requested(id)
                println("Simulazione $(id) interrotta per richiesta stop.")
                return nothing # O un altro valore che indica interruzione
            end
        end
    elseif caser2 == 2
        for k in 1:N
            Rcc[k] = norm([xc1, yc1, zc1] .- [centers[k, 1], centers[k, 2] + sy / 2, centers[k, 3]])
            if is_stop_requested(id)
                println("Simulazione $(id) interrotta per richiesta stop.")
                return nothing # O un altro valore che indica interruzione
            end
        end
    elseif caser2 == 3
        for k in 1:N
            Rcc[k] = norm([xc1, yc1, zc1] .- [centers[k, 1] - sx / 2, centers[k, 2], centers[k, 3]])
            if is_stop_requested(id)
                println("Simulazione $(id) interrotta per richiesta stop.")
                return nothing # O un altro valore che indica interruzione
            end
        end
    elseif caser2 == 4
        for k in 1:N
            Rcc[k] = norm([xc1, yc1, zc1] .- [centers[k, 1] + sx / 2, centers[k, 2], centers[k, 3]])
            if is_stop_requested(id)
                println("Simulazione $(id) interrotta per richiesta stop.")
                return nothing # O un altro valore che indica interruzione
            end
        end
    elseif caser2 == 5
        for k in 1:N
            Rcc[k] = norm([xc1, yc1, zc1] .- [centers[k, 1], centers[k, 2], centers[k, 3] - sz / 2])
            if is_stop_requested(id)
                println("Simulazione $(id) interrotta per richiesta stop.")
                return nothing # O un altro valore che indica interruzione
            end
        end
    else
        for k in 1:N
            Rcc[k] = norm([xc1, yc1, zc1] .- [centers[k, 1], centers[k, 2], centers[k, 3] + sz / 2])
            if is_stop_requested(id)
                println("Simulazione $(id) interrotta per richiesta stop.")
                return nothing # O un altro valore che indica interruzione
            end
        end
    end

    return Rcc
end
