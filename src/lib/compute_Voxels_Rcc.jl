function compute_Voxels_Rcc(centri_m, centri_n)
    N = size(centri_n, 1)
    Rcc = zeros(N)

    c1 = centri_m

    for n in 1:N
        c2 = centri_n[n, :]
        Rcc[n] = norm(c1 - c2, 2)
    end

    return Rcc
end
