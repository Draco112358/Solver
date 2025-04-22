function compute_E_field_Gauss(indx, indy, indz, centriOsservazione, hc, ha, J, sigma, f)

    numCentri = size(centriOsservazione, 1)

    mu0 = 4 * pi * 1e-7
    eps0 = 8.854187816997944e-12
    beta = 2 * pi * f * sqrt(eps0 * mu0)

    Ecx = zeros(ComplexF64, numCentri)
    Ecy = zeros(ComplexF64, numCentri)
    Ecz = zeros(ComplexF64, numCentri)
    for cont = 1:numCentri
        Ecx[cont] = sum(sigma .* squeeze(hc[:, 1, cont])) / (4 * pi * eps0)
        Ecy[cont] = sum(sigma .* squeeze(hc[:, 2, cont]) )/ (4 * pi * eps0)
        Ecz[cont] = sum(sigma .* squeeze(hc[:, 3, cont])) / (4 * pi * eps0)
    end

    Arx = zeros(ComplexF64, numCentri)
    Ary = zeros(ComplexF64, numCentri)
    Arz = zeros(ComplexF64, numCentri)

    for cont = 1:numCentri
        Arx[cont] = sum(J[indx] .* ha[indx, cont]) * 1e-7
        Ary[cont] = sum(J[indy] .* ha[indy, cont]) * 1e-7
        Arz[cont] = sum(J[indz] .* ha[indz, cont]) * 1e-7
    end

    Ex = -1im * 2 * pi * f * Arx + Ecx
    Ey = -1im * 2 * pi * f * Ary + Ecy
    Ez = -1im * 2 * pi * f * Arz + Ecz

    return Ex, Ey, Ez
end

function squeeze(A)
    return reshape(A, size(A)[[i for i in 1:ndims(A) if size(A, i) > 1]])  # Removes singleton dimensions
end