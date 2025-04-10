using LinearAlgebra

function compute_Ar_Gauss(barre::Matrix{Float64}, centriOss::Matrix{Float64}, ordine::Int, beta::Float64)
    numCentri = size(centriOss, 1)
    numBarre = size(barre, 1)

    ha = zeros(ComplexF64, numBarre, numCentri)

    for cont in 1:numBarre
        barra = barre[cont, :]

        xb = barra[[1, 4, 7, 10, 13, 16, 19, 22]]
        yb = barra[[2, 5, 8, 11, 14, 17, 20, 23]]
        zb = barra[[3, 6, 9, 12, 15, 18, 21, 24]]

        x_bar = [minimum(xb), maximum(xb)]
        y_bar = [minimum(yb), maximum(yb)]
        z_bar = [minimum(zb), maximum(zb)]

        for cc in 1:numCentri
            x_o = centriOss[cc, 1]
            y_o = centriOss[cc, 2]
            z_o = centriOss[cc, 3]

            ha[cont, cc] = compute_ha(x_o, x_bar, y_o, y_bar, z_o, z_bar, ordine, beta)
        end
    end

    return ha
end

function compute_ha(xo::Float64, x_vect_bar::Vector{Float64}, yo::Float64, y_vect_bar::Vector{Float64}, zo::Float64, z_vect_bar::Vector{Float64}, ordine::Int, beta::Float64)
    x1 = x_vect_bar[1]
    x2 = x_vect_bar[2]
    y1 = y_vect_bar[1]
    y2 = y_vect_bar[2]
    z1 = z_vect_bar[1]
    z2 = z_vect_bar[2]
    barra = [x1 y1 z1; x2 y1 z1; x1 y2 z1; x2 y2 z1; x1 y1 z2; x2 y1 z2; x1 y2 z2; x2 y2 z2]

    xi1 = barra[[1, 2, 3, 4], 1]
    yi1 = barra[[1, 2, 3, 4], 2]
    zi1 = barra[[1, 2, 3, 4], 3]
    xi2 = barra[[5, 6, 7, 8], 1]
    yi2 = barra[[5, 6, 7, 8], 2]
    zi2 = barra[[5, 6, 7, 8], 3]

    # vectors pointing to the vertices of the quadrilateral i
    ri = zeros(ComplexF64, 8, 3)
    ri[1, :] = [xi1[1], yi1[1], zi1[1]]
    ri[2, :] = [xi1[2], yi1[2], zi1[2]]
    ri[3, :] = [xi1[3], yi1[3], zi1[3]]
    ri[4, :] = [xi1[4], yi1[4], zi1[4]]

    ri[5, :] = [xi2[1], yi2[1], zi2[1]]
    ri[6, :] = [xi2[2], yi2[2], zi2[2]]
    ri[7, :] = [xi2[3], yi2[3], zi2[3]]
    ri[8, :] = [xi2[4], yi2[4], zi2[4]]

    # nuovo approccio
    rmi = 0.125 * sum(ri, dims=1)
    rai = 0.125 * (-ri[1, :] + ri[2, :] + ri[4, :] - ri[3, :] - ri[5, :] + ri[6, :] + ri[8, :] - ri[7, :])
    rbi = 0.125 * (-ri[1, :] - ri[2, :] + ri[4, :] + ri[3, :] - ri[5, :] - ri[6, :] + ri[8, :] + ri[7, :])
    rci = 0.125 * (-ri[1, :] - ri[2, :] - ri[4, :] - ri[3, :] + ri[5, :] + ri[6, :] + ri[8, :] + ri[7, :])
    rabi = 0.125 * (ri[1, :] - ri[2, :] + ri[4, :] - ri[3, :] + ri[5, :] - ri[6, :] + ri[8, :] - ri[7, :])
    rbci = 0.125 * (ri[1, :] + ri[2, :] - ri[4, :] - ri[3, :] - ri[5, :] - ri[6, :] + ri[8, :] + ri[7, :])
    raci = 0.125 * (ri[1, :] - ri[2, :] - ri[4, :] + ri[3, :] - ri[5, :] + ri[6, :] + ri[8, :] - ri[7, :])
    rabci = 0.125 * (-ri[1, :] + ri[2, :] - ri[4, :] + ri[3, :] + ri[5, :] - ri[6, :] + ri[8, :] - ri[7, :])

    rootkx, wekx = qrule(ordine)
    rootky, weky = qrule(ordine)
    rootkz, wekz = qrule(ordine)

    nlkx = length(wekx)
    nlky = length(weky)
    nlkz = length(wekz)

    sum_a1 = 0.0 + 0.0im
    for a1 in 1:nlkx
        sum_b1 = 0.0 + 0.0im
        for b1 in 1:nlky
            sum_c1 = 0.0 + 0.0im
            for c1 in 1:nlkz
                drai = rai + rabi * rootky[b1] + raci * rootkz[c1] + rabci * rootky[b1] * rootkz[c1]
                drbi = rbi + rabi * rootkx[a1] + rbci * rootkz[c1] + rabci * rootkx[a1] * rootkz[c1]
                drci = rci + raci * rootkx[a1] + rbci * rootky[b1] + rabci * rootkx[a1] * rootky[b1]
                draim = norm(drai)
                drbim = norm(drbi)
                drcim = norm(drci)

                r1 = rmi + rai * rootkx[a1] + rbi * rootky[b1] + rci * rootkz[c1] + rabi * rootkx[a1] * rootky[b1] + raci * rootkx[a1] * rootkz[c1] + rbci * rootky[b1] * rootkz[c1] +
                     rabci * rootkx[a1] * rootky[b1] * rootkz[c1]

                x = real(r1[1])
                y = real(r1[2])
                z = real(r1[3])

                delta_x = (xo - x)
                delta_y = (yo - y)
                delta_z = (zo - z)

                R = sqrt(delta_x^2 + delta_y^2 + delta_z^2)

                G = (1 / R) * exp(-1im * beta * R)

                f = draim * drbim * drcim * G

                sum_c1 += wekz[c1] * f
            end   # (c1)
            sum_b1 += weky[b1] * sum_c1
        end  # (b1)
        sum_a1 += wekx[a1] * sum_b1
    end # (a1)

    return sum_a1
end

function qrule(n::Int)
    iter = 2
    m = div(n + 1, 2)
    e1 = n * (n + 1)
    mm = 4 * m - 1
    t = (pi / (4 * n + 2)) * (3:4:mm)
    nn = (1 - (1 - 1 / n) / (8 * n * n))
    xo = nn * cos.(t)
    for kk in 1:iter
        pkm1 = zeros(size(xo))
        pkm1[1:size(xo, 1)] .= 1
        pk = xo
        for k in 2:n
            t1 = xo .* pk
            pkp1 = t1 - pkm1 - (t1 - pkm1) / k + t1
            pkm1 = pk
            pk = pkp1
        end
        den = 1 .- xo .^ 2
        d1 = n * (pkm1 - xo .* pk)
        dpn = d1 ./ den
        d2pn = (2 .* xo .* dpn - e1 .* pk) ./ den
        d3pn = (4 * xo .* d2pn + (2 - e1) .* dpn) ./ den
        d4pn = (6 * xo .* d3pn + (6 - e1) .* d2pn) ./ den
        u = pk ./ dpn
        v = d2pn ./ dpn
        h = -u .* (1 .+ (0.5 * u) .* (v + u .* (v .^ 2 - u .* d3pn ./ (3 * dpn))))
        p = pk + h .* (dpn + (0.5 * h) .* (d2pn + (h / 3) .* (d3pn + 0.25 * h .* d4pn)))
        dp = dpn + h .* (d2pn + (0.5 * h) .* (d3pn + h .* d4pn / 3))
        h = h - p ./ dp
        xo = xo + h
    end
    bp = zeros(Float64, n)
    wf = zeros(Float64, n)
    bp[1:size(xo, 1)] .= -xo .- h
    fx = d1 - h .* e1 .* (pk + (h / 2) .* (dpn + (h / 3) .* (
        d2pn + (h / 4) .* (d3pn + (0.2 * h) .* d4pn))))
    wf[1:size(xo, 1)] .= 2 * (1 .- bp[1:size(xo, 1)] .^ 2) ./ (fx .^ 2)
    if (m + m) > n
        bp[m] = 0.0
    end
    if !((m + m) == n)
        m = m - 1
    end
    jj = 1:m
    n1j = n + 1 .- jj
    bp[n1j] .= -bp[jj]
    wf[n1j] .= wf[jj]
    return bp, wf
end