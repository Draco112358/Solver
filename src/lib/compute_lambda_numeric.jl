using LinearAlgebra

function compute_lambda_numeric(punti_oss, volumi, incidence_selection, vers_punti_oss, ordine_int, beta)

    N = size(volumi[:coordinate], 1)

    vers_celle = zeros(N, 3)
    vers_celle[1:incidence_selection[:mx], 1] .= 1
    vers_celle[incidence_selection[:mx]+1:incidence_selection[:mx]+incidence_selection[:my], 2] .= 1
    vers_celle[incidence_selection[:mx]+incidence_selection[:my]+1:N, 3] .= 1

    M = size(punti_oss, 1)

    Lambda = zeros(ComplexF64, M, N)

    rootkx, wekx = qrule(ordine_int)

    for m = 1:M
        for n = 1:N

            scelta = ""

            if abs(vers_celle[n, 1]) > 1e-5

                if abs(vers_punti_oss[m, 1]) > 1e-5
                    scelta = "xx"
                elseif abs(vers_punti_oss[m, 2]) > 1e-5
                    scelta = "yx"
                else
                    scelta = "zx"
                end

            elseif abs(vers_celle[n, 2]) > 1e-5

                if abs(vers_punti_oss[m, 1]) > 1e-5
                    scelta = "xy"
                elseif abs(vers_punti_oss[m, 2]) > 1e-5
                    scelta = "yy"
                else
                    scelta = "zy"
                end

            else
                if abs(vers_punti_oss[m, 1]) > 1e-5
                    scelta = "xz"
                elseif abs(vers_punti_oss[m, 2]) > 1e-5
                    scelta = "yz"
                else
                    scelta = "zz"
                end
            end

            Lambda[m, n] = compute_hi(volumi[:coordinate][n, :], punti_oss[m, :], scelta, rootkx, wekx, beta) / volumi[:S][n]

        end
    end

    return Lambda
end

function compute_hi(barra, centro_oss, scelta, rootkx, wekx, beta)

    xi1 = [barra[1], barra[4], barra[7], barra[10]]
    yi1 = [barra[2], barra[5], barra[8], barra[11]]
    zi1 = [barra[3], barra[6], barra[9], barra[12]]
    xi2 = [barra[13], barra[16], barra[19], barra[22]]
    yi2 = [barra[14], barra[17], barra[20], barra[23]]
    zi2 = [barra[15], barra[18], barra[21], barra[24]]

    #  vectors pointing to the vertices of the quadrilateral i
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
    rmi = vec(0.125 * sum(ri, dims=1))
    rai = 0.125 * (-ri[1, :] + ri[2, :] + ri[4, :] - ri[3, :] - ri[5, :] + ri[6, :] + ri[8, :] - ri[7, :])
    rbi = 0.125 * (-ri[1, :] - ri[2, :] + ri[4, :] + ri[3, :] - ri[5, :] - ri[6, :] + ri[8, :] + ri[7, :])
    rci = 0.125 * (-ri[1, :] - ri[2, :] - ri[4, :] - ri[3, :] + ri[5, :] + ri[6, :] + ri[8, :] + ri[7, :])
    rabi = 0.125 * (ri[1, :] - ri[2, :] + ri[4, :] - ri[3, :] + ri[5, :] - ri[6, :] + ri[8, :] - ri[7, :])
    rbci = 0.125 * (ri[1, :] + ri[2, :] - ri[4, :] - ri[3, :] - ri[5, :] - ri[6, :] + ri[8, :] + ri[7, :])
    raci = 0.125 * (ri[1, :] - ri[2, :] - ri[4, :] + ri[3, :] - ri[5, :] + ri[6, :] + ri[8, :] - ri[7, :])
    rabci = 0.125 * (-ri[1, :] + ri[2, :] - ri[4, :] + ri[3, :] + ri[5, :] - ri[6, :] + ri[8, :] - ri[7, :])

    nlkx = length(wekx)
    nlky = length(wekx)
    nlkz = length(wekx)

    sum_a1 = 0.0 + 0.0im
    for a1 = 1:nlkx
        sum_b1 = 0.0 + 0.0im
        for b1 = 1:nlky
            sum_c1 = 0.0 + 0.0im
            for c1 = 1:nlkz
                drai = rai + rabi * rootkx[b1] + raci * rootkx[c1] + rabci * rootkx[b1] * rootkx[c1]
                drbi = rbi + rabi * rootkx[a1] + rbci * rootkx[c1] + rabci * rootkx[a1] * rootkx[c1]
                drci = rci + raci * rootkx[a1] + rbci * rootkx[b1] + rabci * rootkx[a1] * rootkx[b1]
                draim = norm(drai, 2)
                drbim = norm(drbi, 2)
                drcim = norm(drci, 2)

                r1 = rmi + rai * rootkx[a1] + rbi * rootkx[b1] + rci * rootkx[c1] + rabi * rootkx[a1] * rootkx[b1] + raci * rootkx[a1] * rootkx[c1] + rbci * rootkx[b1] * rootkx[c1] +
                     rabci * rootkx[a1] * rootkx[b1] * rootkx[c1]

                x = real(r1[1])
                y = real(r1[2])
                z = real(r1[3])

                delta_x = (centro_oss[1] - x)
                delta_y = (centro_oss[2] - y)
                delta_z = (centro_oss[3] - z)

                R = sqrt(delta_x^2 + delta_y^2 + delta_z^2)

                G = if scelta == "xy"
                    delta_z * exp(-1im * beta * R) * (1 / R^3 + 1im * beta / R^2)
                elseif scelta == "xz"
                    -delta_y * exp(-1im * beta * R) * (1 / R^3 + 1im * beta / R^2)
                elseif scelta == "yx"
                    -delta_z * exp(-1im * beta * R) * (1 / R^3 + 1im * beta / R^2)
                elseif scelta == "yz"
                    delta_x * exp(-1im * beta * R) * (1 / R^3 + 1im * beta / R^2)
                elseif scelta == "zx"
                    delta_y * exp(-1im * beta * R) * (1 / R^3 + 1im * beta / R^2)
                elseif scelta == "zy"
                    -delta_x * exp(-1im * beta * R) * (1 / R^3 + 1im * beta / R^2)
                else
                    0.0 + 0.0im
                end

                f = draim * drbim * drcim * G

                sum_c1 += wekx[c1] * f
            end
            sum_b1 += wekx[b1] * sum_c1
        end
        sum_a1 += wekx[a1] * sum_b1
    end
    integral = 1e-7 * sum_a1
    return integral
end

function qrule(n::Int)
    iter = 2
    m = trunc((n + 1) / 2)
    e1 = n * (n + 1)
    mm = 4 * m - 1
    t = (pi / (4 * n + 2)) * (3:4:mm)
    nn = (1 - (1 - 1 / n) / (8 * n * n))
    xo = nn * cos.(t)
    den = []
    d1 = []
    dpn = []
    d2pn = []
    d3pn = []
    d4pn = []
    u = []
    v = []
    h = []
    p = []
    dp = []
    pk = []
    
    for kk = 1:iter
        pkm1 = zeros(size(xo))
        pkm1[1:size(xo, 1)] .= 1
        pk = xo
        for k = 2:n
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
        h = -u .* (1 .+ (0.5 * u) .* (v + u .* (v .* v - u .* d3pn ./ (3 * dpn))))
        p = pk + h .* (dpn + (0.5 * h) .* (d2pn + (h / 3) .* (d3pn + 0.25 * h .* d4pn)))
        dp = dpn + h .* (d2pn + (0.5 * h) .* (d3pn + h .* d4pn / 3))
        h = h - p ./ dp
        xo = xo + h
    end
    bp = zeros(1, n)
    wf = zeros(1, n)
    bp[1:size(xo, 1)] .= -xo .- h
    fx = d1 - h .* e1 .* (pk + (h / 2) .* (dpn + (h / 3) .* (
        d2pn + (h / 4) .* (d3pn + (0.2 * h) .* d4pn))))
    wf[1:size(xo, 1)] .= 2 * (1 .- bp[1:size(xo, 1)] .^ 2) ./ (fx .* fx)
    if (m + m) > n
        bp[Int64(m)] = 0
    end
    if !((m + m) == n)
        m = m - 1
    end
    jj = 1:m
    n1j = (n + 1) .- jj
    bp[Int64.(n1j)] .= -bp[Int64.(jj)]
    wf[Int64.(n1j)] .= wf[Int64.(jj)]

    return vec(bp), vec(wf)
end