function Compute_Lp_Self(bar, curr_dir)
    if curr_dir == 1
        l = abs(mean_length_rev(bar, 1))
        W = abs(mean_length_rev(bar, 3))
        T = abs(mean_length_rev(bar, 2))
    elseif curr_dir == 2
        T = abs(mean_length_rev(bar, 1))
        W = abs(mean_length_rev(bar, 3))
        l = abs(mean_length_rev(bar, 2))
    else
        W = abs(mean_length_rev(bar, 1))
        l = abs(mean_length_rev(bar, 3))
        T = abs(mean_length_rev(bar, 2))
    end

    Lp_Self_Rect = 0.0

    if !isnan(1 / (l * W * T)) && !isinf(1 / (l * W * T))
        # Fast Henry
        w = W / l
        t = T / l
        r = sqrt(w^2 + t^2)
        aw = sqrt(w^2 + 1)
        at = sqrt(t^2 + 1)
        ar = sqrt(w^2 + t^2 + 1)

        mu0 = 4 * π * 1e-7

        Lp_Self_Rect = 2 * mu0 * l / π * (
            1 / 4 * (1 / w * asinh(w / at) + 1 / t * asinh(t / aw) + asinh(1 / r)) +
            1 / 24 * (
                t^2 / w * asinh(w / (t * at * (r + ar))) +
                w^2 / t * asinh(t / (w * aw * (r + ar))) +
                t^2 / w^2 * asinh(w^2 / (t * r * (at + ar))) +
                w^2 / t^2 * asinh(t^2 / (w * r * (aw + ar))) +
                1 / (w * t^2) * asinh(w * t^2 / (at * (aw + ar))) +
                1 / (t * w^2) * asinh(t * w^2 / (aw * (at + ar)))
            ) -
            1 / 6 * (1 / (w * t) * atan(w * t / ar) + t / w * atan(w / (t * ar)) + w / t * atan(t / (w * ar))) -
            1 / 60 * (
                (ar + r + t + at) * t^2 / ((ar + r) * (r + t) * (t + at) * (at + ar)) +
                (ar + r + w + aw) * w^2 / ((ar + r) * (r + w) * (w + aw) * (aw + ar)) +
                (ar + aw + 1 + at) / ((ar + aw) * (aw + 1) * (1 + at) * (at + ar))
            ) -
            1 / 20 * (1 / (r + ar) + 1 / (aw + ar) + 1 / (at + ar))
        )
    end

    return Lp_Self_Rect
end

function mean_length_rev(barra1, curr_dir)
    # Computes the mean cross-section of the PSP described by barra1 and curr_dir
    xi1 = barra1[[1, 4, 7, 10]]
    yi1 = barra1[[2, 5, 8, 11]]
    zi1 = barra1[[3, 6, 9, 12]]
    xi2 = barra1[[13, 16, 19, 22]]
    yi2 = barra1[[14, 17, 20, 23]]
    zi2 = barra1[[15, 18, 21, 24]]

    ri = zeros(8, 3)

    # Vectors pointing to the vertices of the PSP
    ri[1, :] = [xi1[1], yi1[1], zi1[1]]
    ri[2, :] = [xi1[2], yi1[2], zi1[2]]
    ri[3, :] = [xi1[3], yi1[3], zi1[3]]
    ri[4, :] = [xi1[4], yi1[4], zi1[4]]
    ri[5, :] = [xi2[1], yi2[1], zi2[1]]
    ri[6, :] = [xi2[2], yi2[2], zi2[2]]
    ri[7, :] = [xi2[3], yi2[3], zi2[3]]
    ri[8, :] = [xi2[4], yi2[4], zi2[4]]

    if curr_dir == 1
        r1 = 0.25 * (ri[1, :] + ri[3, :] + ri[5, :] + ri[7, :])
        r2 = 0.25 * (ri[2, :] + ri[4, :] + ri[6, :] + ri[8, :])
    elseif curr_dir == 2
        r1 = 0.25 * (ri[1, :] + ri[2, :] + ri[5, :] + ri[6, :])
        r2 = 0.25 * (ri[3, :] + ri[4, :] + ri[7, :] + ri[8, :])
    else
        r1 = 0.25 * (ri[1, :] + ri[2, :] + ri[3, :] + ri[4, :])
        r2 = 0.25 * (ri[5, :] + ri[6, :] + ri[7, :] + ri[8, :])
    end

    return norm(r1 - r2, 2)
end
