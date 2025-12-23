function Compute_Lp_Self(l::Float64, W::Float64, T::Float64)

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
