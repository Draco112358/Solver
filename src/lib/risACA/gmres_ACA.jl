# GMRES Custom Adapted
function gmres_ACA(b, restarted, tol, maxit, x, wk, incidence_selection, P_data, Lp_data, Z_self, Yle, F, resProd, id, chan, portIndex)

    m_size = size(b, 1)
    n = m_size

    # Handling logical inputs for restarted/maxit
    outer = 1
    inner = maxit

    n2b = norm(b)
    if n2b == 0
        return zeros(ComplexF64, n), 0, 0.0, [0, 0], [0.0]
    end

    if isempty(x)
        x = zeros(ComplexF64, n)
    end

    flag = 1
    xmin = x
    minupdated = 0

    tolb = tol * n2b

    # Initial residual
    # ComputeMatrixVectorACA(x, w, incidence_selection, P_data, Lp_data, Z_self, Yle, F, resProd)
    r = b - ComputeMatrixVectorACA(x, wk, incidence_selection, P_data, Lp_data, Z_self, Yle, F, resProd)
    normr = norm(r)

    if normr <= tolb
        return x, 0, normr / n2b, [0, 0], [normr]
    end

    minv_b = b # In this code b is already preconditioned? 
    # Wait, in MATLAB: tn=Q1*(U1\(L1\(P1*[zeros;is]))); This is the preconditioned RHS.
    # So 'b' passed here is indeed preconditioned.
    # ComputeMatrixVectorACA applies the preconditioner F \ (...) at the end.
    # So we are solving preconditioned system A' x = b'

    resvec = zeros(Float64, inner * outer + 1)
    resvec[1] = normr
    normrmin = normr
    J = zeros(ComplexF64, 2, inner)
    U = zeros(ComplexF64, n, inner + 1)
    R = zeros(ComplexF64, inner, inner)
    w = zeros(ComplexF64, inner + 1)

    outitercount = 0
    initercount = 0

    # Householder for r
    u = copy(r)
    beta_h = scalarsignRisACA(r[1]) * normr
    u[1] += beta_h
    u /= norm(u)

    U[:, 1] = copy(u)
    w[1] = -beta_h

    for initer = 1:inner
        initercount = initer

        # Pj * ej
        v = -2 * conj(u[initer]) * u
        v[initer] += 1
        for k = initer-1:-1:1
            v .-= U[:, k] * (2 * dot(U[:, k], v))
        end
        v /= norm(v)

        # Apply A to v
        v = ComputeMatrixVectorACA(v, wk, incidence_selection, P_data, Lp_data, Z_self, Yle, F, resProd)

        # Pj...P1 * Av
        for k = 1:initer
            v .-= U[:, k] * (2 * dot(U[:, k], v))
        end

        if initer != length(v)
            u = copy(v)
            u[1:initer] .= 0
            alpha = norm(u)
            if alpha != 0
                alpha = scalarsignRisACA(v[initer+1]) * alpha
                u[initer+1] += alpha
                u /= norm(u)
                U[:, initer+1] = copy(u)
                v[initer+2:end] .= 0
                v[initer+1] = -alpha
            end
        end

        # Givens
        for colJ = 1:initer-1
            tmpv = v[colJ]
            v[colJ] = conj(J[1, colJ]) * v[colJ] + conj(J[2, colJ]) * v[colJ+1]
            v[colJ+1] = -J[2, colJ] * tmpv + J[1, colJ] * v[colJ+1]
        end

        # Compute Givens Jm
        if !(initer == length(v))
            rho = norm(v[initer:initer+1])
            J[:, initer] = v[initer:initer+1] ./ rho
            w[initer+1] = -J[2, initer] * w[initer]
            w[initer] = conj(J[1, initer]) * w[initer]
            v[initer] = rho
            v[initer+1] = 0
        end

        R[:, initer] = v[1:inner]
        normr = abs(w[initer+1])
        resvec[initer+1] = normr

        if normr <= tolb
            flag = 0
            # update x
            y = R[1:initer, 1:initer] \ w[1:initer]
            # additive = P1...Pj * y
            # logic slightly complex to unwind Householder:
            # Just use the xm logic from custom or standard GMRES.
            # Simplified reconstruction:
            additive = zeros(ComplexF64, n)
            additive[1:initer] = y
            for k = initer:-1:1
                additive .-= U[:, k] * (2 * dot(U[:, k], additive))
            end
            x += additive
            break
        end
    end

    if flag != 0
        y = R[1:initercount, 1:initercount] \ w[1:initercount]
        additive = zeros(ComplexF64, n)
        additive[1:initercount] = y
        for k = initercount:-1:1
            additive .-= U[:, k] * (2 * dot(U[:, k], additive))
        end
        x += additive

        # Final residual check
        r = b - ComputeMatrixVectorACA(x, wk, incidence_selection, P_data, Lp_data, Z_self, Yle, F, resProd)
        normr = norm(r)

        if normr <= tolb
            flag = 0
        end
    end

    relres = normr / n2b
    return x, flag, relres, [1, initercount], resvec
end

function scalarsignRisACA(d)
    sgn = sign(d)
    if sgn == 0
        sgn = 1.0
    end
    return sgn
end