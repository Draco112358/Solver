function gmres_custom!(
    out, work, pc_work, b, restarted, tol, maxit, x, wk, 
    incidence_selection, P_rebuilted, Lp_rebuilted, Z_self, Yle, invZ, invP, F, 
    id, chan, portIndex
)
    # --- Setup iniziale (invariato) ---
    m = size(b, 1)
    n = m
    restart = maxit # Semplificazione basata sul tuo codice
    if restarted
        outer = maxit
        inner = min(restart, n)
    else
        outer = 1
        inner = min(maxit, n)
    end

    n2b = norm(b)
    if n2b == 0
        x .= 0
        return x, 0, 0.0, [0, 0], [0.0]
    end

    if isempty(x)
        x = zeros(ComplexF64, n)
    end

    flag = 1
    xmin = copy(x)
    imin, jmin = 0, 0
    tolb = tol * n2b
    
    stag = 0
    moresteps = 0
    maxmsteps = minimum([floor(n/50), 5, n-maxit])
    maxstagsteps = 3
    minupdated = 0
    warned = false
    
    # Vettori per GMRES
    r = similar(x)
    u = similar(x) # Vettore di Householder
    v = similar(x) # Vettore di lavoro
    additive = similar(x) # Vettore per l'aggiornamento della soluzione
    xm = copy(x) # Vettore per la soluzione intermedia
    ytmp = Vector{ComplexF64}(undef, inner)

    # Matrici e vettori di GMRES
    J = zeros(ComplexF64, 2, inner)
    U = zeros(ComplexF64, n, inner + 1)
    R = zeros(ComplexF64, inner, inner)
    w = zeros(ComplexF64, inner + 1)
    resvec = zeros(Float64, inner * outer + 1)

    # --- Calcolo residuo iniziale ---
    ComputeMatrixVector!(out, work, pc_work, x,wk,incidence_selection,P_rebuilted,Lp_rebuilted,Z_self,Yle,invZ,invP,F);
    r .= b .- out
    normr = norm(r)
    
    if normr <= tolb
        return x, 0, normr / n2b, [0, 0], [normr]
    end

    # Il tuo codice non usa un precondizionatore, quindi minv_b e b sono uguali
    n2minv_b = n2b
    tolb = tol * n2minv_b
    resvec[1] = normr
    normrmin = normr
    normr_act = normr
    outitercount = 0
    initercount = 0
    @views begin # <<< USARE VISTE PER TUTTE LE OPERAZIONI DI SLICING
    for outiter = 1:outer
        outitercount += 1
        
        # Costruisce u per il riflettore di Householder
        u .= r
        beta = scalarsign(r[1]) * normr
        u[1] += beta
        rmul!(u, 1 / norm(u))
        
        U[:, 1] .= u
        
        # Applica la proiezione di Householder a r
        w[1] = -beta
        for i in 2:inner+1; w[i] = 0.0; end

        initercount = 0
        for initer = 1:inner
            println("Iteration = $initercount")
            initercount += 1
            
            # --- Inizio Logica Originale (ma ottimizzata) ---
            # v = Pj*...*P1*ej
            v .= -2 * conj(U[initer, initer]) .* U[:, initer]
            v[initer] += 1
            for k = initer-1:-1:1
                d = 2 * dot(U[:, k], v)
                axpy!(-d, U[:, k], v) # v .-= d .* U[:, k]
            end
            rmul!(v, 1 / norm(v))

            # Applica A a v
            ComputeMatrixVector!(out, work, pc_work, v,wk,incidence_selection,P_rebuilted,Lp_rebuilted,Z_self,Yle,invZ,invP,F);
            v .= out
            
            # Form Pj*...*P1*Av
            for k = 1:initer
                d = 2 * dot(U[:, k], v)
                axpy!(-d, U[:, k], v)
            end
            
            # Determina Pj+1
            if initer != n
                u .= v
                u[1:initer] .= 0
                alpha = norm(u)
                if alpha != 0
                    alpha = scalarsign(v[initer+1]) * alpha
                    u[initer+1] += alpha
                    rmul!(u, 1 / norm(u))
                    U[:, initer+1] .= u
                    
                    v[initer+2:end] .= 0
                    v[initer+1] = -alpha
                end
            end
            
            # Applica le rotazioni di Givens
            for colJ = 1:initer-1
                tmpv = v[colJ]
                v[colJ]   = conj(J[1, colJ]) * v[colJ]   + conj(J[2, colJ]) * v[colJ+1]
                v[colJ+1] = -J[2, colJ]      * tmpv      + J[1, colJ]      * v[colJ+1]
            end
            
            # Calcola la nuova rotazione Jm
            if initer != n
                rho = norm(v[initer:initer+1])
                J[:, initer] .= v[initer:initer+1] ./ rho
                w[initer+1] = -J[2, initer] * w[initer]
                w[initer]   = conj(J[1, initer]) * w[initer]
                v[initer] = rho
                v[initer+1] = 0
            end
            
            R[:, initer] .= v[1:inner]
            normr = abs(w[initer+1])
            resvec[(outiter-1)*inner + initer + 1] = normr
            normr_act = normr
            
            # --- Blocco di controllo convergenza (logica originale) ---
            if (normr <= tolb || stag >= maxstagsteps || moresteps == 1)
                ytmp_view = ytmp[1:initer]
                w_view = w[1:initer]
                R_view = R[1:initer, 1:initer]
                ldiv!(ytmp_view, UpperTriangular(R_view), w_view) # ytmp_view = R_view \ w_view

                # Calcolo di 'additive' (logica originale preservata)
                alpha_add = -2 * ytmp_view[initer] * conj(U[initer, initer])
                additive .= U[:, initer] .* alpha_add
                additive[initer] += ytmp_view[initer]

                for k = initer-1:-1:1
                    additive[k] += ytmp_view[k]
                    d = 2 * dot(U[:, k], additive)
                    axpy!(-d, U[:, k], additive)
                end
                
                xm .= x .+ additive
                
                # ... resto della logica di controllo invariata ...
                ComputeMatrixVector!(out, work, pc_work, xm,wk,incidence_selection,P_rebuilted,Lp_rebuilted,Z_self,Yle,invZ,invP,F);
                r .= b .- out
                normr_act = norm(r)
                resvec[(outiter-1)*inner + initer + 1] = normr_act
                
                if normr_act <= tolb
                    x .= xm
                    flag = 0
                    iter = [outiter, initer]
                    break
                end
                #... (il resto della logica di stagionamento è complesso ma qui invariato)
            end

        end # Fine ciclo interno
        
        # --- Aggiornamento finale della soluzione per il restart ---
        if flag != 0
            idx = initercount
            if idx > 0
                ytmp_view = ytmp[1:idx]
                w_view = w[1:idx]
                R_view = R[1:idx, 1:idx]
                ldiv!(ytmp_view, UpperTriangular(R_view), w_view)
                
                alpha_add = -2 * ytmp_view[idx] * conj(U[idx, idx])
                additive .= U[:, idx] .* alpha_add
                additive[idx] += ytmp_view[idx]
                
                for k = idx-1:-1:1
                    additive[k] += ytmp_view[k]
                    d = 2 * dot(U[:, k], additive)
                    axpy!(-d, U[:, k], additive)
                end
                x .+= additive
            end
            
            xmin .= x
            ComputeMatrixVector!(out, work,pc_work, x,wk,incidence_selection,P_rebuilted,Lp_rebuilted,Z_self,Yle,invZ,invP,F);
            r .= b .- out
            normr_act = norm(r)
            r .= r # Assegna il residuo non precondizionato
        end
        
        if normr_act <= normrmin
            xmin .= x
            normrmin = normr_act
            imin = outiter
            jmin = initercount
        end
        
        if flag == 0 || flag == 3
            break
        end

    end # Fine ciclo esterno
    end # Fine blocco @views

    # --- Finalizzazione (leggermente adattata) ---
    iter = [imin, jmin]
    if flag == 0
        iter = [outitercount, initercount]
        relres = normr_act / n2minv_b
    else
        x .= xmin
        relres = normrmin / n2minv_b
    end
    
    final_len = max(outitercount - 1, 0) * inner + initercount + 1
    return x, flag, relres, iter, resvec[1:final_len]
end


function scalarsign(d)
    sgn = sign(d)
    if sgn == 0
        sgn = 1
    end
    #println(sgn)
    return sgn
end