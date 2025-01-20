using MKL
using SparseArrays, IterativeSolvers, LinearAlgebra, LinearMaps, FLoops
include("compute_Matrix_vector_new.jl")
include("build_Yle_S_new.jl")


function iter_solver_QS_S_type(freq,escalings,incidence_selection,P_data,Lp_data,ports,lumped_elements,GMRES_settings,volumi,use_Zs_in,QS_Rcc_FW, ports_scatter_value, id, chan, commentsEnabled)
    freq = freq .* escalings["freq"]
    # GMRES settings ----------------------------

    Inner_Iter::Int64 = GMRES_settings["Inner_Iter"]
    #Outer_Iter = GMRES_settings.Outer_Iter
    # -------------------------------------------
    mx=incidence_selection["mx"]
    my=incidence_selection["my"]
    mz=incidence_selection["mz"]
    m=mx+my+mz
    incidence_selection["A"] = hcat(incidence_selection["A"]...)
    incidence_selection["Gamma"] = hcat(incidence_selection["Gamma"]...)
    n::Int64 = size(incidence_selection["A"], 2)
    ns::Int64 = size(incidence_selection["Gamma"], 2)
    w = 2 .* pi .* freq
    nfreq = length(w)
    is = zeros(Float64, n, 1)
    S = zeros(ComplexF64, size(ports[:port_nodes], 1), size(ports[:port_nodes], 1), length(freq))
    Vrest = zeros(ComplexF64, m + n + ns, size(ports[:port_nodes], 1))
    invP::SparseMatrixCSC{Float64,Int64} = sparse(1:ns, 1:ns, 1 ./ diag(P_data[:P]), ns, ns)
    R_chiusura=ports_scatter_value;
    keeped_diag=0;
    invCd=zeros(m,1);

    resProd = Array{ComplexF64}(undef, 2 * m)
    # resProd = zeros(ComplexF64, 2*m,1)
    # tn = Array{ComplexF64}(undef, m + ns + n)
    tn = zeros(ComplexF64, m + ns + n)

    for k = 1:nfreq
        if QS_Rcc_FW == 2
            mu0 = 4 * π * 1e-7
            eps0 = 8.854187816997944e-12
            beta = 2 * π * freq[k] / escalings["freq"] * sqrt(eps0 * mu0)
            # Rebuilding Lp_data with phase correction
            Lp_rebuilted = Dict(
                :Lp_x => Lp_data[:Lp_x] .* exp.(-1im * beta * Lp_data[:Rx]),
                :Lp_y => Lp_data[:Lp_y] .* exp.(-1im * beta * Lp_data[:Ry]),
                :Lp_z => Lp_data[:Lp_z] .* exp.(-1im * beta * Lp_data[:Rz])
            )
            diag_Lp = vcat(
                diag(real.(Lp_rebuilted[:Lp_x])),
                diag(real.(Lp_rebuilted[:Lp_y])),
                diag(real.(Lp_rebuilted[:Lp_z]))
            )
            # Rebuilding P_data with phase correction
            P_rebuilted = Dict(
                :P => P_data[:P] .* exp.(-1im * beta * P_data[:R_cc])
            )
        else
            P_rebuilted = P_data
            Lp_rebuilted = Lp_data
            diag_Lp = vcat(
                    diag(real.(Lp_data[:Lp_x])),
                    diag(real.(Lp_data[:Lp_y])),
                    diag(real.(Lp_data[:Lp_z]))
                )
            if keeped_diag == 0
                diag_Lp = vcat(
                    diag(real.(Lp_data[:Lp_x])),
                    diag(real.(Lp_data[:Lp_y])),
                    diag(real.(Lp_data[:Lp_z]))
                )
                keeped_diag = 1
            end
        end
        # Build Yle
        Yle = build_Yle_S_new(
            lumped_elements,
            [],
            ports,
            escalings,
            n,
            w[k] / escalings["freq"],
            R_chiusura,
            lumped_elements[:type],
            lumped_elements[:R],
            lumped_elements[:L],
            lumped_elements[:C]
        )
        # Compute inverse Cd for dielectric indices
        invCd = zeros(ComplexF64, length(volumi["Cd"]))
        if !isempty(volumi["indici_dielettrici"])
            invCd[Int64.(volumi["indici_dielettrici"])] .= 1 ./ (1im * w[k] * volumi["Cd"][Int64.(volumi["indici_dielettrici"])])
        end

        # Check use_Zs_in and compute Z_self
        if use_Zs_in == 1
            # Compute Zs
            Zs = real.(sqrt(1im * w[k]) .* volumi["Zs_part"])
            indR = findall(x -> volumi["R"][x] > Zs[x], 1:length(volumi["R"]))
            indZs = setdiff(1:length(volumi["R"]), indR)
            # Initialize Z_self
            Z_self = zeros(ComplexF64, length(volumi["R"]))
            Z_self[indR] .= volumi["R"][indR]
            Z_self[indZs] .= Zs[indZs]
            Z_self .+= invCd
        else
            Z_self = volumi["R"] .+ invCd
        end

        invZ=sparse(1:m, 1:m,1 ./ (Z_self+1im*w[k]*diag_Lp), m, m);
        # --------------------- preconditioner ------------------------
        SS::SparseArrays.SparseMatrixCSC{ComplexF64,Int64} = Yle + (transpose(incidence_selection["A"]) * (invZ * incidence_selection["A"])) + 1im * w[k] * (incidence_selection["Gamma"] * invP) * transpose(incidence_selection["Gamma"])
        F::SparseArrays.UMFPACK.UmfpackLU{ComplexF64,Int64} = lu(SS)
        # --------------------------------------------------------------
        for c1::Int64 = 1:size(ports[:port_nodes], 1)
            n1::Int64 = convert(Int64, ports[:port_nodes][c1, 1])
            n2::Int64 = convert(Int64, ports[:port_nodes][c1, 2])
            is[n1] = escalings["Is"]
            # is[n1] = 1 * escalings["Is"]
            is[n2] = -1.0 * escalings["Is"]
            tn = precond_3_3_Kt!(F, invZ, invP, incidence_selection["A"], incidence_selection["Gamma"], m, ns, vec(is), tn, resProd)
            #println("end precond 3 3 kt")
            #ComputeMatrixVector(x , wk, incidence_selection, FFTCP, FFTCLp, DZ, Yle, expansions, invZ, invP, lu, PLIVector, PVector, PLI2Vector, P2Vector, chi2Vector)
            products_law = x ->   ComputeMatrixVectorNew(x,w[k],incidence_selection,P_rebuilted,Lp_rebuilted,Z_self,Yle,invZ,invP,F,resProd);
            prodts = LinearMap{ComplexF64}(products_law, n + m + ns, n + m + ns)
            x0 = Vrest[:,c1]
            #prob = LinearProblem(prodts, vec(tn))
            if QS_Rcc_FW == 1

                V, info = IterativeSolvers.gmres!(x0, prodts, tn; reltol=GMRES_settings["tol"][k], restart=Inner_Iter, maxiter=Inner_Iter, initially_zero=false, log=true, verbose=true)
                if commentsEnabled
                    println(info)
                end
            else
                println("start gmres Rcc")
                V, info = IterativeSolvers.gmres!(x0, prodts, tn; reltol=GMRES_settings["tol"][k], restart=Inner_Iter, maxiter=Inner_Iter, initially_zero=false, log=true, verbose=true)
                if commentsEnabled
                    println(info)
                end
            end

            Vrest[:, c1] = V
            is[n1] = 0
            is[n2] = 0
            for c2::Int64 = c1:size(ports[:port_nodes], 1)
                n3::Int64 = convert(Int64, ports[:port_nodes][c2, 1])
                n4::Int64 = convert(Int64, ports[:port_nodes][c2, 2])
                if c1 == c2
                    S[c1, c2, k] = (2 * (V[m+ns+n3] - V[m+ns+n4]) - R_chiusura) / R_chiusura
                else
                    S[c1, c2, k] = (2 * (V[m+ns+n3] - V[m+ns+n4])) / R_chiusura
                end
                S[c2, c1, k] = S[c1, c2, k]
            end
        end
        if !isnothing(chan)
            publish_data(Dict("freqNumber" => k, "id" => id), "solver_feedback", chan)
        end
        if (commentsEnabled == true)
            partial_res = dump_json_data(s2z(S, ports_scatter_value), S, s2y(S, ports_scatter_value), size(ports[:port_nodes], 1), id; partial=true, freqIndex=k)
            publish_data(partial_res, "solver_results", chan)
        end
    end
    out::Dict = Dict()
    out[:S] = S
    out[:Z] = s2z(S, R_chiusura)
    out[:Y] = s2y(S, R_chiusura)
    out[:f] = freq ./ escalings["freq"]
    return out

end


function precond_3_3_Kt!(F, invZ, invP, A, Gamma, n1, n2, X3, Y, resProd)
    n3 = length(X3)
    i1 = range(1, stop=n1)
    i2 = range(n1 + 1, stop=n1 + n2)
    i3 = range(n1 + n2 + 1, stop=n1 + n2 + n3)

    # Y = zeros(ComplexF64, n1 + n2 + n3, 1)
    M5 = F \ X3

    # Yi1 = @view Y[i1]
    A_view = @view resProd[1:size(A, 1)]
    invZ_view = @view resProd[size(resProd, 1)-size(invZ, 1)+1:end]
    mul!(A_view, A, M5)
    mul!(invZ_view, invZ, A_view)
    Y[i1] .= lmul!(-1.0, invZ_view)
    # Yi2 = @view Y[i2]
    Gamma_view = @view resProd[size(resProd, 1)-size(Gamma, 2)+1:end]
    mul!(Gamma_view, transpose(Gamma), M5)
    invP_view = @view resProd[1:size(invP, 1)]
    mul!(invP_view, invP, Gamma_view)
    Y[i2] .= invP_view
    # Yi3 = @view Y[i3]
    Y[i3] .= M5

    return Y
end


function s2z(S, Zo)
    num_ports = size(S)[1]
    nfreq = size(S)[3]
    Z = zeros(ComplexF64, num_ports, num_ports, nfreq)
    Id = Matrix{Int64}(I, num_ports, num_ports)
    for cont in range(1, stop=nfreq)
        Z[:, :, cont] = Zo * ((Id - 1.0 * S[:, :, cont]) \ (Id + S[:, :, cont]))
        #Z[:, :, cont] = ((Id - 1.0 * S[:, :, cont]) \ (Id + S[:, :, cont]))
    end
    return Z
end

function s2y(S, Zo)
    num_ports = size(S)[1]
    nfreq = size(S)[3]
    Y = zeros(ComplexF64, num_ports, num_ports, nfreq)
    Id = Matrix{Int64}(I, num_ports, num_ports)
    for cont in range(1, stop=nfreq)
        Y[:, :, cont] = Zo * ((Id + S[:, :, cont]) \ (Id - 1.0 * S[:, :, cont]))
        #Y[:, :, cont] = ((Id + S[:, :, cont]) \ (Id - 1.0 * S[:, :, cont]))
    end
    return Y
end