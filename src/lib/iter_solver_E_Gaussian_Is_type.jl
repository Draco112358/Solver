using MKL
using SparseArrays, IterativeSolvers, LinearAlgebra, LinearMaps, FLoops, DelimitedFiles
using Base.Threads, JLD2
include("compute_Matrix_vector_new2.jl")
include("gmres_custom_new2.jl")
include("build_Yle_new.jl")
include("compute_Ec_Gauss.jl")
include("compute_Ar_Gauss3.jl")
include("compute_lambda_numeric.jl")
include("compute_E_field_Gauss.jl")
include("compute_H_field_Gauss.jl")
include("save_matrix_for_matlab.jl")


const μ0  = 4π * 1e-7
const ε0  = 8.854_187_816_997_944e-12
const TWO_PI = 2π
const K0 = sqrt(ε0 * μ0)        # √(μ0 ε0) è costante


function iter_solver_E_Gaussian_Is_type(
	freq,
	escalings,
	incidence_selection,
	P_data,
	Lp_data,
	ports,
	lumped_elements,
	GMRES_settings,
	volumi,
	superfici,
	use_Zs_in,
	QS_Rcc_FW,
	ports_scatter_value,
	Vs,
	Is,
	centri_oss,
	centri_oss_3D,
	id,
	chan,
	commentsEnabled
)
	@time begin
		num_oss = size(centri_oss, 1)
		num_oss_3D = size(centri_oss_3D, 1)
		ordine_int = 4
		unitario = ones(num_oss_3D)
		zeroo = zeros(num_oss_3D)

		freq .= freq .* escalings[:freq]
		# GMRES settings ----------------------------
		Inner_Iter::Int64 = GMRES_settings["Inner_Iter"]
		#Outer_Iter = GMRES_settings.Outer_Iter
		# -------------------------------------------
		mx = incidence_selection[:mx]
		my = incidence_selection[:my]
		mz = incidence_selection[:mz]

		indx = 1:mx
		indy = mx+1:mx+my
		indz = mx+my+1:mx+my+mz

		m = mx + my + mz
		n::Int64 = size(incidence_selection[:A], 2)
		ns::Int64 = size(incidence_selection[:Gamma], 2)
		w = 2 .* pi .* freq
		nfreq = length(w)
		is = zeros(ComplexF64, n)
		out = Dict(
			"Vp" => zeros(ComplexF64, size(ports[:port_start], 1), length(freq)),
			"Ex" => zeros(ComplexF64, num_oss, length(freq)),
			"Ey" => zeros(ComplexF64, num_oss, length(freq)),
			"Ez" => zeros(ComplexF64, num_oss, length(freq)),
			"Ex_3D" => zeros(ComplexF64, num_oss_3D, length(freq)),
			"Ey_3D" => zeros(ComplexF64, num_oss_3D, length(freq)),
			"Ez_3D" => zeros(ComplexF64, num_oss_3D, length(freq)),
			"Hx_3D" => zeros(ComplexF64, num_oss_3D, length(freq)),
			"Hy_3D" => zeros(ComplexF64, num_oss_3D, length(freq)),
			"Hz_3D" => zeros(ComplexF64, num_oss_3D, length(freq)),
			"f" => zeros(length(freq)),
		)
		Vrest = zeros(m + n + ns)

		invP::SparseMatrixCSC{Float64, Int64} = spdiagm(1.0 ./ diag(P_data[:P]))
		R_chiusura = ports_scatter_value
		keeped_diag = 0
		invCd = zeros(ComplexF64, m)
		not_switched = true
		# resProd = Array{ComplexF64}(undef, 2 * m)
		tn = zeros(ComplexF64, m + ns + n)
		out["f"] = freq / escalings[:freq]

		if freq[1] == 0
			freq[1] = 1e-8 / escalings[:freq]
		end
		w = 2 * pi * freq
		mu0 = 4 * pi * 1e-7
		eps0 = 8.854187816997944e-12
	end
	
	diag_Lp = zeros(Float64, size(Lp_data[:Lp_x],1)+size(Lp_data[:Lp_y],1)+size(Lp_data[:Lp_z],1))
	m_gmres   = size(incidence_selection[:A],1)
	out_gmres = similar(Vrest, ComplexF64, m_gmres+ns+size(incidence_selection[:Gamma],1))  # vettore risultato
	work = (Y1 = similar(Vrest, ComplexF64, m_gmres),
        Y2 = similar(Vrest, ComplexF64, ns),
        Y3 = similar(Vrest, ComplexF64, size(incidence_selection[:Gamma],1))
	)
	pc_work = PCWork(ComplexF64, m_gmres, ns, size(incidence_selection[:Gamma],1))
	incidence_selection[:A_t]= transpose(incidence_selection[:A])
	incidence_selection[:Gamma_t] = transpose(incidence_selection[:Gamma])	
	for k ∈ 1:nfreq
		@time begin
			β, w, keeped_diag, invP = handle_scaling_and_rebuilding!(
        k, freq, escalings, not_switched,
        volumi, P_data, Lp_data,
        Vrest, m, ns, QS_Rcc_FW,
        keeped_diag, diag_Lp, invP)
		end
		# Build Yle
		Yle = @time build_Yle_new(
			lumped_elements,
			[],
			ports,
			escalings,
			n,
			w[k] / escalings[:freq],
			R_chiusura,
			lumped_elements[:type],
			lumped_elements[:R],
			lumped_elements[:L],
			lumped_elements[:C],
		)
		@time begin
			# Compute inverse Cd for dielectric indices
			length_cd = !isempty(volumi[:indici_dielettrici]) ? length(volumi[:Cd]) : length(volumi[:R])
			invCd = zeros(ComplexF64, length_cd)
			if !isempty(volumi[:indici_dielettrici])
				invCd[Int64.(volumi[:indici_dielettrici])] .= 1 ./ (1im * w[k] * volumi[:Cd][Int64.(volumi[:indici_dielettrici])])
			end

			# Check use_Zs_in and compute Z_self
			if use_Zs_in == 1
				# Compute Zs
				Zs = real.(sqrt(1im * w[k] / escalings[:freq]) .* volumi[:Zs_part])
				indR = findall(x -> volumi[:R][x] > Zs[x], 1:length(volumi[:R]))
				indZs = setdiff(1:length(volumi[:R]), indR)
				Z_self = zeros(ComplexF64, length(volumi[:R]))
				Z_self[indR] .= volumi[:R][indR]
				Z_self[indZs] .= Zs[indZs]
				Z_self .+= invCd
			else
				Z_self = volumi[:R] .+ invCd
			end

			invZ = sparse(1:m, 1:m, 1 ./ (Z_self + 1im * w[k] * diag_Lp), m, m)
			# --------------------- preconditioner ------------------------
			SS::SparseArrays.SparseMatrixCSC{ComplexF64, Int64} = Yle + (transpose(incidence_selection[:A]) * (invZ * incidence_selection[:A])) + 1im * w[k] * (incidence_selection[:Gamma] * invP) * transpose(incidence_selection[:Gamma])
			F::SparseArrays.UMFPACK.UmfpackLU{ComplexF64, Int64} = lu(SS)
			# --------------------------------------------------------------
			for c1 in 1:size(ports[:port_nodes], 1)
				n1::Int64 = convert(Int64, ports[:port_nodes][c1, 1])
				n2::Int64 = convert(Int64, ports[:port_nodes][c1, 2])
				is[n1] = Is[c1, k] * escalings[:Is]
				is[n2] = -1.0 * Is[c1, k] * escalings[:Is]
			end

			tn = precond_3_3_vector_new(F, invZ, invP, incidence_selection[:A], incidence_selection[:Gamma], ns, Vs[:, k], is)
		end

		V, flag, relres, iter, resvec = @time gmres_custom_new2!(out_gmres, work, pc_work, tn, false, GMRES_settings["tol"][k], Inner_Iter, ComplexF64.(Vrest), w[k], incidence_selection, P_data, Lp_data, Z_self, Yle, invZ, invP, F, id, chan, 1)
		if flag == 99
			return nothing
		end

		tot_iter_number = (iter[1] - 1) * Inner_Iter + iter[2] + 1
		if commentsEnabled
			if (flag == 0)
				println("Flag $flag - Iteration = $k - Convergence reached, number of iterations:$tot_iter_number")
			end

			if (flag == 1)
				println("Flag $flag - Iteration = $k - Convergence not reached, number of iterations:$Inner_Iter")
			end
		end
		Vrest = V
		for c1 in 1:size(ports[:port_nodes], 1)
			n1::Int64 = convert(Int64, ports[:port_nodes][c1, 1])
			n2::Int64 = convert(Int64, ports[:port_nodes][c1, 2])
			out["Vp"][c1, k] = V[m+ns+n1] - V[m+ns+n2]
		end
		I = V[1:m] ./ escalings[:Is]
		J = I ./ volumi[:S]

		#send_rabbitmq_feedback(Dict("freqNumber" => k, "id" => id), "solver_feedback")

		sigma = (V[m+1:m+ns] ./ superfici["S"]) / escalings[:Cd]
		beta = 2 * pi * freq[k] / escalings[:freq] * sqrt(ε0 * μ0)
		hc = @time compute_Ec_Gauss(Float64.(superfici["estremi_celle"]), map(v -> Float64.(v), superfici["normale"]), centri_oss, ordine_int, complex(beta, 0.0), id, chan)
		if isnothing(hc)
			return nothing
		end
		#send_rabbitmq_feedback(Dict("electric_fields_results_step" => 1, "electric_fields_results_name" => "hc_3D", "id" => id), "solver_feedback")
		hc_3D = @time compute_Ec_Gauss(Float64.(superfici["estremi_celle"]), map(v -> Float64.(v), superfici["normale"]), centri_oss_3D, ordine_int, complex(beta, 0.0), id, chan)
		if isnothing(hc_3D)
			return nothing
		end
		#send_rabbitmq_feedback(Dict("electric_fields_results_step" => 2, "electric_fields_results_name" => "ha", "id" => id), "solver_feedback")
		ha = @time compute_Ar_Gauss(transpose(Float64.(volumi[:coordinate].parent)), centri_oss, ordine_int, complex(beta, 0.0), id, chan)
		if isnothing(ha)
			return nothing
		end
		#saveComplexMatrix("ha_Opt3.mat", ha, varname="haOpt3")
		#send_rabbitmq_feedback(Dict("electric_fields_results_step" => 3, "electric_fields_results_name" => "ha_3D", "id" => id), "solver_feedback")
		ha_3D = @time compute_Ar_Gauss(transpose(Float64.(volumi[:coordinate].parent)), centri_oss_3D, ordine_int, complex(beta, 0.0), id, chan)
		if isnothing(ha_3D)
			return nothing
		end
		#send_rabbitmq_feedback(Dict("electric_fields_results_step" => 4, "electric_fields_results_name" => "Lambda_x", "id" => id), "solver_feedback")
		Lambda_x = @time compute_lambda_numeric(centri_oss_3D, volumi, incidence_selection, [unitario zeroo zeroo], ordine_int, complex(beta, 0.0), id, chan)
		if isnothing(Lambda_x)
			return nothing
		end
		#send_rabbitmq_feedback(Dict("electric_fields_results_step" => 5, "electric_fields_results_name" => "Lambda_y", "id" => id), "solver_feedback")
		Lambda_y = @time compute_lambda_numeric(centri_oss_3D, volumi, incidence_selection, [zeroo unitario zeroo], ordine_int, complex(beta, 0.0), id, chan)
		if isnothing(Lambda_y)
			return nothing
		end
		#send_rabbitmq_feedback(Dict("electric_fields_results_step" => 6, "electric_fields_results_name" => "Lambda_z", "id" => id), "solver_feedback")
		Lambda_z = @time compute_lambda_numeric(centri_oss_3D, volumi, incidence_selection, [zeroo zeroo unitario], ordine_int, complex(beta, 0.0), id, chan)
		if isnothing(Lambda_z)
			return nothing
		end
		#send_rabbitmq_feedback(Dict("electric_fields_results_step" => 7, "electric_fields_results_name" => "Loading Results", "id" => id), "solver_feedback")

		out["Ex"][:, k], out["Ey"][:, k], out["Ez"][:, k] = compute_E_field_Gauss(indx, indy, indz, centri_oss, hc, ha, J, sigma, freq[k] / escalings[:freq])
		out["Ex_3D"][:, k], out["Ey_3D"][:, k], out["Ez_3D"][:, k] = compute_E_field_Gauss(indx, indy, indz, centri_oss_3D, hc_3D, ha_3D, J, sigma, freq[k] / escalings[:freq])
		out["Hx_3D"][:, k], out["Hy_3D"][:, k], out["Hz_3D"][:, k] = compute_H_field_Gauss(Lambda_x, Lambda_y, Lambda_z, I)
	end
	return out
end

function precond_3_3_Kt!(F, invZ, invP, A, Gamma, n1, n2, X3, Y, resProd)
	n3 = length(X3)
	i1 = 1:n1
	i2 = n1+1:n1+n2
	i3 = n1+n2+1:n1+n2+n3

	M5 = F \ X3

	A_view = @view resProd[1:size(A, 1)]
	invZ_view = @view resProd[size(resProd, 1)-size(invZ, 1)+1:end]
	mul!(A_view, A, M5)
	mul!(invZ_view, invZ, A_view)
	Y[i1] .= lmul!(-1.0, invZ_view)
	Gamma_view = @view resProd[size(resProd, 1)-size(Gamma, 2)+1:end]
	mul!(Gamma_view, transpose(Gamma), M5)
	invP_view = @view resProd[1:size(invP, 1)]
	mul!(invP_view, invP, Gamma_view)
	Y[i2] .= invP_view
	Y[i3] .= M5

	return Y
end

function s2z(S, Zo)
	num_ports = size(S)[1]
	nfreq = size(S)[3]
	Z = zeros(ComplexF64, num_ports, num_ports, nfreq)
	Id = Matrix{Int64}(I, num_ports, num_ports)
	for cont in 1:nfreq
		Z[:, :, cont] = Zo * ((Id - 1.0 * S[:, :, cont]) \ (Id + S[:, :, cont]))
	end
	return Z
end

function s2y(S, Zo)
	num_ports = size(S)[1]
	nfreq = size(S)[3]
	Y = zeros(ComplexF64, num_ports, num_ports, nfreq)
	Id = Matrix{Int64}(I, num_ports, num_ports)
	for cont in 1:nfreq
		Y[:, :, cont] = Zo * ((Id + S[:, :, cont]) \ (Id - 1.0 * S[:, :, cont]))
	end
	return Y
end

function precond_3_3_vector_new(lu, invZ, invP, A, Gamma, n2, X1, X3)
	n1 = length(X1)
	n3 = length(X3)

	i1 = 1:n1
	i2 = n1+1:n1+n2
	i3 = n1+n2+1:n1+n2+n3

	Y = zeros(ComplexF64, n1 + n2 + n3)

	M1 = invZ*X1
	M2 = lu \ (transpose(A)*M1)
	M5 = lu \ X3

	@views Y[i1] .= invZ*X1 .- invZ*(A*M2) .- invZ*(A*M5)
	@views Y[i2] .= invP*(transpose(Gamma)*M2) .+ invP*(transpose(Gamma)*M5)
	@views Y[i3] .= M2 .+ M5

	return Y
end

function prod_real_complex(A, x)
	# A is a N x N real matrix and x is a complex vector
	N = size(A, 1)
	y = zeros(ComplexF64, N)
	y .= A * real(x) .+ 1im * (A * imag(x))
	return y
end

# function prod_complex_real(A, x)
# 	# A is a N x N complex matrix and x is a real vector
# 	N = size(A, 1)
# 	y = zeros(ComplexF64, N)
# 	y .= real(A) * x .+ 1im * (imag(A) * x)
# 	return y
# end

function prod_complex_real(A::AbstractMatrix{<:Complex}, 
                                x::AbstractVector{<:Real})
    n = size(A,1)
    @assert size(A,2) == n == length(x)
    # Vettori d’appoggio reali (non inizializzati)
    re  = Vector{Float64}(undef, n)
    im  = Vector{Float64}(undef, n)
    # dgemv! → nessuna nuova allocazione
    mul!(re, real(A), x)   # re  = real(A)*x
    mul!(im, imag(A), x)   # im  = imag(A)*x
    y = Vector{ComplexF64}(undef, n)
    @inbounds @simd for i in eachindex(y)
        y[i] = ComplexF64(re[i], im[i])
    end
    return y
end

function convert_transpose_to_float64(t_matrix::Transpose{Real, Matrix{Real}})
    return transpose(Float64.(t_matrix.parent))
end

function handle_scaling_and_rebuilding(
    k::Int,
    freq::Vector{Float64},
    escalings::Dict{Symbol, Real},
    not_switched::Bool,
    volumi::Dict{Symbol, AbstractArray},
    P_data::Dict{Symbol, Union{Matrix{ComplexF64},Matrix{Float64}}},
    Lp_data::Dict{Symbol, Union{Matrix{ComplexF64},Matrix{Float64}}},
    Vrest::Vector{Float64}, # Assuming ComplexF64 based on usage
    m::Int64,
    ns::Int64,
    QS_Rcc_FW::Int64,
    keeped_diag::Int64
)
    # Prima parte: Gestione dello scaling
    if freq[k] / escalings[:freq] >= 1e8 && not_switched
		keeped_diag = 0
		not_switched = false

		freq = freq ./ escalings[:freq]
		w = 2 .* pi .* freq
		volumi[:R] = volumi[:R] ./ escalings[:R]
		volumi[:Zs_part] = volumi[:Zs_part] ./ escalings[:R]
		if !isempty(volumi[:indici_dielettrici])
			volumi[:Cd] = volumi[:Cd] ./ escalings[:Cd]
		end
		P_data[:P] = P_data[:P] ./ escalings[:P]
		invP = spdiagm(1.0 ./ diag(P_data[:P]))

		Lp_data[:Lp_x] = Lp_data[:Lp_x] ./ escalings[:Lp]
		Lp_data[:Lp_y] = Lp_data[:Lp_y] ./ escalings[:Lp]
		Lp_data[:Lp_z] = Lp_data[:Lp_z] ./ escalings[:Lp]

		Vrest[1:m] = Vrest[1:m] ./ escalings[:Is]
		Vrest[m+1:m+ns] = Vrest[m+1:m+ns] ./ escalings[:Cd]

		escalings[:Lp] = 1.0
		escalings[:R] = 1.0
		escalings[:Cd] = 1.0
		escalings[:P] = 1.0
		escalings[:Is] = 1.0
		escalings[:freq] = 1.0
		escalings[:Yle] = 1.0
		escalings[:time] = 1.0
	end
	if QS_Rcc_FW == 2
		mu0 = 4 * π * 1e-7
		eps0 = 8.854187816997944e-12
		beta = 2 * π * freq[k] / escalings[:freq] * sqrt(eps0 * mu0)
		# Rebuilding Lp_data with phase correction
		Lp_rebuilted = Dict(
			:Lp_x => Lp_data[:Lp_x] .* exp.(-1im * beta * Lp_data[:Rx]),
			:Lp_y => Lp_data[:Lp_y] .* exp.(-1im * beta * Lp_data[:Ry]),
			:Lp_z => Lp_data[:Lp_z] .* exp.(-1im * beta * Lp_data[:Rz]),
		)
		@show size(Lp_rebuilted[:Lp_x])
		@show size(Lp_rebuilted[:Lp_y])
		@show size(Lp_rebuilted[:Lp_z])
		diag_Lp = vcat(
			diag(real.(Lp_rebuilted[:Lp_x])),
			diag(real.(Lp_rebuilted[:Lp_y])),
			diag(real.(Lp_rebuilted[:Lp_z])),
		)
		# Rebuilding P_data with phase correction
		P_rebuilted = Dict(
			:P => P_data[:P] .* exp.(-1im * beta * P_data[:R_cc]),
		)
	else
		P_rebuilted = P_data
		Lp_rebuilted = Lp_data
		if keeped_diag == 0
			diag_Lp = vcat(
				diag(real.(Lp_data[:Lp_x])),
				diag(real.(Lp_data[:Lp_y])),
				diag(real.(Lp_data[:Lp_z])),
			)
			keeped_diag = 1
		end
	end
    
    return freq, escalings, not_switched, volumi, P_data, Lp_data, Vrest, keeped_diag, P_rebuilted, Lp_rebuilted, diag_Lp, w, invP, eps0, mu0
end


function promote_to_complex!(d::Dict, keys)
    for k in keys
        d[k] = ComplexF64.(d[k])   # crea la nuova matrice e la rimpiazza
    end
end

function handle_scaling_and_rebuilding!(
        k::Int, freq::Vector{Float64}, esc::Dict{Symbol,<:Real},
        not_switched::Bool,
        volumi::Dict{Symbol,AbstractArray},
        P_data::Dict{Symbol,<:Union{Matrix{Float64},Matrix{ComplexF64}}},
        Lp_data::Dict{Symbol,<:Union{Matrix{Float64},Matrix{ComplexF64}}},
        Vrest::Vector{Float64},
        m::Int, ns::Int,
        QS_Rcc_FW::Int,
        keeped_diag::Int,
        diag_Lp::Vector{Float64},
        invP::SparseMatrixCSC{Float64,Int64}
    )
    # ── 1. Mega-scaling una sola volta ────────────────────────────
    if freq[k] / esc[:freq] ≥ 1e8 && not_switched
        not_switched = false
        keeped_diag  = 0
        invfreq = 1 / esc[:freq]
        @threads for i in eachindex(freq)
            freq[i] *= invfreq
        end
        # ================ blocco di scaling parallelo ================
        invR, invCd = 1/esc[:R], 1/esc[:Cd]
        invP_fac, invLp = 1/esc[:P], 1/esc[:Lp]
        invIs = 1/esc[:Is]
        @threads for i in eachindex(volumi[:R])
            volumi[:R][i]       *= invR
            volumi[:Zs_part][i] *= invR
        end
        if !isempty(volumi[:indici_dielettrici])
            @threads for i in eachindex(volumi[:Cd])
                volumi[:Cd][i] *= invCd
            end
        end
        @threads for i in eachindex(P_data[:P])
            P_data[:P][i] *= invP_fac
        end
        @threads for i in eachindex(Lp_data[:Lp_x])
            Lp_data[:Lp_x][i] *= invLp
        end
        @threads for i in eachindex(Lp_data[:Lp_y])
            Lp_data[:Lp_y][i] *= invLp
        end
        @threads for i in eachindex(Lp_data[:Lp_z])
            Lp_data[:Lp_z][i] *= invLp
        end
        @threads for i in 1:m
            Vrest[i] *= invIs
        end
        @threads for i in 1:ns
            Vrest[m+i] *= invCd
        end
        # ============================================================
        # invdiag parallelo
        diagP = diag(P_data[:P])                 # view
        invdiag = similar(diagP)
        @threads for i in eachindex(diagP)
            invdiag[i] = 1.0 / diagP[i]
        end
        invP = spdiagm(0 => invdiag)             # una sola allocazione
        esc[:Lp] = esc[:R] = esc[:Cd] = esc[:P] =
        esc[:Is] = esc[:freq] = esc[:Yle] = esc[:time] = 1.0
    end
    # ── 2. Costanti di passo ───────────────────────────────────────
    w = TWO_PI * freq[k]
    β = w * K0
    # ── 3. Correzione di fase (se richiesta) ───────────────────────
    if QS_Rcc_FW == 2
        promote_to_complex!(Lp_data, [:Lp_x, :Lp_y, :Lp_z])
        promote_to_complex!(P_data,  [:P])
        Lpx = Lp_data[:Lp_x];  Lpy = Lp_data[:Lp_y];  Lpz = Lp_data[:Lp_z]
        Rx  = Lp_data[:Rx];    Ry  = Lp_data[:Ry];    Rz  = Lp_data[:Rz]
        Pm  = P_data[:P];      Rcc = P_data[:R_cc]
        # tre loop indipendenti → tre thread-loops
        @threads for i in eachindex(Lpx); Lpx[i] *= cis(-β * Rx[i]); end
        @threads for i in eachindex(Lpy); Lpy[i] *= cis(-β * Ry[i]); end
        @threads for i in eachindex(Lpz); Lpz[i] *= cis(-β * Rz[i]); end
        @threads for i in eachindex(Pm);  Pm[i]  *= cis(-β * Rcc[i]); end
    end
    # ── 4. Salvataggio diagonale una sola volta ────────────────────
    if keeped_diag == 0
        n  = size(Lp_data[:Lp_x],1)
        n2 = size(Lp_data[:Lp_y],1)
        n3 = size(Lp_data[:Lp_z],1)
        @threads for i in 1:n
            diag_Lp[i] = real(Lp_data[:Lp_x][i,i])
        end
        @threads for i in 1:n2
            diag_Lp[n+i] = real(Lp_data[:Lp_y][i,i])
        end
        @threads for i in 1:n3
            diag_Lp[n+n2+i] = real(Lp_data[:Lp_z][i,i])
        end
        keeped_diag = 1
    end
    return β, w, keeped_diag, invP
end