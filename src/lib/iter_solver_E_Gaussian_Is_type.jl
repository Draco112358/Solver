using MKL
using SparseArrays, IterativeSolvers, LinearAlgebra, LinearMaps, FLoops
using Base.Threads
include("compute_Matrix_vector_new2.jl")
include("gmres_custom_new2.jl")
include("build_Yle_S_new2.jl")
include("compute_Ec_Gauss.jl")
include("compute_Ar_Gauss.jl")
include("compute_lambda_numeric.jl")
include("compute_E_field_Gauss.jl")
include("compute_H_field_Gauss.jl")

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
	commentsEnabled,
)
	num_oss = size(centri_oss, 1)
	num_oss_3D = size(centri_oss_3D, 1)
	ordine_int = 3
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
	resProd = Array{ComplexF64}(undef, 2 * m)
	tn = zeros(ComplexF64, m + ns + n)
	out["f"] = freq / escalings[:freq]

	if freq[1] == 0
		freq[1] = 1e-8 / escalings[:freq]
	end
	w = 2 * pi * freq
	mu0 = 4 * pi * 1e-7
	eps0 = 8.854187816997944e-12

	for k ∈ 1:nfreq
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
		# Build Yle
		Yle = build_Yle_S_new2(
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
		V, flag, relres, iter, resvec = gmres_custom_new2(tn, false, GMRES_settings["tol"][k], Inner_Iter, Vrest, w[k], incidence_selection, P_rebuilted, Lp_rebuilted, Z_self, Yle, invZ, invP, F, resProd, id, chan, 1)
		if flag == 99
			return false
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

		if !isnothing(chan)
            publish_data(Dict("freqNumber" => k, "id" => id), "solver_feedback", chan)
        end

		sigma = (V[m+1:m+ns] ./ superfici["S"]) / escalings[:Cd]
		beta = 2 * pi * freq[k] / escalings[:freq] * sqrt(eps0 * mu0)

		hc = compute_Ec_Gauss(superfici["estremi_celle"], superfici["normale"], centri_oss, ordine_int, beta, id, chan)
		if hc == false
			return false
		end
		println("hc")
		if !isnothing(chan)
            publish_data(Dict("electric_fields_results_step" => 1, "electric_fields_results_name" => "hc_3D", "id" => id), "solver_feedback", chan)
        end
		hc_3D = compute_Ec_Gauss(superfici["estremi_celle"], superfici["normale"], centri_oss_3D, ordine_int, beta, id, chan)
		if hc_3D == false
			return false
		end
		println("hc_3D")
		if !isnothing(chan)
            publish_data(Dict("electric_fields_results_step" => 2, "electric_fields_results_name" => "ha", "id" => id), "solver_feedback", chan)
        end
		ha = compute_Ar_Gauss(volumi[:coordinate], centri_oss, ordine_int, beta, id, chan)
		if ha == false
			return false
		end
		println("ha")
		if !isnothing(chan)
            publish_data(Dict("electric_fields_results_step" => 3, "electric_fields_results_name" => "ha_3D", "id" => id), "solver_feedback", chan)
        end
		ha_3D = compute_Ar_Gauss(volumi[:coordinate], centri_oss_3D, ordine_int, beta, id, chan)
		if ha_3D == false
			return false
		end
		println("ha_3D")
		if !isnothing(chan)
            publish_data(Dict("electric_fields_results_step" => 4, "electric_fields_results_name" => "Lambda_x", "id" => id), "solver_feedback", chan)
        end
		Lambda_x = compute_lambda_numeric(centri_oss_3D, volumi, incidence_selection, [unitario zeroo zeroo], ordine_int, beta, id, chan)
		if Lambda_x == false
			return false
		end
		println("Lambda_x")
		if !isnothing(chan)
            publish_data(Dict("electric_fields_results_step" => 5, "electric_fields_results_name" => "Lambda_y", "id" => id), "solver_feedback", chan)
        end
		Lambda_y = compute_lambda_numeric(centri_oss_3D, volumi, incidence_selection, [zeroo unitario zeroo], ordine_int, beta, id, chan)
		if Lambda_y == false
			return false
		end
		println("Lambda_y")
		if !isnothing(chan)
            publish_data(Dict("electric_fields_results_step" => 6, "electric_fields_results_name" => "Lambda_z", "id" => id), "solver_feedback", chan)
        end
		Lambda_z = compute_lambda_numeric(centri_oss_3D, volumi, incidence_selection, [zeroo zeroo unitario], ordine_int, beta, id, chan)
		if Lambda_z == false
			return false
		end
		println("Lambda_z")
		if !isnothing(chan)
            publish_data(Dict("electric_fields_results_step" => 7, "electric_fields_results_name" => "Loading Results", "id" => id), "solver_feedback", chan)
        end

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

	M1 = prod_real_complex(invZ, X1)
	M2 = lu \ prod_real_complex(transpose(A), M1)
	M5 = lu \ X3

	@views Y[i1] .= prod_real_complex(invZ, X1) .- prod_real_complex(invZ, prod_real_complex(A, M2)) .- prod_real_complex(invZ, prod_real_complex(A, M5))
	@views Y[i2] .= prod_real_complex(invP, prod_real_complex(transpose(Gamma), M2)) .+ prod_real_complex(invP, prod_real_complex(transpose(Gamma), M5))
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

function prod_complex_real(A, x)
	# A is a N x N complex matrix and x is a real vector
	N = size(A, 1)
	y = zeros(ComplexF64, N)
	y .= real(A) * x .+ 1im * (imag(A) * x)
	return y
end
