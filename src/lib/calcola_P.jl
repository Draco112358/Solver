include("Song_P_improved_Ivana_strategy2.jl")

function calcola_P(superfici, escalings, QS_Rcc_FW, id)

	eps0 = 8.854187816997944e-12

	epsilon1 = 5e-3
	epsilon2 = 1e-3
	epsilon3 = 1e-3
	epsilon4 = 3e-1

	use_suppression = 1
	superfici["estremi_celle"] = hcat(superfici["estremi_celle"]...)
	superfici["centri"] = hcat(superfici["centri"]...)
	nsup = size(superfici["estremi_celle"], 1)
	# Inside calcola_P, after nsup is determined
	# Allocate arrays to store the calculated properties for each surface
	x_vs = Vector{Vector{Float64}}(undef, nsup)
	y_vs = Vector{Vector{Float64}}(undef, nsup)
	z_vs = Vector{Vector{Float64}}(undef, nsup)

	xc_s = Vector{Float64}(undef, nsup)
	yc_s = Vector{Float64}(undef, nsup)
	zc_s = Vector{Float64}(undef, nsup)

	a_s = Vector{Float64}(undef, nsup)
	b_s = Vector{Float64}(undef, nsup)
	c_s = Vector{Float64}(undef, nsup)

	sup_s = Vector{Float64}(undef, nsup)
	sup_xz_planes = Vector{Int}(undef, nsup)
	sup_yz_planes = Vector{Int}(undef, nsup)

	round_precision = 14 # Define this once

	# Loop once to pre-calculate all these properties for all surfaces
	Threads.@threads for i in 1:nsup
		# Use @view to avoid copying the surface array
		current_surface = @view superfici["estremi_celle"][i, :]

		x_vs[i] = unique(round.(current_surface[[1, 10]], digits=round_precision))
		y_vs[i] = unique(round.(current_surface[[2, 11]], digits=round_precision))
		z_vs[i] = unique(round.(current_surface[[3, 12]], digits=round_precision))

		xc_s[i] = 0.5 * (x_vs[i][end] + x_vs[i][1])
		yc_s[i] = 0.5 * (y_vs[i][end] + y_vs[i][1])
		zc_s[i] = 0.5 * (z_vs[i][end] + z_vs[i][1])

		a_val = abs(x_vs[i][end] - x_vs[i][1])
		b_val = abs(y_vs[i][end] - y_vs[i][1])
		c_val = abs(z_vs[i][end] - z_vs[i][1])

		sup1_yz_plane = 0
		sup1_xz_plane = 0

		if (a_val <= b_val && a_val <= c_val)
			sup1_yz_plane = 1
			a_val = 1.0 # Ensure float type consistency
		elseif (b_val <= a_val && b_val <= c_val)
			sup1_xz_plane = 1
			b_val = 1.0
		else
			c_val = 1.0
		end

		a_s[i] = a_val
		b_s[i] = b_val
		c_s[i] = c_val
		sup_s[i] = a_val * b_val * c_val
		sup_xz_planes[i] = sup1_xz_plane
		sup_yz_planes[i] = sup1_yz_plane
	end
	R_cc = []
	println("Calcolo P initialization")
	if QS_Rcc_FW >= 2
		R_cc = zeros(nsup, nsup)
		block_size1 = 200  # ad esempio, 10 iterazioni per blocco

		for m_block in 1:block_size1:nsup
			m_end = min(m_block + block_size1 - 1, nsup)
			Threads.@threads for m in m_block:m_end
				for n in m:nsup
					dist = norm(view(superfici["centri"], m, :) .- view(superfici["centri"], n, :))
					R_cc[m, n] = dist
					R_cc[n, m] = dist
					if n % 20 == 0
						yield()  # Permette alla task dell'heartbeat di essere schedulata
					end
				end
			end
			# Dopo aver processato un blocco, cedi il controllo per permettere la gestione dei heartbeat e altre operazioni
			sleep(0)
			if is_stop_requested(id)
				println("Simulazione $(id) interrotta per richiesta stop.")
				return nothing # O un altro valore che indica interruzione
			end
			println("block: ", round(m_end / block_size1), " / ", round(nsup / block_size1))
		end
	end
	P = zeros(nsup, nsup)
	println("Calcolo P Song_P_improved_Ivana_strategy")

	# Scegli una dimensione di blocco adatta (da regolare in base alle tue esigenze)
	block_size2 = 200  # ad esempio, processa 1000 iterazioni per blocco

	for m_block in 1:block_size2:nsup
		m_end = min(m_block + block_size2 - 1, nsup)
		Threads.@threads for m in m_block:m_end
			for n in m:nsup
				# integ, _ = Song_P_improved_Ivana_strategy(
				# 	superfici["estremi_celle"][m, :],
				# 	superfici["estremi_celle"][n, :],
				# 	epsilon1, epsilon2, epsilon3, epsilon4,
				# 	use_suppression)
				integ, _ = Song_P_improved_Ivana_strategy(
					x_vs[m], y_vs[m], z_vs[m], xc_s[m], yc_s[m], zc_s[m], a_s[m], b_s[m], c_s[m], sup_s[m], sup_yz_planes[m], sup_xz_planes[m],
					x_vs[n], y_vs[n], z_vs[n], xc_s[n], yc_s[n], zc_s[n], a_s[n], b_s[n], c_s[n], sup_s[n], sup_yz_planes[n], sup_xz_planes[n],
					epsilon1, epsilon2, epsilon3, epsilon4, use_suppression)
				P[m, n] = 1 / (4 * π * eps0 * superfici["S"][m] * superfici["S"][n]) *
						  integ * escalings[:P]
				P[n, m] = P[m, n]
				if n % 20 == 0
					yield()  # Permette alla task dell'heartbeat di essere schedulata
				end
			end
		end
		# Al termine di ogni blocco, cedi il controllo per consentire l'invio dei heartbeat
		sleep(0)
		if is_stop_requested(id)
			println("Simulazione $(id) interrotta per richiesta stop.")
			return nothing # O un altro valore che indica interruzione
		end
		println("block: ", round(m_end / block_size2), " / ", round(nsup / block_size2))
	end

	return Dict(
		:P => P,
		:R_cc => R_cc,
	)
end
