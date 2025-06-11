include("Compute_Lp_Self2.jl")
include("Song_improved_Ivana_strategy2.jl")

function calcola_Lp2(volumi, incidence_selection, escalings, QS_Rcc_FW, id)
	epsilon1 = 5e-3
	epsilon2 = 1e-3
	epsilon3 = 1e-3
	epsilon4 = 3e-1

	use_suppression = true

	mx = incidence_selection[:mx]
	my = incidence_selection[:my]
	mz = incidence_selection[:mz]
	ntot = mx + my + mz # Numero totale di volumi

	# --- INIZIO PRE-CALCOLO DELLE PROPRIETÀ DEI VOLUMI ---
    # Questa parte viene eseguita una sola volta per tutti i volumi.
    round_precision = 14

    # Definisci gli indici una sola volta
    indici_x_min = [1, 7, 13, 19]
    indici_y_min = [2, 5, 14, 17]
    indici_z_min = [3, 6, 9, 12]
    indici_x_max = [4, 10, 16, 22]
    indici_y_max = [8, 11, 20, 23]
    indici_z_max = [15, 18, 21, 24]

    x_vs = Vector{Vector{Float64}}(undef, ntot)
    y_vs = Vector{Vector{Float64}}(undef, ntot)
    z_vs = Vector{Vector{Float64}}(undef, ntot)
    
    xc_s = Vector{Float64}(undef, ntot)
    yc_s = Vector{Float64}(undef, ntot)
    zc_s = Vector{Float64}(undef, ntot)

    a_s = Vector{Float64}(undef, ntot)
    b_s = Vector{Float64}(undef, ntot)
    c_s = Vector{Float64}(undef, ntot)
    V_s = Vector{Float64}(undef, ntot)

	# --- INIZIO PRE-CALCOLO DEI PARAMETRI PER Compute_Lp_Self2 ---
    l_s = Vector{Float64}(undef, ntot)
    W_s = Vector{Float64}(undef, ntot)
    T_s = Vector{Float64}(undef, ntot)
    # --- FINE INIZIO PRE-CALCOLO ---

    Base.Threads.@threads for i in 1:ntot
        current_volume_coords = @view volumi[:coordinate][i, :] # Usa @view per evitare copie

        x_vs[i] = unique(round.(current_volume_coords[[indici_x_min; indici_x_max]], digits=round_precision))
        y_vs[i] = unique(round.(current_volume_coords[[indici_y_min; indici_y_max]], digits=round_precision))
        z_vs[i] = unique(round.(current_volume_coords[[indici_z_min; indici_z_max]], digits=round_precision))

        xc_s[i] = 0.5 * (x_vs[i][end] + x_vs[i][1])
        yc_s[i] = 0.5 * (y_vs[i][end] + y_vs[i][1])
        zc_s[i] = 0.5 * (z_vs[i][end] + z_vs[i][1])

        a_s[i] = abs(x_vs[i][end] - x_vs[i][1])
        b_s[i] = abs(y_vs[i][end] - y_vs[i][1])
        c_s[i] = abs(z_vs[i][end] - z_vs[i][1])

        V_s[i] = a_s[i] * b_s[i] * c_s[i]
		# Pre-calcolo per Compute_Lp_Self2
        curr_dir_temp = 0 # Questo deve essere gestito correttamente in base al contesto del volume
        if i <= mx
            curr_dir_temp = 1
        elseif i <= mx + my
            curr_dir_temp = 2
        else
            curr_dir_temp = 3
        end

        if curr_dir_temp == 1
            l_s[i] = abs(mean_length_rev(current_volume_coords, 1))
            W_s[i] = abs(mean_length_rev(current_volume_coords, 3))
            T_s[i] = abs(mean_length_rev(current_volume_coords, 2))
        elseif curr_dir_temp == 2
            T_s[i] = abs(mean_length_rev(current_volume_coords, 1))
            W_s[i] = abs(mean_length_rev(current_volume_coords, 3))
            l_s[i] = abs(mean_length_rev(current_volume_coords, 2))
        else # curr_dir_temp == 3
            W_s[i] = abs(mean_length_rev(current_volume_coords, 1))
            l_s[i] = abs(mean_length_rev(current_volume_coords, 3))
            T_s[i] = abs(mean_length_rev(current_volume_coords, 2))
        end
    end
    # --- FINE PRE-CALCOLO ---

	Rx = nothing
	Ry = nothing
	Rz = nothing

	if QS_Rcc_FW >= 2
		Rx = Matrix{Float64}(undef, mx, mx)
		Ry = Matrix{Float64}(undef, my, my)
		Rz = Matrix{Float64}(undef, mz, mz)

		# Dimensione del blocco (da regolare in base alle esigenze)
		block_size = 200

		# Suddividiamo il ciclo sul parametro m in blocchi
		for m_block in 1:block_size:mx
			m_end = min(m_block + block_size - 1, mx)

			Threads.@threads for m in m_block:m_end
				for n in m:mx
					dist = norm(view(volumi[:centri], m, :) .- view(volumi[:centri], n, :))
					Rx[m, n] = dist
					Rx[n, m] = dist
					if n % 20 == 0
						yield()  # Permette alla task dell'heartbeat di essere schedulata
					end
				end
			end

			# Alla fine di ogni blocco, cede il controllo per permettere ad altre operazioni
			sleep(0)
			if is_stop_requested(id)
				println("Simulazione $(id) interrotta per richiesta stop.")
				return nothing # O un altro valore che indica interruzione
			end
			println("block Lp1 : ", round(m_end / block_size), " / ", round(mx / block_size))
		end

		for m_block in 1:block_size:my
			m_end = min(m_block + block_size - 1, my)

			Threads.@threads for m in m_block:m_end
				for n in m:my
					dist = norm(view(volumi[:centri], m + mx, :) .- view(volumi[:centri], n + mx, :))
					Ry[m, n] = dist
					Ry[n, m] = dist
					if n % 20 == 0
						yield()  # Permette alla task dell'heartbeat di essere schedulata
					end
				end
			end

			# Alla fine di ogni blocco, cede il controllo per permettere ad altre operazioni
			sleep(0)
			if is_stop_requested(id)
				println("Simulazione $(id) interrotta per richiesta stop.")
				return nothing # O un altro valore che indica interruzione
			end
			println("block Lp2 : ", round(m_end / block_size), " / ", round(my / block_size))
		end

		for m_block in 1:block_size:mz
			m_end = min(m_block + block_size - 1, mz)

			Threads.@threads for m in m_block:m_end
				for n in m:mz
					dist = norm(view(volumi[:centri], m + mx + my, :) .- view(volumi[:centri], n + mx + my, :))
					Rz[m, n] = dist
					Rz[n, m] = dist
					if n % 20 == 0
						yield()  # Permette alla task dell'heartbeat di essere schedulata
					end
				end
			end

			# Alla fine di ogni blocco, cede il controllo per permettere ad altre operazioni
			sleep(0)
			if is_stop_requested(id)
				println("Simulazione $(id) interrotta per richiesta stop.")
				return nothing # O un altro valore che indica interruzione
			end
			println("block Lp3 : ", round(m_end / block_size), " / ", round(mz / block_size))
		end
	end

	Lp_x = Matrix{Float64}(undef, mx, mx)
	Lp_y = Matrix{Float64}(undef, my, my)
	Lp_z = Matrix{Float64}(undef, mz, mz)

	# Dimensione del blocco (da regolare in base alle esigenze)
	block_size = 200

	# Suddividiamo il ciclo sul parametro m in blocchi
	for m_block in 1:block_size:mx
		m_end = min(m_block + block_size - 1, mx)

		Threads.@threads for m in m_block:m_end
			Lp_x[m, m] = Compute_Lp_Self2(l_s[m], W_s[m], T_s[m]) * escalings[:Lp]

			for n in m+1:mx
				integ, _ = Song_improved_Ivana_strategy2(
                    x_vs[m], y_vs[m], z_vs[m], xc_s[m], yc_s[m], zc_s[m], a_s[m], b_s[m], c_s[m], V_s[m],
                    x_vs[n], y_vs[n], z_vs[n], xc_s[n], yc_s[n], zc_s[n], a_s[n], b_s[n], c_s[n], V_s[n],
                    epsilon1, epsilon2, epsilon3, epsilon4, use_suppression,
                )
				value = 1e-7 / (volumi[:S][m] * volumi[:S][n]) * integ * escalings[:Lp]
				Lp_x[m, n] = value
				Lp_x[n, m] = value
				if n % 20 == 0
					yield()  # Permette alla task dell'heartbeat di essere schedulata
				end
			end
		end

		# Alla fine di ogni blocco, cede il controllo per permettere ad altre operazioni
		sleep(0)
		if is_stop_requested(id)
			println("Simulazione $(id) interrotta per richiesta stop.")
			return nothing # O un altro valore che indica interruzione
		end
		println("block Lp4 : ", round(m_end / block_size), " / ", round(mx / block_size))
	end

	# Suddividiamo il ciclo sul parametro m in blocchi
	for m_block in 1:block_size:my
		m_end = min(m_block + block_size - 1, my)

		Threads.@threads for m in m_block:m_end
			idx_m_y = m + mx # Indice assoluto per i volumi Y
			Lp_y[m, m] = Compute_Lp_Self2(l_s[idx_m_y], W_s[idx_m_y], T_s[idx_m_y]) * escalings[:Lp]

			for n in m+1:my
				idx_n_y = n + mx # Indice assoluto per i volumi Y
				integ, _ = Song_improved_Ivana_strategy2(
                    x_vs[idx_m_y], y_vs[idx_m_y], z_vs[idx_m_y], xc_s[idx_m_y], yc_s[idx_m_y], zc_s[idx_m_y], a_s[idx_m_y], b_s[idx_m_y], c_s[idx_m_y], V_s[idx_m_y],
                    x_vs[idx_n_y], y_vs[idx_n_y], z_vs[idx_n_y], xc_s[idx_n_y], yc_s[idx_n_y], zc_s[idx_n_y], a_s[idx_n_y], b_s[idx_n_y], c_s[idx_n_y], V_s[idx_n_y],
                    epsilon1, epsilon2, epsilon3, epsilon4, use_suppression,
                )
				value = 1e-7 / (volumi[:S][m+mx] * volumi[:S][n+mx]) * integ * escalings[:Lp]
				Lp_y[m, n] = value
				Lp_y[n, m] = value
				if n % 20 == 0
					yield()  # Permette alla task dell'heartbeat di essere schedulata
				end
			end
		end

		# Alla fine di ogni blocco, cede il controllo per permettere ad altre operazioni
		sleep(0)
		if is_stop_requested(id)
			println("Simulazione $(id) interrotta per richiesta stop.")
			return nothing # O un altro valore che indica interruzione
		end
		println("block Lp5 : ", round(m_end / block_size), " / ", round(my / block_size))
	end

	# Suddividiamo il ciclo sul parametro m in blocchi
	for m_block in 1:block_size:mz
		m_end = min(m_block + block_size - 1, mz)

		Threads.@threads for m in m_block:m_end
			idx_m_z = m + mx + my # Indice assoluto per i volumi Z
			Lp_z[m, m] = Compute_Lp_Self2(l_s[idx_m_z], W_s[idx_m_z], T_s[idx_m_z]) * escalings[:Lp]

			for n in m+1:mz
				idx_n_z = n + mx + my # Indice assoluto per i volumi Z
				integ, _ = Song_improved_Ivana_strategy2(
                    x_vs[idx_m_z], y_vs[idx_m_z], z_vs[idx_m_z], xc_s[idx_m_z], yc_s[idx_m_z], zc_s[idx_m_z], a_s[idx_m_z], b_s[idx_m_z], c_s[idx_m_z], V_s[idx_m_z],
                    x_vs[idx_n_z], y_vs[idx_n_z], z_vs[idx_n_z], xc_s[idx_n_z], yc_s[idx_n_z], zc_s[idx_n_z], a_s[idx_n_z], b_s[idx_n_z], c_s[idx_n_z], V_s[idx_n_z],
                    epsilon1, epsilon2, epsilon3, epsilon4, use_suppression,
                )
				value = 1e-7 / (volumi[:S][m+mx+my] * volumi[:S][n+mx+my]) * integ * escalings[:Lp]
				Lp_z[m, n] = value
				Lp_z[n, m] = value
                if n % 20 == 0
					yield()  # Permette alla task dell'heartbeat di essere schedulata
				end
			end
		end

		# Alla fine di ogni blocco, cede il controllo per permettere ad altre operazioni
		sleep(0)
		if is_stop_requested(id)
			println("Simulazione $(id) interrotta per richiesta stop.")
			return nothing # O un altro valore che indica interruzione
		end
		println("block Lp6 : ", round(m_end / block_size), " / ", round(mz / block_size))
	end

	return Dict(
		:Lp_x => Lp_x,
		:Lp_y => Lp_y,
		:Lp_z => Lp_z,
		:Rx => Rx,
		:Ry => Ry,
		:Rz => Rz,
	)
end
