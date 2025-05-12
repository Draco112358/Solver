include("Compute_Lp_Self2.jl")
include("Song_improved_Ivana_strategy2.jl")

function calcola_Lp2(volumi, incidence_selection, escalings, QS_Rcc_FW)
	epsilon1 = 5e-3
	epsilon2 = 1e-3
	epsilon3 = 1e-3
	epsilon4 = 3e-1

	use_suppression = true

	mx = incidence_selection[:mx]
	my = incidence_selection[:my]
	mz = incidence_selection[:mz]

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
			Lp_x[m, m] = Compute_Lp_Self2(volumi[:coordinate][m, :], 1) * escalings[:Lp]

			for n in m+1:mx
				integ, _ = Song_improved_Ivana_strategy2(
					volumi[:coordinate][m, :], volumi[:coordinate][n, :],
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
		println("block Lp4 : ", round(m_end / block_size), " / ", round(mx / block_size))
	end

	# Suddividiamo il ciclo sul parametro m in blocchi
	for m_block in 1:block_size:my
		m_end = min(m_block + block_size - 1, my)

		Threads.@threads for m in m_block:m_end
			Lp_y[m, m] = Compute_Lp_Self2(volumi[:coordinate][m+mx, :], 2) * escalings[:Lp]

			for n in m+1:my
				integ, _ = Song_improved_Ivana_strategy2(
					volumi[:coordinate][m+mx, :], volumi[:coordinate][n+mx, :],
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
		println("block Lp5 : ", round(m_end / block_size), " / ", round(my / block_size))
	end

	# Suddividiamo il ciclo sul parametro m in blocchi
	for m_block in 1:block_size:mz
		m_end = min(m_block + block_size - 1, mz)

		Threads.@threads for m in m_block:m_end
			Lp_z[m, m] = Compute_Lp_Self2(volumi[:coordinate][m+mx+my, :], 3) * escalings[:Lp]

			for n in m+1:mz
				integ, _ = Song_improved_Ivana_strategy2(
					volumi[:coordinate][m+mx+my, :], volumi[:coordinate][n+mx+my, :],
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
