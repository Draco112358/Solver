using LinearAlgebra

function get_punti_oss_3D(r::Float64, N_circ::Int, baricentro::Vector{Float64})
    theta = LinRange(0, pi, N_circ)      # theta from 0 to pi
    phi = LinRange(0, 2*pi, N_circ)      # phi from 0 to 2pi
  
    phi1 = [p for p in phi, t in theta]
    theta1 = [t for p in phi, t in theta]
  
    x_grid = baricentro[1] .+ r * sin.(theta1) .* cos.(phi1)
    y_grid = baricentro[2] .+ r * sin.(theta1) .* sin.(phi1)
    z_grid = baricentro[3] .+ r * cos.(theta1)
  
    # Generiamo tutte le possibili combinazioni di theta e phi senza ripetizioni
    theta_grid = [t for t in theta, p in phi]  # Genera la griglia di theta
    phi_grid = [p for t in theta, p in phi]    # Genera la griglia di phi
  
    # Convertiamo la griglia in una lista di punti unici
    theta_vals = vec(theta_grid)  # Otteniamo i valori di theta come vettore colonna
    phi_vals = vec(phi_grid)      # Otteniamo i valori di phi come vettore colonna
  
    # Creazione della matrice N x 3 per le coordinate cartesiane
    N = length(theta_vals)  # Numero totale di punti
    coordinates = zeros(N, 3)  # Matrice per le coordinate cartesiane (x, y, z)
    distanze = zeros(N)
    theta_vect = zeros(N) # It seems this variable is not used in the MATLAB code
  
    # Calcoliamo le coordinate cartesiane per ogni coppia di (theta, phi)
    for i in 1:N
      x = r * sin(theta_vals[i]) * cos(phi_vals[i])
      y = r * sin(theta_vals[i]) * sin(phi_vals[i])
      z = r * cos(theta_vals[i])
  
      coordinates[i, :] = baricentro + [x, y, z]  # assemblo le coordinate
      distanze[i] = norm([x, y, z] .- baricentro)
    end
  
    return coordinates, distanze, theta_vals, transpose(x_grid), transpose(y_grid), transpose(z_grid)
  end