function genera_punti_circonferenza(r::Float64, N::Int, baricentro::Vector{Float64}, cases::Int)
    angolo = collect(LinRange(0, 2*pi, N))
  
    A = r * cos.(angolo)
    B = r * sin.(angolo)
    C = zeros(N)
  
    if cases == 1 # circonferenza nel piano xy
      punti = transpose(baricentro) .+ hcat(A, B, C)
    elseif cases == 2 # circonferenza nel piano zx
      punti = transpose(baricentro) .+ hcat(B, C, A)
    else # circonferenza nel piano yz (corrected based on variable order in MATLAB's else)
      punti = transpose(baricentro) .+ hcat(C, A, B)
    end
  
    return punti
  end