function genera_segnale_sinusoidale(f0::Float64, time, time_delay_vs::Float64)
  vs = zeros(length(time))
  for i in eachindex(time)
    if time[i] >= time_delay_vs
      vs[i] = cos(2 * pi * f0 * (time[i] - time_delay_vs))
    end
  end
  return vs
end