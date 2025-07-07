using MKL
function genera_segnale_Gaussiano_modulato(f0::Float64, dev_stand::Float64, time, time_delay_vs::Float64)
    vs = zeros(length(time))
    for i in eachindex(time)
      if time[i] >= time_delay_vs
        vs[i] = cos(2 * pi * f0 * (time[i] - time_delay_vs)) * exp(-(time[i] - time_delay_vs)^2 / (2 * dev_stand^2))
      end
    end
    return vs
  end