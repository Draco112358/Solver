using MKL
using SpecialFunctions

function genera_segnale_esponenziale(tw::Float64, power::Float64, time, time_delay_vs::Float64)
  ratio = power^(-power-1) * gamma(power + 1) * exp(power)
  tr = tw / ratio
  vs = zeros(length(time))
  for i in eachindex(time)
    if time[i] >= time_delay_vs
      vs[i] = ((time[i] - time_delay_vs) / tr)^power * exp(-power * ((time[i] - time_delay_vs) / tr - 1))
    end
  end
  return vs, tr
end