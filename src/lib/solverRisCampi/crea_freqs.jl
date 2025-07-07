using MKL
function crea_freqs(t)
    fintem = t[end] - t[1]  # Lunghezza finestra temporale
    ncampt = length(t)      # N. di camp. nel tempo
  
    frefond = 1 / fintem
    f = (0:floor(ncampt / 2)) * frefond
    return f
  end