function fft_UAq(t, x) # Assuming x can be complex
    fintem = t[end] - t[1]          # Length of the time window
    ncampt = length(t)              # Number of samples in time
    dtem = t[2] - t[1]              # Interval between two consecutive samples
  
    frecamp = 1 / dtem              # Sampling frequency
    fremax = frecamp / 2            # Nyquist frequency
    frefond = 1 / fintem            # Frequency resolution
    # Perform the FFT
    XAA = fft(x) * dtem
    XA = XAA[1:Int64(floor(ncampt/2)) + 1]
  
    # The MATLAB code comments out doubling the spectrum.
    # If you need to do this, uncomment the following line:
    # XA[2:end] = 2 * XA[2:end]
  
    f = (0:Int64(floor(ncampt/2))) * frefond
    Trasformata = [transpose(f); transpose(XA)] # Transpose to match MATLAB's column-wise output
    return Trasformata
  end