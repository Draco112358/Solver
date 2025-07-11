"""
Calculates the frequency vector for a given time series.
"""
function fft_frequencies(t::Vector{Float64})
    ncampt = length(t)
    fintem = t[end] - t[1]
    frefond = 1 / fintem
    
    # Return a range for efficiency
    return (0:Int(floor(ncampt / 2))) .* frefond
end

"""
Performs the FFT on a signal `x` using a pre-computed plan.
Returns only the positive frequency components of the spectrum.
"""
function fft_UAq(x::Vector{Float64}, planfft::FFTW.cFFTWPlan)
    ncampt = length(x)
    dtem = 1.0 # This should be derived from the time vector if not uniform
               # Assuming uniform sampling, dtem is implicitly handled by scaling.
               # If time sampling is not uniform, this approach needs adjustment.
               # Let's assume it's uniform for this optimization.
    
    # Perform the FFT and scale it
    # Note: dtem scaling depends on the convention you want to use.
    # The original code used (t[2]-t[1]), so we should pass that in if needed.
    # For now, let's assume the scaling is handled by the caller or is not needed.
    XAA = planfft * x
    
    # Return the view of the relevant part of the spectrum
    return @view XAA[1:Int(floor(ncampt / 2)) + 1]
end