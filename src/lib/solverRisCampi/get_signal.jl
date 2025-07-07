include("genera_segnale_esponenziale.jl")
include("genera_segnale_Gaussiano_modulato.jl")
include("genera_segnale_sinusoidale.jl")
include("build_trapezoidal_pulse.jl")
using MKL

function getSignalbasedOn(signal_type, times)
    if signal_type["type"] == "exponential"
        f0=0
        dev_stand=0;

        time_delay_vs=parse(Float64,signal_type["params"]["time_delay_vs"]);
        tw = parse(Float64,signal_type["params"]["tw"])
        power=parse(Float64,signal_type["params"]["power"]);
        vs, tr = genera_segnale_esponenziale(tw, power, times, time_delay_vs)
        return vs
    elseif signal_type["type"] == "gaussian_modulated"
        tr=0;
        power=0;
        
        time_delay_vs=parse(Float64,signal_type["params"]["time_delay_vs"]);
        f0=parse(Float64,signal_type["params"]["f0"]);
        dev_stand=parse(Float64,signal_type["params"]["dev_stand"]);
        vs = genera_segnale_Gaussiano_modulato(f0, dev_stand, times, time_delay_vs)
        return vs
    elseif signal_type["type"] == "sinusoidal"
        tr=0;
        dev_stand=0;
        power=0;
        time_delay_vs=parse(Float64,signal_type["params"]["time_delay_vs"]);
        f0=parse(Float64,signal_type["params"]["f0"]);
        vs = genera_segnale_sinusoidale(f0, times, time_delay_vs)
        return vs
    elseif signal_type["type"] == "trapezoidal_pulse"
        Amplitude=parse(Float64,signal_type["params"]["A"]);
        initial_delay_time=parse(Float64,signal_type["params"]["initial_delay_time"]);
        high_level_time=parse(Float64,signal_type["params"]["high_level_time"]);
        raise_time=parse(Float64,signal_type["params"]["raise_time"]);
        falling_time=parse(Float64,signal_type["params"]["falling_time"]);
        vs = build_trapezoidal_pulse(initial_delay_time, Amplitude, high_level_time, raise_time, falling_time, times)
        return vs
    end
end