using MKL
using Interpolations

function build_trapezoidal_pulse(initial_delay_time::Float64, A::Float64, high_level_time::Float64, raise_time::Float64, falling_time::Float64, time_vect)

    # Find indices
    index_start_rise_time = findfirst(time -> time >= initial_delay_time, time_vect)
    index_end_rise_time = findfirst(time -> time >= initial_delay_time + raise_time, time_vect)
    index_start_falling_time = findfirst(time -> time >= initial_delay_time + raise_time + high_level_time, time_vect)
    index_end_falling_time = findfirst(time -> time >= initial_delay_time + raise_time + high_level_time + falling_time, time_vect)

    num_samples = length(time_vect)
    signal_time_sampled = zeros(num_samples) # Julia uses 1D arrays by default

    # High level
    if index_end_rise_time !== nothing && index_start_falling_time !== nothing
        start_idx = index_end_rise_time
        end_idx = index_start_falling_time
        if start_idx <= end_idx
            signal_time_sampled[start_idx:end_idx] .= A
        end
    end

    # Rising edge
    if index_start_rise_time !== nothing && index_end_rise_time !== nothing
        start_idx = index_start_rise_time
        end_idx = index_end_rise_time
        if start_idx <= end_idx
            # Handle boundary conditions for interpolation
            t_start_interp = time_vect[max(1, start_idx - 1)]
            t_end_interp = time_vect[min(num_samples, end_idx + 1)]
            y_start_interp = (start_idx > 1) ? 0.0 : 0.0 # Assume starts from 0 if index is the first
            y_end_interp = A

            itp_rise = linear_interpolation([t_start_interp, t_end_interp], [y_start_interp, y_end_interp])
            for i in start_idx:end_idx
                signal_time_sampled[i] = itp_rise(time_vect[i])
            end
        end
    end

    # Falling edge
    if index_start_falling_time !== nothing && index_end_falling_time !== nothing
        start_idx = index_start_falling_time
        end_idx = index_end_falling_time
        if start_idx <= end_idx
            # Handle boundary conditions for interpolation
            t_start_interp = time_vect[max(1, start_idx - 1)]
            t_end_interp = time_vect[min(num_samples, end_idx + 1)]
            y_start_interp = A
            y_end_interp = (end_idx < num_samples) ? 0.0 : 0.0 # Assume goes to 0 if index is the last

            itp_fall = linear_interpolation([t_start_interp, t_end_interp], [y_start_interp, y_end_interp])
            for i in start_idx:end_idx
                signal_time_sampled[i] = itp_fall(time_vect[i])
            end
        end
    end

    return signal_time_sampled
end