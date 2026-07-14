"""
    interpolate_noise(time_grid, noise, t)

Linearly interpolate a discrete noise realization at time `t`.

# Arguments
- `time_grid`: uniformly spaced time grid.
- `noise`: noise values sampled on `time_grid`.
- `t`: interpolation time.

# Returns
- `ComplexF64`: interpolated noise value at `t`.

# Throws
- `AssertionError`: if `t` lies outside `time_grid`.
- `AssertionError`: if `noise` and `time_grid` have different lengths.
"""
@inline function interpolate_noise(
                                time_grid::TimeGrid, 
                                noise::Vector{ComplexF64}, 
                                t::Float64,
                            )
                            
    @assert t>=time_grid[1]&&t<=time_grid[end]  "Time value out of bounds for interpolation."
    @assert length(noise)==length(time_grid)    "Mismatch of size of discrete noise and grid."

    dt = time_grid[2]-time_grid[1]
    ind::Int = round(Int,(t-time_grid[1])/dt)+1
    if time_grid[ind]==t
        return noise[ind]
    elseif t-time_grid[ind] > 0
        return noise[ind+1]*(t-time_grid[ind])/dt+noise[ind]*(time_grid[ind+1]-t)/dt
    else
        return noise[ind]*(t-time_grid[ind-1])/dt+noise[ind-1]*(time_grid[ind]-t)/dt
    end
end