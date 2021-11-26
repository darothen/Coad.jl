export Timestamp

using Printf

"""
    struct Timestamp

A simple time-keeping struct with space for minutes and (fractional) seconds.
"""
mutable struct Timestamp
    minutes::Int64
    seconds::Real

    function Timestamp(minutes::Int64, seconds::Real)
        new(minutes, seconds)
    end

    function Timestamp(seconds::Real)
        seconds_over = mod(seconds, 60)
        minutes = floor(Int64, (seconds - seconds_over) / 60)
        new(minutes, seconds_over)
    end
end

function Base.show(io::IO, ts::Timestamp)
    seconds_fmt = @sprintf "%05.2f" ts.seconds
    print(io, "$(ts.minutes) mins $(seconds_fmt) secs")
end 

function Base.:+(x::Timestamp, y::Timestamp)
    total_minutes = x.minutes + y.minutes
    total_seconds = x.seconds + y.seconds
    if total_seconds ≥ 60
        total_minutes += 1
        total_seconds -= 60
    end
    return Timestamp(total_minutes, total_seconds)
end

function Base.:+(x::Timestamp, y::Real) 
    ts2 = Timestamp(y)
    return x + ts2
end
# TODO: Figure out why we don't automatically get commutation here...
Base.:+(x::Real, y::Timestamp) = y + x

function Base.:/(x::Timestamp, y::Real)
    total_seconds = x.minutes*60 + x.seconds
    divided_seconds = total_seconds / y
    return Timestamp(divided_seconds)
end

function Base.:%(x::Timestamp, y::Real)
    total_seconds = x.minutes*60 + x.seconds
    modulo = total_seconds % y
    return Timestamp(modulo)
end

function Base.:(==)(x::Timestamp, y::Timestamp)
    return (x.minutes == y.minutes) && (x.seconds == y.seconds)
end

function Base.:(<)(x::Timestamp, y::Timestamp)
    x_tot_s = x.minutes*60 + x.seconds
    y_tot_s = y.minutes*60 + y.seconds
    return x_tot_s < y_tot_s
end

function Base.:(<=)(x::Timestamp, y::Timestamp)
    return (x == y) || (x < y)
end