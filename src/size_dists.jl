export ExponentialDist, nc

"""
    ExponentialDist{FT} <: AbstractSizeDist{FT}

Exponential size distribution defined by two moments: total water content, L̅ and
mean droplet mass, x̅
"""
struct ExponentialDist{FT} <: AbstractSizeDist{FT}
    L̅::FT # Total water content, kg/m3
    x̅::FT # Mean initial droplet mass, kg

    # TODO: Add special constructor to enforce positive args
    function ExponentialDist(;L̅::FT, x̅::FT) where {FT <: Real}
        new{FT}(L̅, x̅)
    end
end

function Base.show(io::IO, d::ExponentialDist{FT}) where {FT <: Real}
    print(io, "ExponentalDist(L̅=$(d.L̅ * 1e3) g/m³, x̅=$(d.x̅) kg)")
end 

"""
    nc(dist::ExponentialDist{FT}, x::FT) where {FT <: Real}

Compute the number concentration for droplets of the given mass for a 
specific size distribution.
"""
function nc(dist::ExponentialDist{FT}, x::FT) where {FT <: Real}
    return (dist.L̅ / dist.x̅^2) * exp(-x / dist.x̅)
end