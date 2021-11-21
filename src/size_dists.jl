
struct ExponentialDist{R} <: AbstractSizeDist
    L̅::R, # Total water content, kg/m3
    x̅::R # Mean initial droplet mass, kg
end

function nc(dist::ExponentialDist{R}, x::R)
    return (dist.L / dist.x̅^2) * exp(-x / x̅)
end