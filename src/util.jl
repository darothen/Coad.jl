const ρw = 1000.0 # Density of water, kg/m³
const golovin_b = (1500.0)*1e-3 # Golovin collision coefficient, input cm³/g/s to -> m³/kg/s

"""
    mass_from_r(r; ρ=ρw)

Compute mass (kg) of a droplet with given radius (m) and density (kg/m³)
"""
mass_from_r(r; ρ=ρw) = (4*π/3)*ρ*(r^3)

"""
    r_from_mass(x; ρ=ρw)

Compute radius (m) of a droplet with given mass (kg) and density (kg/m³)
"""
r_from_mass(x; ρ=ρw) = (3*x/4/π/ρ)^(1/3)

# Initial cloud droplet distribution

## Collision Kernels

"""
    golovin_kernel(xᵢ, xⱼ)

Compute the Golovin collision kernel for droplets of two given masses (in kg)
"""
golovin_kernel(xᵢ, xⱼ) = golovin_b * (xᵢ + xⱼ)


function _interior_hydro_kernel(E_coal, E_coll, r_sum, tv_diff)
    E_coal * E_coll * π * r_sum*r_sum * abs(tv_diff)
end


"""
    hydrodynamic_kernel(xᵢ, xⱼ)

Compute the hydrodynamic collision kernel for droplets of the two given masses
(in kg) assuming unity collision and coalescence efficiency.
"""
function hydrodynamic_kernel(xᵢ, xⱼ)
    tvᵢ = terminal_v(xᵢ)
    tvⱼ = terminal_v(xⱼ)
    tv_diff = tvⱼ - tvᵢ
    r_sum = r_from_mass(xᵢ) + r_from_mass(xⱼ)
    _interior_hydro_kernel(1.0, 1.0, r_sum, tv_diff)
end

"""
    long_kernel(xᵢ, xⱼ)

Compute the hydrodynamic collision kernel with parameterization of collision
efficiency from Long (1974) for two droplets of given masses (in kg)
"""
function long_kernel(xᵢ, xⱼ)
    rᵢ = r_from_mass(xᵢ) 
    rⱼ = r_from_mass(xⱼ)
    tvᵢ = terminal_v(xᵢ)
    tvⱼ = terminal_v(xⱼ)

    tv_diff = tvⱼ - tvᵢ
    r_sum = rᵢ + rⱼ 

    r_small = min(rᵢ, rⱼ) * 1e6 # m -> micron
    r_large = max(rᵢ, rⱼ) * 1e6

    # Collection efficiency cut-off in limit of very large drops
    if r_large >= 50
        E_coll = 1.
    else
        E_coll = 4.5e-4 * r_large^2 * (1.0 - 3.0/(max(3.0, r_small) + 1e-2))
    end

    # Limit collection efficiency to 0 <= E_coll <= 1.0
    E_coll = min(E_coll, 1.)
    E_coll = max(0., E_coll)

    _interior_hydro_kernel(1.0, E_coll, r_sum, tv_diff)
end

"""
    terminal_v(x)

Compute terminal velocity of a droplet of given mass (in kg) following the
parameterization of Beard (1976)
"""
function terminal_v(x)
    r = r_from_mass(x)
    d = 2 * r * 1e6 # diameter, m -> μm
    x = x * 1e3 # mass, kg -> g

    if d ≤ 134.43
        α = 4.5795e5
        x_to_beta = x^(2 / 3)
    elseif 134.43 < d && d ≤ 1511.64
        α = 4962.0
        x_to_beta = x^(1 / 3)
    elseif 1511.64 < d && d ≤ 3477.84
        α = 1732.0
        x_to_beta = x^(1 / 6)
    else
        α = 917.0
        x_to_beta = 1.0
    end

    tv = 1e-2 * α * x_to_beta # cm/s -> m/s
end

"""
    kernels(x::AbstractArray{FT, 1}, kernel::Symbol) where {FT <: AbstractFloat}

Pre-compute and cache pair-wise collision kernels for a 1D binned discretization
of mass space.
"""
function kernels(x::AbstractArray{FT, 1}, kernel::Symbol) where {FT <: AbstractFloat}

  m = length(x)
  ck = zeros(FT, m, m)
  cck = zeros(FT, m, m)

  # Compute collision kernel for all potential bin interactions
  for j ∈ 1:m
    for i ∈ 1:j

      if kernel == :golovin
        cck[j, i] = golovin_kernel(x[i], x[j])
      elseif kernel == :hydro
        cck[j, i] = hydrodynamic_kernel(x[i], x[j])
      elseif kernel == :long
        cck[j, i] = long_kernel(x[i], x[j])
      else # default to Golovin
        cck[j, i] = golovin_kernel(x[i], x[j])
      end # kernel 

      # Copy over diagonal
      cck[i, j] = cck[j, i]
    end # i
  end # j

  # cache 2d interpolation on kernel vals
  for i ∈ 1:m
    for j ∈ 1:m
      jm = max(j - 1, 1)
      im = max(i - 1, 1)
      jp = min(j + 1, m)
      ip = min(i + 1, m)
      ck[i, j] = 0.125 * (
        cck[i, jm] + cck[im, j] + cck[ip, j] + cck[i, jp]
      )  + 0.5 * cck[i, j]
      if i == j
        ck[i, j] = 0.5 * ck[i, j]
      end
    end
  end

  return ck

end

