const ρw = 1000.0 # Density of water, kg/m³
const golovin_b = (1500.0)*1e-3 # Golovin collision coefficient, input cm³/g/s to -> m³/kg/s

mass_from_r(r; ρ=ρw) = (4*π/3)*ρ*(r^3)
r_from_mass(x; ρ=ρw) = (3*x/4/π/ρ)^(1/3)

# Initial cloud droplet distribution
nc(x; L=L, x̅=x̅) = (L / x̅^2) * exp(-x / x̅)
# NOTE: this probably needs to be some sort of struct with an abstract type?
Exponential(x, L, x̅) = nc(x, L, x̅)

## Collision Kernels

# inputs in kg
golovin_kernel(xᵢ, xⱼ) = golovin_b * (xᵢ + xⱼ)

function _interior_hydro_kernel(E_coal, E_coll, r_sum, tv_diff)
    E_coal * E_coll * π * r_sum*r_sum * abs(tv_diff)
end

function hydrodynamic_kernel(xᵢ, xⱼ)
    tvᵢ = terminal_v(xᵢ)
    tvⱼ = terminal_v(xⱼ)
    tv_diff = tvⱼ - tvᵢ
    r_sum = r_from_mass(xᵢ) + r_from_mass(xⱼ)
    _interior_hydro_kernel(1.0, 1.0, r_sum, tv_diff)
end

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

# Other Functions
function terminal_v(x)
    # Beard (1976)
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

# Pre-compute collison kernels
function kernels(x, kernel)

  m = length(x)
  ck = zeros(Float64, m, m)
  cck = zeros(Float64, m, m)

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

# Pre-compute Courant limits
function courant(x)
  
  # Set up return matrices; first inspect length of our mass grid, 'x'
  m = length(x)
  c = zeros(Float64, m, m) # Courant numbers for limiting advection flux
  ima = zeros(Int16, m, m) # This is basically keeping track of the left bound of the 

  # Compute Courant limits
  for i ∈ 1:m
    for j ∈ i:m

      x0 = x[i] + x[j]  # Summed / total mass from collision
      for k ∈ j:m
        if (x[k] ≥ x0) && (x[k-1] < x0) # There is probably an easier way to exploit the size of the collision here than linear searching for bounding masses
          if (c[i, j] < (1 - 1e-8))
            kk = k - 1
            c[i, j] = log(x0 / x[k-1])
          else
            c[i, j] = 0.0
            kk = k
          end
          ima[i, j] = min(m - 1, kk)
          break
        end
      end # k loop

      # Copy over diagonal for symmetry
      c[j, i] = c[i, j]
      ima[j, i] = ima[i, j]

    end # j loop
  end # i loop

  return c, ima
end

@inline function find_bounds(g; gmin=1e-60)

  m = length(g)

  # TODO: refactor since this can be wrapped in a single array function
  # This basically sets a "focus" in the array where we have mass that needs
  # to get collided / advected around, so we don't waste cycles on empty bins.
  # In practice seems to be a limiter on numerical issues.
  i0 = 1
  for i ∈ 1:m-1
    i0 = i
    if g[i] > gmin
      break
    end
  end
  i1 = m-1
  for i ∈ m-1:-1:1
    i1 = i
    if g[i] > gmin
      break
    end
  end

  return i0, i1
end

