include("util.jl")

using CPUTime
using Plots
using Printf

const r₁ = (0.1)*1e-6 # minimum droplet size bin, micron->m
const rm = (100_000.0)*1e-6 # maximum droplet size bin, micron->m
# NOTE: We use subscript "1" because Julia is a 1-indexed language. By convention
# .     we reserve the nought subscript for general parameters.
const gmin = 1e-60 # lower bound on permissible bin mass density

## Parameters / Configuration
# Size Distribtion
r̅ = (10.0) * 1e-6 # mode radius of initial distribution, input micron convert to m
x̅ = mass_from_r(r̅) # mean initial droplet mass (kg)
L = (1.0)*1e-3 # total water content, g/m3 convert to kg/m3

# Grid Spacing
n = 4
α = 2^(1/n) # Mass bin scaling ratio, 2^(1/n) where 'n' is the number of bins 
            # between mass doubling

m = ceil(Integer, 1 + 3*log(rm / r₁)/log(α)) # number of mass bins to populate
Δy = Δlnr = log(α) / 3  # constant grid distance of logarithmic grid

# Collision kernel
# Options:
#   1) :golovin - Golovin (1963)
#   2) :hydro - basic hydrodynamic kernel with unit coal/coll efficiency
#   3) :long - Long collision efficiency
kernel = :long 

# Time Integration
tmax = 3601 # seconds
Δt = 5.0 # s
Δt_plot = 10 # minutes
nt = ceil(Integer, tmax / Δt)

# Other configs
debug = false
do_plots = true

# Report configuration to terminal
println("""
1D Stochastic Collection Equation Solver
========================================

Initial Conditions
------------------
L = $(L*1e3) g/m³ - Total LWC
r̅ = $(r̅*1e6) μm - Mode radius
x̅ = $x̅ kg - Mean droplet mass

Grid Setup
----------
$m bins spanning radii ($(r₁*1e6) μm - $(rm*1e6) μm)
α factor = $α (mass doubling every $n bins)

Time Integration
----------------
Integrating to $tmax seconds by Δt=$Δt

- DEBUG is $debug
- DO_PLOTS is $do_plots

""")

## Model Setup

# Mass and Radius Grid
#=
I find it a bit easier to reason through things in droplet-radius-space. We can
adopt the same convention of constructing a mass grid with an α scaling factor
(e.g. for α=√2 every mass bin is twice as big as the bin twice before it) and 
write simple, similar recurrence relationship in radius space:

  rᵢ = α rᵢ₋₁ 

Given r₀ as the smallest bin, we can then solve that

  rᵢ = r₁ α^{(i-1)/3} s.t. i ∈ ℤ > 0

Computing the number of grid cells...
  ln rm > ln(r1) + ln(α) (m-1)/3
  ln rm - ln r1 > ln α (m-1)/3
  3 ln(rm / r1) / (ln α) > m - 1

=# 
rᵢ = r₁*(α.^((collect(1:m) .- 1)./3)) # meter
xᵢ = mass_from_r.(rᵢ) # kg

# Initial droplet distribution; g(y, t) = 3x² * n(x, t), eq. 2 from B98 
gᵢ = (3 .* xᵢ.^2) .* nc.(xᵢ, L=L, x̅=x̅) # kg / m3

# Courant Numbers
#=
The final take-home here is that these courant limits are purely a function of
the mass grid discretization. On the one hand, they can be pre-computed and then
cached for retrieval as needed. This fact is used in combination with the fact
that the same holds true for the collision kernel to pre-compute all the limits
in the original code, as well as pre-estimate where collisions will land in the
actual coad subroutine
=#

println("""
COMPUTE COURANT LIMITS ON GRID""")
CPUtic()
c, ima = courant(xᵢ)
# Correct Courant numbers for grid spacing
c = c / (3*Δy)
elapsed = CPUtoq()
@printf "%3.1f seconds elapsed\n" elapsed

# Collision Kernel

println("""
COMPUTE COLLISION KERNELS ON GRID""")
CPUtic()
ck = kernels(xᵢ, kernel)
# Update kernel with contsant timestep and log grid distance
ck = ck.*(Δt*Δy)
elapsed = CPUtoq()
@printf "%3.1f seconds elapsed\n" elapsed

## Plot Initial Check

# Plot initial conditions and then begin the time loop
if do_plots
  println("Plotting initial conditions")
  p = plot(
    rᵢ*1e6, gᵢ*1e3, label="t = 0 min", 
    xaxis=:log, xlabel="r (μm)",
    # xlim=(0.5, 5000), xticks=[1, 10, 100, 1000],
    # ylim=(0, 0.9), yticks=0:0.06:0.9,
    ylabel="g (g / m³)",
  )
  display(p)
end

## TIME LOOP

#=
All of the time looping logic here is super old-school and can be totally
re-written for simplicity.
=#
# time integration
tlmin = 1e-6
t = 0.0
lmin = 0.0
total_runtime = 0.0

@inline function coad!(g, x, c, ck, ima; debug=false, gmin=gmin)

  # Lower and Upper integration limit i0, i1
  i0, i1 = find_bounds(g, gmin=gmin)

  if debug
    @printf "DEBUG bnds_check %6d %8d %8d\n" t i0 i1
  end

  # Main collision/coalescence loop
  for i ∈ i0:i1
    for j ∈ i:i1
      k = ima[i, j]  # Get pre-computed index of coalescence bin edge
      kp = k + 1

      # PORT - handle a weird edge condition in the initialization of ima?
      if k < 1
        continue
      end
    
      x0 = ck[i, j] * g[i] * g[j]
      x0 = min(x0, g[i] * x[j])
      if j != k # Not sure what's going on here.
        x0 = min(x0, g[j] * x[i])
      end
    
      gsi = x0 / x[j]
      gsj = x0 / x[i]
      gsk = gsi + gsj
      g[i] = g[i] - gsi
      g[j] = g[j] - gsj
      gk = g[k] + gsk

      if gk > gmin

        # ORIGINAL
        # x1 = log(gᵢ[kp] / gk + 1e-60)

        # MODIFIED - Apply a limiter to avoid negative args to log
        #   We note that this may not strictly obey the formulation of the flux
        #   algorithm, but this tends to work okay in practice.
        log_arg = g[kp] / gk + 1e-60
        log_arg = max(1e-60, log_arg)
        x1 = log(log_arg)

        flux = gsk / x1 * (exp(0.5 * x1) - exp(x1 * (0.5 - c[i, j])))
        flux = min(flux, gk)
        g[k] = gk - flux
        g[kp] = g[kp] + flux
      end
  
    end # j
  end # i

end

println("""
BEGIN TIME INTEGRATION
""")
CPUtic()
for i ∈ 1:nt
  global t = t + Δt
  global tlmin = tlmin + Δt

  # Collision
  # subroutine coad
  coad!(gᵢ, xᵢ, c, ck, ima, debug=debug, gmin=gmin)

  # Plotting
  if tlmin ≥ 60
    global tlmin = tlmin - 60
    global lmin = lmin + 1

    if (lmin % Δt_plot) < 1 & do_plots
      display(plot!(p, rᵢ*1e6, gᵢ*1e3, label = "t = $lmin min"))
    end

    # Mass balance? Not sure what's going on here. Maybe numerical checking?
    x0 = 0.0
    x1 = 1.0
    imax = 0
    for i ∈ 1:m
      x0 = x0 + gᵢ[i] * Δy
      x1 = max(x1, gᵢ[i])
      if abs(x1 - gᵢ[i]) < 1e-9
        imax = i
      end
    end
    @printf "    %4d mins |" lmin
    @printf " mass %10.3e  max %10.3e imax %3d" x0 x1 imax
    @printf "\n"
  end

  # Mass balance? Not sure what's going on here. Maybe numerical checking?
  x0 = 0.0
  x1 = 1.0
  imax = 0
  for i ∈ 1:m
    x0 = x0 + gᵢ[i] * Δy
    x1 = max(x1, gᵢ[i])
    if abs(x1 - gᵢ[i]) < 1e-9
      imax = i
    end
  end
  @printf "    %4d s |" t
  @printf " mass %10.3e  max %10.3e  imax %3d" x0 x1 imax
  @printf "\n"

  local elapsed = CPUtoq()
  global total_runtime += elapsed
  CPUtic()
end
@printf "Total time - %5.2f seconds\n" total_runtime

println("End")
xxx = readline()