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
α = 2^(1/4) # Mass bin scaling ratio, 2^(1/n) where 'n' is the number of bins 
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
# subroutine courant -- inplace

function courant(x)
  
  # Set up return matrices; first inspect length of our mass grid, 'x'
  m = length(x)
  c = zeros(Float64, m, m) # Courant numbers for limiting advection flux
  ima = zeros(Int16, m, m) # This is basically keeping track of the left bound of the 

  # Grid spacing in ln radius -> directly compute from the mass grid that was
  # passed in, since Δy = log(α) / 3
  α = x[2] / x[1]
  threeΔy = log(α)  # multiply term by 3 since we factor it out in the Courant calc

  # Compute Courant limits
  for i ∈ 1:m
    for j ∈ i:m

      x0 = x[i] + x[j]  # Summed / total mass from collision
      for k ∈ j:m
        if (x[k] ≥ x0) && (x[k-1] < x0) # There is probably an easier way to exploit the size of the collision here than linear searching for bounding masses
          if (c[i, j] < (1 - 1e-8))
            kk = k - 1
            c[i, j] = log(x0 / x[k-1]) / threeΔy
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
println("""
COMPUTE COURANT LIMITS ON GRID""")
CPUtic()
c, ima = courant(xᵢ)
elapsed = CPUtoq()
@printf "%3.1f seconds elapsed\n" elapsed


# Collision Kernel - we just use Golovin for now
# subroutine trkern
# cache kernel

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

function find_bounds(g; gmin=gmin)

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

function coad!(g, x, c, ck, ima; debug=false, gmin=gmin)

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