include("util.jl")

using LaTeXStrings
using Plots
using Printf

debug = false

const n = 251
# const emin = (1e-9) * 1e-6 # input mg convert to kg
const r₁ = (0.1)*1e-6  # minimum droplet size bin, micron->m
# NOTE: We use subscript "1" because Julia is a 1-indexed language. By convention
# .     we reserve the nought subscript for general parameters.
const tmax = 3601 # seconds
const Δt_plot = 30 # minutes
const gmin = 1e-60 # lower bound on permissible bin mass density

## Arrays
# Common arrays
# cour
c = zeros(Float64, n, n) # Courant numbers for limiting advection flux
ima = zeros(Int16, n, n) # This is basically keeping track of the left bound of the 
# bin that a particular collision on the mass grid will land in

# grid
g = zeros(Float64, n)
r = zeros(Float64, n)
e = zeros(Float64, n)

#kern 
ck = zeros(Float64, n, n)
ec = zeros(Float64, n, n)
cck = zeros(Float64, n, n)

## Parameters
r̅ = (10.0) * 1e-6 # mode radius of initial distribution, input micron convert to m
x̅ = mass_from_r(r̅) # mean initial droplet mass (kg)
L = (1.0)*1e-3 # total water content, g/m3 convert to kg/m3

α = 2^(1/4) # Mass bin scaling ratio, 2^(1/n) where 'n' is the number of bins 
            # between mass doubling
Δy = Δlnr = log(α) / 3  # constant grid distance of logarithmic grid

# Mass and Radius Grid
#=
I find it a bit easier to reason through things in droplet-radius-space. We can
adopt the same convention of constructing a mass grid with an α scaling factor
(e.g. for α=√2 every mass bin is twice as big as the bin twice before it) and 
write simple, similar recurrence relationship in radius space:

  rᵢ = α rᵢ₋₁ 

Given r₀ as the smallest bin, we can then solve that

  rᵢ = r₁ α^{(i-1)/3} s.t. i ∈ ℤ > 0

=# 
rᵢ = r₁*(α.^((collect(1:n) .- 1)./3)) # meter
xᵢ = mass_from_r.(rᵢ) # kg

# Initial droplet distribution; g(y, t) = 3x² * n(x, t), eq. 2 from B98 
gᵢ = (3 .* xᵢ.^2) .* nc.(xᵢ, L=L, x̅=x̅) # kg / m3

# for i ∈ 1:n
#   e[i] = xᵢ[i]
#   r[i] = rᵢ[i]
#   g[i] = g2[i]
# end

# Sanity check
if debug
  for i ∈ 1:n
    @printf "%10d %24g %24g %24g\n" i xᵢ[i]*1e6 rᵢ[i] gᵢ[i]
  end
end

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
for i ∈ 1:n
  for j ∈ i:n

    local x0 = xᵢ[i] + xᵢ[j]  # Summed / total mass from collision
    for k ∈ j:n
      if (xᵢ[k] ≥ x0) && (xᵢ[k-1] < x0) # There is probably an easier way to exploit the size of the collision here than linear searching for bounding masses
        if (c[i, j] < (1 - 1e-8))
          kk = k - 1
          c[i, j] = log(x0 / xᵢ[k-1]) / (3 * Δy)
        else
          c[i, j] = 0.0
          kk = k
        end
        # @printf "ima | %3d %3d | %3d %3d \n" i j n-1 kk
        ima[i, j] = min(n - 1, kk)
        # ima[i, j] = kk
        break
      end
    end # k loop

    # Copy over diagonal for symmetry
    c[j, i] = c[i, j]
    ima[j, i] = ima[i, j]

  end # j loop
end # i loop

for i ∈ 2:n
  isq = Integer(i*i)
  isq_2 = Integer(floor(isq/2))
  if (isq > n) 
    break
  end
  @printf "%8d %8d %10g %8d\n" isq isq_2 c[isq, isq_2] ima[isq, isq_2]
end
for j ∈ 15:n
  @printf "%8d %8d %10g %8d\n" 20 j c[20, j] ima[20, j]
end

# Collision Kernel - we just use Golovin for now
# subroutine trkern
# cache kernel
for j ∈ 1:n
  for i ∈ 1:j
    cck[j, i] = golovin_kernel(xᵢ[i], xᵢ[j])
    cck[i, j] = cck[j, i]
  end
end

# cache 2d interpolation on kernel vals
for i ∈ 1:n
  for j ∈ 1:n
    jm = max(j - 1, 1)
    im = max(i - 1, 1)
    jp = min(j + 1, n)
    ip = min(i + 1, n)
    ck[i, j] = 0.125 * (
      cck[i, jm] + cck[im, j] + cck[ip, j] + cck[i, jp]
    )  + 0.5 * cck[i, j]
    if i == j
      ck[i, j] = 0.5 * ck[i, j]
    end
  end
end

for i ∈ 2:n
  isq = Integer(i*i)
  isq_2 = Integer(floor(isq/2))
  if (isq > n) 
    break
  end
  @printf "%8d %8d %10g %10g\n" isq isq_2 cck[isq, isq_2] ck[isq, isq_2]
end

# Plot initial conditions and then begin the time loop
p = plot(
  rᵢ*1e6, gᵢ*1e3, label="t = 0 min", 
  xaxis=:log, xlim=(0.5, 5000), xlabel="r (μm)",
  xticks=[1, 10, 100, 1000],
  ylim=(0, 0.9), ylabel="g (g / m³)",
  yticks=0:0.06:0.9
)
display(p)

## TIME LOOP

Δt = 10.0 # s
nt = ceil(Integer, tmax / Δt)

# Update kernel with contsant timestep and log grid distance
for i ∈ 1:n
  for j ∈ 1:n
    ck[i, j] = ck[i, j]*Δt*Δy
  end
end

for i ∈ 2:n
  isq = Integer(i*i)
  isq_2 = Integer(floor(isq/2))
  if (isq > n) 
    break
  end
  @printf "%8d %8d %10g \n" isq isq_2 ck[isq, isq_2]
end


#=
All of the time looping logic here is super old-school and can be totally
re-written for simplicity.
=#
# time integration
tlmin = 1e-6
t = 0.0
lmin = 0.0
for i ∈ 1:nt
  global t = t + Δt
  global tlmin = tlmin + Δt

  # Collision
  # subroutine coad

  # Lower and Upper integration limit i0, i1
  # TODO: refactor since this can be wrapped in a single array function
  # This basically sets a "focus" in the array where we have mass that needs
  # to get collided / advected around, so we don't waste cycles on empty bins.
  # In practice seems to be a limiter on numerical issues.
  i0 = 1
  for i ∈ 1:n-1
    i0 = i
    if gᵢ[i] > gmin
      break
    end
  end
  i1 = n-1
  for i ∈ n-1:-1:1
    i1 = i
    if gᵢ[i] > gmin
      break
    end
  end

  @printf "bnds_check %6d %8d %8d\n" t i0 i1

  # Main collision/coalescence loop
  for i ∈ i0:i1
    for j ∈ i:i1
      k = ima[i, j]  # Get pre-computed index of coalescence bin edge
      kp = k + 1

      # PORT - handle a weird edge condition in the initialization of ima?
      if k < 1
        continue
      end
    
      local x0 = ck[i, j] * gᵢ[i] * gᵢ[j]
      x0 = min(x0, gᵢ[i] * xᵢ[j])
      if j != k # Not sure what's going on here.
        local x0 = min(x0, gᵢ[j] * xᵢ[i])
      end
    
      gsi = x0 / xᵢ[j]
      gsj = x0 / xᵢ[i]
      gsk = gsi + gsj
      gᵢ[i] = gᵢ[i] - gsi
      # if g[i] < 0
      #   @printf "WARNING - g[%d] = %e < 0 | %e \n" i g[i] gsi
      # end
      gᵢ[j] = gᵢ[j] - gsj
      # if g[j] < 0
      #   @printf "WARNING - g[%d] = %e < 0 | %e \n" j g[j] gsj
      # end
      gk = gᵢ[k] + gsk

      # @printf "a | (%3d, %3d) %13.6e %13.6e %13.6e\n" i j gsi gsj gsk

      if gk > gmin
        x1 = log(gᵢ[kp] / gk + 1e-60)
        # @printf "x | %13.6e %13.6e %13.6e %13.6e\n" g[kp] g[k] gk x1
        flux = gsk / x1 * (exp(0.5 * x1) - exp(x1 * (0.5 - c[i, j])))
        # @printf "b | %13.6e %13.6e %13.6e %13.6e\n" x1 gk flux gsk
        flux = min(flux, gk)
        # @printf " %3d %3d %e %e\n" i j gk flux
        gᵢ[k] = gk - flux
        gᵢ[kp] = gᵢ[kp] + flux
      end
    
    end # j
  end # i

  # Plotting
  if tlmin ≥ 60
    global tlmin = tlmin - 60
    global lmin = lmin + 1

    if (lmin % Δt_plot) < 1
      display(plot!(p, rᵢ*1e6, gᵢ*1e3, label = "t = $lmin min"))
    end

    # Mass balance? Not sure what's going on here. Maybe numerical checking?
    x0 = 0.0
    x1 = 1.0
    imax = 0
    for i ∈ 1:n
      x0 = x0 + gᵢ[i] * Δy
      x1 = max(x1, gᵢ[i])
      if abs(x1 - gᵢ[i]) < 1e-9
        imax = i
      end
    end
    @printf "    %4d mins |" lmin
    @printf " mass %3.2e  max %3.2e  imax %3d" x0 x1 imax
    @printf "\n"
  end

  # Mass balance? Not sure what's going on here. Maybe numerical checking?
  x0 = 0.0
  x1 = 1.0
  imax = 0
  for i ∈ 1:n
    x0 = x0 + gᵢ[i] * Δy
    x1 = max(x1, gᵢ[i])
    if abs(x1 - gᵢ[i]) < 1e-9
      imax = i
    end
  end
  @printf "    %4d s |" t
  @printf " mass %3.2f  max %3.2f  imax %3d" x0 x1 imax
  @printf "\n"

end

println("End")
xxx = readline()