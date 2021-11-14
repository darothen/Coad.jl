using Plots
using Printf

const n = 201
const emin = (1e-9) * 1e-6 # input mg convert to kg
const tmax = 1801
const ρw = 1000.0 # Density of water, kg/m³
# TODO: Figure out why we need the scaling here....
#       it should be 1500 cm3/g/s = 1.5 m3/kg/s but we have the extra 1e3
const golovin_b = 1.5/1e3 # Golovin collision coefficient, input cm³/g/s to ≡ m³/kg/s
const gmin = 1e-60 # lower bound on permissible bin mass density

## Arrays
arr1d = zeros(Float64, n)
arr2d = zeros(Float64, n, n)

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
rq0 = (10.0) * 1e-6 # mode radius of initial distribution, input micron convert to m
xmw = (1000.0) * 1e-3 # total water content, input g/m3 convert to kg/m3
scal = 4 # scaling factor of ax

xn0 = (4 * π / 3) * ρw * (rq0)^3 # mean initial droplet mass (kg)
xn1 = xmw / xn0 # total initial droplet number concentration, 1/m3

dlnr = log(2) / 3 / scal # constant grid distance of logarithmic grid
ax = 2^(1 / scal) # growth factor for consecutive masses

# Mass and Radius Grid
e[1] = emin * 0.5 * (ax + 1) # this should really be "x" - it's the mass grid; kg
# r[1] = 1000*(3*e[1]/4/π)^(1/3) # this is just computing radius from water droplet mass; m
# This is convolving water density with scaling factors for units. Need to make cleaner.
r[1] = (3 * e[1] / 4 / π / ρw)^(1 / 3)
for i ∈ 2:n
  e[i] = ax * e[i-1]
  r[i] = (3 * e[i] / 4 / π / ρw)^(1 / 3)
end

# Initial Mass Distribution
x0 = xn1 / xn0
for i = 1:n
  x1 = e[i]
  # g = 3x^2 * n(x, t)
  g[i] = (3 * x1^2) * (x0 * exp(-x1 / xn0))
end

# Sanity check
for i ∈ 1:n
  @printf "%10d %24g %24g %24g\n" i e[i]*1e6 r[i] g[i]
end
# xxx = readline()

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

    local x0 = e[i] + e[j]  # Summed / total mass from collision
    for k ∈ j:n
      if (e[k] ≥ x0) && (e[k-1] < x0) # There is probably an easier way to exploit the size of the collision here than linear searching for bounding masses
        if (c[i, j] < (1 - 1e-8))
          kk = k - 1
          c[i, j] = log(x0 / e[k-1]) / (3 * dlnr)
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
    cck[j, i] = golovin_b * (e[i] + e[j])
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
p = plot(r * 1e6, g, xaxis = :log, label = "t = 0 min", xlim=(0.5, 5000))
display(p)

## TIME LOOP

Δt = 1.0 # s
nt = Integer(tmax / Δt)

# Update kernel with contsant timestep and log grid distance
for i ∈ 1:n
  for j ∈ 1:n
    ck[i, j] = ck[i, j]*Δt*dlnr
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
    if g[i] > gmin
      break
    end
  end
  i1 = n-1
  for i ∈ n-1:-1:1
    i1 = i
    if g[i] > gmin
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
    
      local x0 = ck[i, j] * g[i] * g[j]
      x0 = min(x0, g[i] * e[j])
      if j != k # Not sure what's going on here.
        local x0 = min(x0, g[j] * e[i])
      end
    
      gsi = x0 / e[j]
      gsj = x0 / e[i]
      gsk = gsi + gsj
      g[i] = g[i] - gsi
      # if g[i] < 0
      #   @printf "WARNING - g[%d] = %e < 0 | %e \n" i g[i] gsi
      # end
      g[j] = g[j] - gsj
      # if g[j] < 0
      #   @printf "WARNING - g[%d] = %e < 0 | %e \n" j g[j] gsj
      # end
      gk = g[k] + gsk

      # @printf "a | (%3d, %3d) %13.6e %13.6e %13.6e\n" i j gsi gsj gsk

      if gk > gmin
        x1 = log(g[kp] / gk + 1e-60)
        # @printf "x | %13.6e %13.6e %13.6e %13.6e\n" g[kp] g[k] gk x1
        flux = gsk / x1 * (exp(0.5 * x1) - exp(x1 * (0.5 - c[i, j])))
        # @printf "b | %13.6e %13.6e %13.6e %13.6e\n" x1 gk flux gsk
        flux = min(flux, gk)
        # @printf " %3d %3d %e %e\n" i j gk flux
        g[k] = gk - flux
        g[kp] = g[kp] + flux
      end
    
    end # j
  end # i

  # Plotting
  if tlmin ≥ 60
    global tlmin = tlmin - 60
    global lmin = lmin + 1

    if (lmin % 10) < 1
      display(plot!(p, r * 1e6, g, xaxis = :log, label = "t = $lmin min"))
    end

    # Mass balance? Not sure what's going on here. Maybe numerical checking?
    x0 = 0.0
    x1 = 1.0
    imax = 0
    for i ∈ 1:n
      x0 = x0 + g[i] * dlnr
      x1 = max(x1, g[i])
      if abs(x1 - g[i]) < 1e-9
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
    x0 = x0 + g[i] * dlnr
    x1 = max(x1, g[i])
    if abs(x1 - g[i]) < 1e-9
      imax = i
    end
  end
  @printf "    %4d s |" t
  @printf " mass %3.2f  max %3.2f  imax %3d" x0 x1 imax
  @printf "\n"

end

