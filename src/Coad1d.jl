
const r₁ = (0.1)*1e-6 # minimum droplet size bin, micron->m
const rm = (100_000.0)*1e-6 # maximum droplet size bin, micron->m
# NOTE: We use subscript "1" because Julia is a 1-indexed language. By convention
# .     we reserve the nought subscript for general parameters.
const gmin = 1e-60 # lower bound on permissible bin mass density

"""
    Coad1D()

A non-linear system of equations type.

"""
struct Coad1D <: AbstractCoadModel
    # Grid spacing
    α,
    # Collison kernel
    kernel,

    # Grid definition
    xᵢ, rᵢ, Δlnr,

    # Initial conditions - don't set by default
    gᵢ, 

    # Cached computatonal values
    c, ima, # courant limits and collison idx LUT
    ck, # smoothed collision kernel LUT

    # Do we need to include anything else here?    
    function Coad1D(n, k)
        α = 2^(1/n)

        m = ceil(Integer, 1 + 3*log(rm / r₁)/log(α)) # number of mass bins to populate
        Δlnr = log(α) / 3  # constant grid distance of logarithmic grid
        
        rᵢ = r₁*(α.^((collect(1:m) .- 1)./3)) # meter
        xᵢ = mass_from_r.(rᵢ) # kg
        gᵢ = similar(xᵢ) .* 0

        c, ima = courant(xᵢ)
        ck = kernels(xᵢ, k)

        return new(α, k, xᵢ, rᵢ, Δlnr, gᵢ, c, ima, ck)
    end
end

function set!(model::Coad1D, dist::AbstractSizeDist)
  # Initial droplet distribution; g(y, t) = 3x² * n(x, t), eq. 2 from B98 
  model.gᵢ = (3 .* model.xᵢ.^2) .* nc.(dist, model.xᵢ) # kg / m3
end

@inline function step!(model::Coad1D, Δt)

  # Lower and Upper integration limit i0, i1
  i0, i1 = find_bounds(model.gᵢ, gmin=gmin)

  # if debug
  #   @printf "DEBUG bnds_check %6d %8d %8d\n" t i0 i1
  # end

  # Main collision/coalescence loop
  for i ∈ i0:i1
    for j ∈ i:i1
      k = model.ima[i, j]  # Get pre-computed index of coalescence bin edge
      kp = k + 1

      # PORT - handle a weird edge condition in the initialization of ima?
      if k < 1
        continue
      end
    
      # We did not pre-scale the collison kernels LUT by Δt and Δy in this 
      # version so we have to account for them here.
      x0 = model.ck[i, j] * model.gᵢ[i] * model.gᵢ[j] * Δt * model.Δlnr
      x0 = min(x0, model.gᵢ[i] * model.xᵢ[j])
      if j != k # Not sure what's going on here.
        x0 = min(x0, model.gᵢ[j] * model.xᵢ[i])
      end
    
      gsi = x0 / model.xᵢ[j]
      gsj = x0 / model.xᵢ[i]
      gsk = gsi + gsj
      model.gᵢ[i] = model.gᵢ[i] - gsi
      model.gᵢ[j] = model.gᵢ[j] - gsj
      gk = model.gᵢ[k] + gsk

      if gk > gmin

        # ORIGINAL
        # x1 = log(gᵢ[kp] / gk + 1e-60)

        # MODIFIED - Apply a limiter to avoid negative args to log
        #   We note that this may not strictly obey the formulation of the flux
        #   algorithm, but this tends to work okay in practice.
        log_arg = model.gᵢ[kp] / gk + 1e-60
        log_arg = max(1e-60, log_arg)
        x1 = log(log_arg)

        # We did not pre-scale the cached Courant limits for this version so we
        # have to apply them when estimating the flux directly.
        flux = gsk / x1 * (exp(0.5 * x1) - exp(x1 * (0.5 - (model.c[i, j]/3/Δy))))
        flux = min(flux, gk)
        model.gᵢ[k] = gk - flux
        model.gᵢ[kp] = model.gᵢ[kp] + flux
      end
  
    end # j
  end # i

end