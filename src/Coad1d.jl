
const r₁ = (0.1)*1e-6 # minimum droplet size bin, micron->m
const rm = (100_000.0)*1e-6 # maximum droplet size bin, micron->m
# NOTE: We use subscript "1" because Julia is a 1-indexed language. By convention
# .     we reserve the nought subscript for general parameters.
const gmin = 1e-60 # lower bound on permissible bin mass density

"""
    Coad1D()

A non-linear system of equations type.

Typing inspired by: https://github.com/CliMA/RRTMGP.jl/blob/master/src/optics/RTE.jl
"""
struct Coad1D{
  FT <: AbstractFloat,
  I <: Int, 
  FTA1D <: AbstractArray{FT, 1},
  FTA2D <: AbstractArray{FT, 2},
  IA2D <: AbstractArray{I, 2}
} <: AbstractCoadModel
    # Grid spacing
    α::FT
    # Collison kernel
    kernel::Symbol

    # Grid definition
    xᵢ::FTA1D
    rᵢ::FTA1D
    Δlnr::FT

    # Initial conditions - don't set by default
    gᵢ::FTA1D

    # Cached computatonal values
    c::FTA2D
    ima::IA2D # courant limits and collison idx LUT
    ck::FTA2D # smoothed collision kernel LUT

    # Do we need to include anything else here?    
    function Coad1D(;n::Int64, kernel::Symbol)
        α = 2^(1/n)

        m = ceil(Integer, 1 + 3*log(rm / r₁)/log(α)) # number of mass bins to populate
        Δlnr = log(α) / 3  # constant grid distance of logarithmic grid

        println("Instatiating model with $m mass bins...")
        rᵢ = r₁*(α.^((collect(1:m) .- 1)./3)) # meter
        xᵢ = mass_from_r.(rᵢ) # kg
        gᵢ = similar(xᵢ) .* 0

        println("  Pre-caching grid Courant limits...")
        c, ima = courant(xᵢ)
        println("  Pre-caching $kernel kernel lookup table...")
        ck = kernels(xᵢ, kernel)

        FT = eltype(xᵢ)
        I = eltype(ima)
        FTA1D = typeof(xᵢ)
        FTA2D = typeof(c)
        IA2D = typeof(ima)

        return new{
          FT, I, FTA1D, FTA2D, IA2D
        }(
          α, kernel, xᵢ, rᵢ, Δlnr, gᵢ, c, ima, ck
        )
    end
end

function Base.show(io::IO, model::Coad1D)
  m = length(model.xᵢ)
  print(io, "Coad1D Model ($m radii ∈ ($(r₁*1e6), $(rm*1e6)) μm, kernel=$(model.kernel))")
end 

# Pre-compute Courant limits
@inline function courant(x::AbstractArray{FT}) where {FT <: AbstractFloat} 
  
  # Set up return matrices; first inspect length of our mass grid, 'x'
  m = length(x)
  # FT = eltype(x)
  c = zeros(FT, m, m) # Courant numbers for limiting advection flux
  ima = zeros(Int64, m, m) # This is basically keeping track of the left bound of the 

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


function set!(model::Coad1D, dist::AbstractSizeDist{FT}) where {FT <: Real}
  # Initial droplet distribution; g(y, t) = 3x² * n(x, t), eq. 2 from B98 
  model.gᵢ .= (3 .* model.xᵢ.^2) .* nc.(Ref(dist), model.xᵢ) # kg / m3
end

@inline function step!(model::Coad1D, Δt::Float64)

  # Lower and Upper integration limit i0, i1
  i0, i1 = find_bounds(model.gᵢ, gmin=gmin)

  # if debug
  #   @printf "DEBUG bnds_check %6d %8d %8d\n" t i0 i1
  # end

  # Main collision/coalescence loop
  @inbounds for i ∈ i0:i1
    @inbounds for j ∈ i:i1
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
        flux = gsk / x1 * (exp(0.5 * x1) - exp(x1 * (0.5 - (model.c[i, j]/3/model.Δlnr))))
        flux = min(flux, gk)
        model.gᵢ[k] = gk - flux
        model.gᵢ[kp] = model.gᵢ[kp] + flux
      end
  
    end # j
  end # i

end

"""
    find_bounds(g; gmin)

Given a binned discretization of a droplet mass distribution, find the range of
bin indices which have non-negligble mass. This function is intended to help
accelerate other computations which require analyzing a process on the full mass
grid by minimizing the number of bins which need to be inspected.
"""
@inline function find_bounds(g; gmin=1e-60)

  m = length(g)

  # TODO: refactor since this can be wrapped in a single array function
  # This basically sets a "focus" in the array where we have mass that needs
  # to get collided / advected around, so we don't waste cycles on empty bins.
  # In practice seems to be a limiter on numerical issues.
  i0 = 1
  @inbounds for i ∈ 1:m-1
    i0 = i
    if g[i] > gmin
      break
    end
  end
  i1 = m-1
  @inbounds for i ∈ m-1:-1:1
    i1 = i
    if g[i] > gmin
      break
    end
  end

  return i0, i1
end