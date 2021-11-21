
const r₁ = (0.1)*1e-6 # minimum droplet size bin, micron->m
const rm = (100_000.0)*1e-6 # maximum droplet size bin, micron->m
# NOTE: We use subscript "1" because Julia is a 1-indexed language. By convention
# .     we reserve the nought subscript for general parameters.
const gmin = 1e-60 # lower bound on permissible bin mass density

"""
    Coad1D(
            f!::F!,
            j!::J!,
            x_init::A,
        ) where {F! <: Function, J! <: Function, A <: AbstractArray}

A non-linear system of equations type.

"""
struct Coad1D <: AbstractCoadModel
    # Grid spacing
    α,
    # Collison kernel
    kernel,

    # Grid definition
    xᵢ, rᵢ, Δlnr
    # Do we need to include anything else here?    
    function Coad1D(n, k)
        α = 2^(1/n)

        m = ceil(Integer, 1 + 3*log(rm / r₁)/log(α)) # number of mass bins to populate
        Δlnr = log(α) / 3  # constant grid distance of logarithmic grid
        
        rᵢ = r₁*(α.^((collect(1:m) .- 1)./3)) # meter
        xᵢ = mass_from_r.(rᵢ) # kg

        return new(α, k, xᵢ, rᵢ, Δlnr)
    end
end

