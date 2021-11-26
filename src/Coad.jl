module Coad

using CPUTime
using Printf

"""
    AbstractCoadModel

A super type for simulating collision/coalescence problems.
"""
abstract type AbstractCoadModel end

"""
    AbstractSizeDist{FT}

A super type for encapsulating size distribution calculations.
"""
abstract type AbstractSizeDist{FT} end

include("util.jl")
include("size_dists.jl")
include("time.jl")
include("coad1d.jl")

export
  Coad1D, set!, step!, ExponentialDist, nc,

  mass_from_r, r_from_mass

end # module