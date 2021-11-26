module Coad

#   #Integrator1D,

#   set!#, run!

using CPUTime
using Printf

"""
    AbstractCoadModel

A super type for simulating collision/coalescence problems.
"""
abstract type AbstractCoadModel end

abstract type AbstractSizeDist{FT} end

include("util.jl")
include("size_dists.jl")
include("time.jl")
include("coad1d.jl")

export
  Coad1D, set!, step!, ExponentialDist, nc,

  mass_from_r, r_from_mass


end # module