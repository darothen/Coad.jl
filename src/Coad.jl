module Coad

export
  Coad1d,

  Exponential,

  Intergator1D,

  set!, run!

using CPUTime
using Printf

"""
    AbstractCoadModel

A super type for simulating collision/coalescence problems.
"""
abstract type AbstractCoadModel end

abstract type AbstractSizeDist end

include("util.jl")
include("size_dists.jl")
include("coad1d.jl")


end # module