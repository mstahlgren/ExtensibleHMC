abstract type AbstractMass{N} end

Base.length(::AbstractMass{N}) where N = N

include("unit.jl")
#include("repdense.jl")
include("diag.jl")
#include("tridiag.jl")
#include("diagplus.jl")