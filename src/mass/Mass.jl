abstract type AbstractMass{D} end

Base.length(m::AbstractMass) = prod(size(m))

include("unit.jl")
include("diag.jl")
#include("repdense.jl")
#include("tridiag.jl")
#include("diagplus.jl")