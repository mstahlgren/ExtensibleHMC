abstract type AbstractMass{D} end

Base.length(m::AbstractMass) = prod(size(m))

include("unit.jl")
#include("repdense.jl")
include("diag.jl")
#include("tridiag.jl")
#include("diagplus.jl")