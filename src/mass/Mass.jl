abstract type AbstractMass{D} end

Base.size(m::AbstractMass) = m.size

Base.length(m::AbstractMass) = prod(size(m))

include("unit.jl")
#include("repdense.jl")
include("diag.jl")
#include("tridiag.jl")
#include("diagplus.jl")