abstract type AbstractMass{D} end

Base.length(m::AbstractMass) = prod(size(m))

include("unit.jl")
include("diag.jl")
include("robust.jl")
include("tridiag.jl")