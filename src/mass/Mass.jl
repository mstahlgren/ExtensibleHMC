abstract type AbstractMass{S} end

Base.size(::AbstractMass{S}) where S = S

Base.length(::AbstractMass{S}) where S = prod(S)

include("unit.jl")
#include("repdense.jl")
include("diag.jl")
#include("tridiag.jl")
#include("diagplus.jl")