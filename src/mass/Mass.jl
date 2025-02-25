abstract type AbstractMass end

Base.length(m::AbstractMass) = m.R * m.C

include("unit.jl")
include("repdense.jl")
include("diag.jl")
include("tridiag.jl")
include("diagplus.jl")