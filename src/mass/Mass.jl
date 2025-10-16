import Random: randn!

abstract type AbstractMass end

function resize end

include("unit.jl")
include("diag.jl")
include("robust.jl")