abstract type Mass{T} end

Base.eltype(::Mass{T}) where T = T

Base.:*(m::Mass, x) = m.McL * x

Base.:\(m::Mass, x) = m.M⁻¹ * x

include("unitmass.jl")
include("densemass.jl")