abstract type Mass end

Base.:*(m::Mass, x) = m.McL * x

Base.:\(m::Mass, x) = m.M⁻¹ * x

include("unitmass.jl")
include("densemass.jl")