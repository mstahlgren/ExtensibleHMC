struct UnitMass <: Mass end

Base.:*(::UnitMass, x) = x

Base.:\(::UnitMass, x) = x

