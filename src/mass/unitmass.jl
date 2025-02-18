struct UnitMass{T} <: Mass{T} end

UnitMass(T = Float64) = UnitMass{T}()

Base.:*(::UnitMass, x) = x

Base.:\(::UnitMass, x) = x

