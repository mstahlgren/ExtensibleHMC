struct UnitMass{N} <: AbstractMass end

UnitMass(N::Int) = UnitMass{N}()

(m::UnitMass)(samples, ν = 0.0) = m

Base.rand(m::UnitMass{N}) where N = randn(N)

Base.:\(::UnitMass, x) = x

Base.:+(m::UnitMass, y::UnitMass) = m

LinearAlgebra.logabsdet(::UnitMass) = 0
