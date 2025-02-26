struct UnitMass{S} <: AbstractMass{S} end

UnitMass(size...) = UnitMass{size}()

(m::UnitMass)(samples, ν = 0.0) = m

Base.rand(m::UnitMass) = m |> size |> randn

Base.:\(::UnitMass, x) = x

Base.:+(m::UnitMass, ::UnitMass) = m

LinearAlgebra.logabsdet(::UnitMass) = 0

testing(::UnitMass{S}) where S = S