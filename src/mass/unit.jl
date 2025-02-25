struct UnitMass{N} <: AbstractMass{N} end

UnitMass(N::Int) = UnitMass{N}()

(m::UnitMass)(samples, ν = 0.0) = m

Base.rand(m::UnitMass) = m |> length |> randn

Base.:\(::UnitMass, x) = x

Base.:+(m::UnitMass, ::UnitMass) = m

LinearAlgebra.logabsdet(::UnitMass) = 0
