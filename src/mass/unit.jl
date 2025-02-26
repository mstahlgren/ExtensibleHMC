struct UnitMass{D} <: AbstractMass{D}
    S::NTuple{D, Int}
end

UnitMass(size...) = UnitMass{length(size)}(size)

(m::UnitMass)(samples, ν = 0.0) = m

Base.rand(m::UnitMass) = m |> size |> randn

Base.:\(::UnitMass, x) = x

Base.:+(m::UnitMass, ::UnitMass) = m

LinearAlgebra.logabsdet(::UnitMass) = 0
