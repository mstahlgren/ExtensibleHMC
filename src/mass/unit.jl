struct UnitMass{D} <: AbstractMass{D}
    size::NTuple{D, Int}
end

UnitMass(size...) = UnitMass{length(size)}(size)

(m::UnitMass)(samples, ν) = m

Base.size(m::UnitMass) = m.size

Base.rand(m::UnitMass) = m |> size |> randn

Base.:\(::UnitMass, x) = x

Base.:+(m::UnitMass, ::UnitMass) = m

LinearAlgebra.logabsdet(::UnitMass) = 0
