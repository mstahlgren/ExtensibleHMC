import LinearAlgebra: I, Diagonal, LowerTriangular, Cholesky

export UnitMass, DiagMass, DenseMass

abstract type Mass end

struct UnitMass <: Mass end

struct DiagMass <: Mass
    value::Diagonal{Float64}
end

DiagMass(x::Vector{T}) where T = x |> Diagonal |> DiagMass

struct DenseMass <: Mass
    value::LowerTriangular{Float64}
end

Base.:*(::UnitMass, x) = x

Base.:*(m::DiagMass, x) = m.value * x

Base.:*(m::DenseMass, x) = m.value * x

Base.:\(::UnitMass, x) = x

Base.:\(m::DiagMass, x) = (m.value.^2)\x

Base.:\(m::DenseMass, x) = Cholesky(m.value)\x