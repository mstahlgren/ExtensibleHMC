import LinearAlgebra: I, Diagonal, LowerTriangular, Cholesky

abstract type Mass end

struct UnitMass end

struct DiagMass
    value::Diagonal{Float64}
end

struct DenseMass
    value::LowerTruangular{Float64}
end

Base.:*(::UnitMass, x) = x

Base.:*(m::DiagMass, x) = m.value * x

Base.:*(m::DenseMass, x) = m.value * x

Base.:\(::UnitMass, x) = x

Base.:\(m::DiagMass, x) = m.value\x

Base.:\(m::DenseMass, x) = Cholesky(m.value)\x