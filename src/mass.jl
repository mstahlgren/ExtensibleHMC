import LinearAlgebra: Diagonal

export UnitMass, DiagMass


abstract type Mass{T} end

Base.eltype(::Mass{T}) where T = T


struct UnitMass{T} <: Mass{T} end

UnitMass(T = Float64) = UnitMass{T}()


struct DiagMass{T, M} <: Mass{T}
    sqrmass::M
    invmass::M
end

function DiagMass(v)
    sqrmass, invmass = Diagonal(v .|> sqrt .|> inv), Diagonal(v)
    DiagMass{eltype(v), typeof(sqrmass)}(sqrmass, invmass)
end


Base.:*(::UnitMass, x) = x

Base.:*(m::DiagMass, x) = m.sqrmass * x

Base.:\(::UnitMass, x) = x

Base.:\(m::DiagMass, x) = m.invmass * x
