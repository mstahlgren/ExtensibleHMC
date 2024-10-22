import LinearAlgebra: Diagonal

export UnitMass, DiagMass


abstract type Mass{T, N} end

Base.eltype(::Mass{T,N}) where {T,N} = T

Base.length(::Mass{T,N}) where {T,N} = N


struct UnitMass{T, N} <: Mass{T, N} end

UnitMass(N::Int, T = Float64) = UnitMass{T, N}()


struct DiagMass{T, N, M} <: Mass{T, N}
    sqrmass::M
    invmass::M
end

function DiagMass(v)
    sqrmass, invmass = Diagonal(v .|> sqrt .|> inv), Diagonal(v)
    DiagMass{eltype(v), length(v), typeof(sqrmass)}(sqrmass, invmass)
end


Base.:*(::UnitMass, x) = x

Base.:*(m::DiagMass, x) = m.sqrmass * x

Base.:\(::UnitMass, x) = x

Base.:\(m::DiagMass, x) = m.invmass * x
