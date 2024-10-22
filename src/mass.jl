import LinearAlgebra: Diagonal

export UnitMass, DiagMass

abstract type Mass end

struct UnitMass <: Mass end

struct DiagMass{T} <: Mass
    sqrmass::T
    invmass::T
end

function DiagMass(w::Welford)
    sqrmass = w |> var .|> sqrt .|> inv |> Diagonal
    invmass = w |> var |> Diagonal
    NewMass{typeof(sqrmass)}(sqrmass, invmass)
end

Base.:*(::UnitMass, x) = x

Base.:*(m::DiagMass, x) = m.sqrmass * x

Base.:\(::UnitMass, x) = x

Base.:\(m::DiagMass, x) = m.invmass * x
