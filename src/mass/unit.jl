struct UnitMass <: AbstractMass
    R::Int
    C::Int
end

(m::UnitMass)(samples, ν = 0.0) = UnitMass(m.R, m.C)

Base.rand(m::UnitMass) = randn(m.R, m.C)

Base.:\(::UnitMass, x) = x

Base.:+(m::UnitMass, y::UnitMass) = UnitMass(m.R, m.C)
