struct UnitMass <: Mass
    R::Int
    C::Int
end

(m::UnitMass)(samples, ν = 0.0) = UnitMass(M.R, M.C)

Base.rand(m::UnitMass) = randn(m.R, M.C)

Base.:\(::UnitMass, x) = x

Base.:+(x::UnitMass, y::UnitMass) = UnitMass(x.R, x.C)
