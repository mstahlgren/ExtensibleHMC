struct UnitMass <: AbstractMass
    size::Int
end

(m::UnitMass)(samples, ν) = m

Base.rand(::UnitMass, buffer) = randn!(pop!(buffer))

Base.:\(::UnitMass, x) = x
