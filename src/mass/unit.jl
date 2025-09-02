struct UnitMass{D} <: AbstractMass{D}
    size::NTuple{D, Int}
end

UnitMass(size...) = UnitMass{length(size)}(size)

(m::UnitMass)(samples, ν) = m

Base.rand(::UnitMass, buffer) = randn!(pop!(buffer))

Base.:\(::UnitMass, x) = x
