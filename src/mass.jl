abstract type Mass{T} end

struct UnitMass{T} <: Mass{T} end

struct DenseMass{T, S, M} <: Mass{T}
    McL::S
    M⁻¹::M
end

UnitMass(T = Float64) = UnitMass{T}()

function DenseMass(cov)
    c = cholesky(Hermitian(cholesky(cov) |> inv)).L
    DenseMass{eltype(cov), typeof(c), typeof(cov)}(c, cov)
end

Base.eltype(::Mass{T}) where T = T

Base.:*(m::Mass, x) = m.McL * x

Base.:*(::UnitMass, x) = x

Base.:\(m::Mass, x) = m.M⁻¹ * x

Base.:\(::UnitMass, x) = x
