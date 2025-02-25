struct DiagMass{N} <: AbstractMass{N}
    M⁻¹::Vector{Float64}
    nₛ::Int
end

DiagMass(N) = DiagMass{N}(ones(N), 1)

function (m::DiagMass)(S, ν = 0.0)
    V = var(values(S), dims = 2, corrected = false)
    n₁ = round(ν * m.N)
    nₛ = Int(N₀′) + length(S)
    typeof(m)((n₁ .* m.M⁻¹ .+ length(S) .* V) / nₛ, nₛ)
end

Base.:+(x::T, y::T) where T <: DiagMass = T((x.M⁻¹ .+ y.M⁻¹) ./ 2, x.N + y.N)

Base.rand(m::DiagMass) = sqrt.(1 ./ m.M⁻¹) .* randn(length(m))

Base.:\(m::DiagMass, x) = m.M⁻¹ .* x

LinearAlgebra.logabsdet(m::DiagMass) = -sum(log.(m.M⁻¹))