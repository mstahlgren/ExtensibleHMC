import Statistics: var

struct DiagMass <: AbstractMass
    M⁻¹::Vector{Float64}
    N::Int
end

DiagMass(m) = DiagMass(ones(m), 1)

(m::DiagMass)(S, ν) = begin
    V = var(values(S), dims = 2, corrected = false)
    n₋ = round(ν * m.N); nₛ = Int(n₋) + length(S); η = n₋ / nₛ
    DiagMass(η .* m.M⁻¹ .+ V, nₛ)
end

Base.rand(m::DiagMass, buffer) = begin
    x = randn!(pop!(buffer))
    x .*= sqrt.(1 ./ m.M⁻¹)
    return x
end

Base.:\(m::DiagMass, x) = m.M⁻¹ .* x
