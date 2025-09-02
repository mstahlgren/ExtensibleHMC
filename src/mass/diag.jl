struct DiagMass{D} <: AbstractMass{D}
    M⁻¹::Array{Float64, D}
    N::Int
end

DiagMass(size...) = DiagMass{length(size)}(ones(size...), 1)

(m::DiagMass)(S, ν) = begin
    V = reshape(var(values(S), dims = 2, corrected = false), size(m))
    n₋ = round(ν * m.N); nₛ = Int(n₋) + length(S); η = n₋ / nₛ
    typeof(m)(η .* m.M⁻¹ .+ V, nₛ)
end

Base.rand(m::DiagMass, buffer) = begin
    x = randn!(pop!(buffer))
    x .*= sqrt.(1 ./ m.M⁻¹)
    return x
end

Base.:\(m::DiagMass, x) = m.M⁻¹ .* x
