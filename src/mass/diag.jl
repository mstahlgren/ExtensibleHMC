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

Base.size(m::DiagMass) = size(m.M⁻¹)

Base.:\(m::DiagMass, x) = m.M⁻¹ .* x

Base.:+(x::T, y::T) where T <: DiagMass = T((x.M⁻¹ .+ y.M⁻¹) ./ 2, x.N + y.N)

LinearAlgebra.logabsdet(m::DiagMass) = -sum(log.(m.M⁻¹))