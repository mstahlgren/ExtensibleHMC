struct DiagMass{D} <: AbstractMass{D}
    M⁻¹::Array{Float64, D}
    N::Int
end

DiagMass(size...) = DiagMass{length(size)}(ones(size...), 1)

(m::DiagMass)(S, ν) = begin
    V = reshape(var(values(S), dims = 2, corrected = false), size(m))
    n₋ = round(ν * m.N)
    nₛ = Int(n₋) + length(S)
    typeof(m)((n₋ .* m.M⁻¹ .+ length(S) .* V) / nₛ, nₛ)
end

Base.size(m::DiagMass) = size(m.M⁻¹)

Base.rand(m::DiagMass) = sqrt.(1 ./ m.M⁻¹) .* randn(size(m))

Base.:\(m::DiagMass, x) = m.M⁻¹ .* x

Base.:+(x::T, y::T) where T <: DiagMass = T((x.M⁻¹ .+ y.M⁻¹) ./ 2, x.N + y.N)

LinearAlgebra.logabsdet(m::DiagMass) = -sum(log.(m.M⁻¹))