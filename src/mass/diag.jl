struct DiagMass <: AbstractMass
    R::Int
    C::Int
    N::Int
    M⁻¹::Matrix{Float64}
end

DiagMass(R, C) = DiagMass(R, C, 1, ones(R, C))

function (m::DiagMass)(samples, ν = 0.0)
    expanded = reduce(hcat, vec(s.value) for s in samples)
    variance = reshape(var(expanded, dims = 2, corrected = false), m.R, m.C)
    N₀′ = round(ν * m.N); N₁ = Int(N₀′ + length(samples))
    DiagMass(m.R, m.C, N₁, (N₀′ .* m.M⁻¹ .+ length(samples) .* variance) / N₁)
end

Base.:+(x::DiagMass, y::DiagMass) = ColDiag(x.R, x.C, x.N + y.N, (x.M⁻¹ + y.M⁻¹) ./ 2)

Base.rand(m::DiagMass) = sqrt.(1.0 ./ m.M⁻¹) .* randn(m.R, m.C)

Base.:\(m::DiagMass, x) = m.M⁻¹ .* x

