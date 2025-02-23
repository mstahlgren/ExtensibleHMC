struct ColDiag <: Mass
    R::Int
    C::Int
    N::Int
    M⁻¹::Matrix{Float64}
end

ColDiag(R, C) = ColDiag(R, C, 1, ones(R, C))

function (m::ColDiag)(samples, ν = 0.0)
    expanded = reduce(hcat, vec(s.value) for s in samples)
    variance = reshape(var(expanded, dims = 2, corrected = false), R, C)
    N₀′ = round(ν * M.N); N₁ = Int(N₀′ + length(samples))
    ColDiag(M.R, M.C, N₁, (N₀′ .* m.M⁻¹ .+ length(samples) .* variance) / N₁)
end

Base.:+(x::ColDiag, y::ColDiag) = ColDiag(x.R, x.C, x.N + y.N, (x.M⁻¹ + y.M⁻¹) ./ 2)

Base.rand(m::ColDiag) = sqrt.(1.0 ./ m.M⁻¹) .* randn(m.R, m.C)

Base.:\(m::ColDiag, x) = m.M⁻¹ .* x

