struct ColDiag{S, M} <: Mass
    R::Int
    C::Int
    N::Int
    M⁻¹::M
end

ColDiag(R, C) = ColDiag(R, C, 1, ones(R, C))

function (m::ColDiag)(samples, ν = 0.0) where T <: AbstractMatrix
    expanded = reduce(hcat, vec(s.value) for s in samples)
    variance = reshape(var(expanded, dims = 2, corrected = false), R, C)
    N₀′ = round(ν * M.N); N₁ = Int(N₀′ + length(samples))
    ColDiag(M.R, M.C, N₁, (N₀′ .* m.M⁻¹ .+ length(samples) .* variance) / N₁)
end

Base.:+(x::ColDiag, y::ColDiag) = ColDiag(x.R, x.C, x.N + y.N, (x.M⁻¹ + y.M⁻¹) ./ 2)

Base.rand(m::ColDiag) = sqrt.(1.0 ./ m.M⁻¹) .* randn(m.R, m.C)

Base.:\(::ColDiag, x) = M⁻¹ .* x

