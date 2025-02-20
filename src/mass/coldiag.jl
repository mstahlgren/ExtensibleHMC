struct ColDiag{S, M} <: Mass
    R::Int
    C::Int
    N::Int
    M⁻¹::M
end

ColDiag(R, C) = ColDiag(R, C, 1, ones(R, C))

function ColDiag(S::Samples, M = ColDiag(size(S[1].value)...), ν = 0.0)
    expanded = reduce(hcat, vec(s.value) for s in S)
    variance, frac = var(expanded, dims = 2, corrected = false), round(ν * M.N)
    ColDiag(M.R, M.C, Int(frac + length(S)), (frac / M.N) .* M⁻¹ .+ variance)
end

Base.:*(::ColDiag, x) = sqrt.(1 ./ M⁻¹)' .* x

Base.:\(::ColDiag, x) = M⁻¹ .* x