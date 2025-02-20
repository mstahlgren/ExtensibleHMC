struct DenseMass{S, M} <: Mass
    R::Int
    C::Int
    N::Int
    McL::S
    M⁻¹::M
end

DenseMass(R, C) = DenseMass(R, C, 1, LowerTriangular(Matrix(1.0I, R, R)), Matrix(1.0I, R, R))

function (m::DenseMass)(S::Samples{T}, ν = 0.0) where T <: AbstractMatrix
    expanded = reduce(vcat, vec(s.value)' for s in S)
    covariance = sum(cov(view(M, :, i:i+N-1), corrected = false) for i in 1:N:size(expanded, 2))
    N₀′ = round(ν * m.N)
    N₁ = Int(N₀′ + length(S))
    M⁻¹₁ = (N₀′ .* m.M⁻¹ .+ covariance) ./ N₁
    McL₁ = cholesky(Hermitian(cholesky(M⁻¹₁) |> inv)).L
    DenseMass(M.R, M.C, N₁, McL₁, M⁻¹₁)
end

Base.:+(x::DenseMass, y::DenseMass) = begin
    M⁻¹ = (x.M⁻¹ + y.M⁻¹) ./ 2
    McL₁ = cholesky(Hermitian(cholesky(M⁻¹) |> inv)).L
    DenseMass(x.R, x.C, x.N + y.N, McL₁, M⁻¹)
end

Base.rand(m::DenseMass, x) = m.McL * randn(m.R, m.C)

Base.:\(m::Mass, x) = m.M⁻¹ * x
