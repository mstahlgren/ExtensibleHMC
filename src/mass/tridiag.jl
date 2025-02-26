struct TriDiagMass <: AbstractMass
    R::Int
    C::Int
    N::Int
    L::Bidiagonal{Float64, Vector{Float64}}
    M⁻¹::SymTridiagonal{Float64, Vector{Float64}}
end

function TriDiagMass(R, C)
    D, S = ones(R * C), zeros(R * C - 1)
    TriDiagMass(R, C, 1, Bidiagonal(copy(D), copy(S), :L), SymTridiagonal(D, S))
end

function factorise(m::SymTridiagonal, buff = Bidiagonal(similar(m.dv), similar(m.ev), :L))
    buff.dv[1] = sqrt(m.dv[1])
    for i in 2:length(buff.dv)
        buff.ev[i-1] = m.ev[i-1] / buff.dv[i-1]
        buff.dv[i] = sqrt(m.dv[i] - buff.ev[i-1]^2)
    end
    buff
end

function (m::TriDiagMass)(samples, ν = 0.0)
    W = reduce(hcat, vec(s.value') for s in samples)
    W .-= mean(W, dims = 2)
    N = length(m)
    dv⁻¹ = sum(W.^2, dims = 2) |> vec
    ev⁻¹ = sum(view(W, 1:N-1, :) .* view(W, 2:N, :), dims = 2) |> vec
    N₀ = round(ν * m.N)
    Nₛ = Int(N₀ + length(samples))
    M⁻¹ = SymTridiagonal((N₀ .* m.M⁻¹.dv .+ dv⁻¹) ./ Nₛ, (N₀ .* m.M⁻¹.ev .+ ev⁻¹) ./ Nₛ)
    McL = factorise(M⁻¹)
    TriDiagMass(m.R, m.C, Nₛ, McL, M⁻¹)
end

function Base.:+(x::TriDiagMass, y::TriDiagMass)
    M⁻¹ = SymTridiagonal((x.M⁻¹.dv .+ y.M⁻¹.dv) ./ 2, (x.M⁻¹.ev .+ y.M⁻¹.ev) ./ 2)
    McL = factorise(M⁻¹)
    TriDiagMass(x.R, x.C, x.N + y.N, McL, M⁻¹)
end

# Verify with https://math.stackexchange.com/questions/1845132/sample-points-from-a-multivariate-normal-distribution-using-only-the-precision-m
# This is likely wrong as L should be used with \
Base.rand(m::TriDiagMass) = reshape(m.L * randn(length(m)), m.C, m.R)' |> collect

Base.:\(m::TriDiagMass, x) = reshape(m.M⁻¹ * vec(x'), m.C, m.R)' |> collect

LinearAlgebra.logabsdet(m::TriDiagMass) = -2 * sum(log.(m.L.dv))