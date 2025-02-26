struct DiagPlus <: AbstractMass
    R::Int
    C::Int
    N::Int
    L::Bidiagonal{Float64, Vector{Float64}}
    M⁻¹::Vector{Float64}
    ρ⁻¹::Float64
end

function DiagPlus(R, C)
    D, S = ones(R * C), zeros(R * C - 1)
    DiagPlus(R, C, 1, Bidiagonal(copy(D), S, :L), D, 0.0)
end

function factorise(dv::Vector, ρ)
    buff = Bidiagonal(similar(dv), similar(dv, length(dv)-1), :L)
    ρ², buff.dv[1] = ρ^2, sqrt(dv[1])
    for i in 2:length(dv)
        buff.ev[i-1] = ρ / buff.dv[i-1]
        buff.dv[i] = sqrt(dv[i] - ρ²)
    end
    buff
end

function (m::DiagPlus)(samples, ν = 0.0)
    W = reduce(hcat, vec(s.value') for s in samples)
    W .-= mean(W, dims = 2)
    N = length(m)
    dv⁻¹ = sum(W.^2, dims = 2) |> vec
    ρ⁻¹ = mean(view(W, 1:N-1, :) .* view(W, 2:N, :))
    N₀ = round(ν * m.N)
    Nₛ = Int(N₀ + length(samples))
    M⁻¹ = (N₀ .* m.M⁻¹ .+ dv⁻¹) ./ Nₛ
    ρₛ⁻¹ = (N₀ .* m.ρ⁻¹ .+ length(samples) * ρ⁻¹) ./ Nₛ
    McL = factorise(M⁻¹, ρₛ⁻¹)
    DiagPlus(m.R, m.C, Nₛ, McL, M⁻¹, ρₛ⁻¹)
end

function Base.:+(x::DiagPlus, y::DiagPlus)
    M⁻¹, ρ = (x.M⁻¹ .+ y.M⁻¹) ./ 2, x.ρ⁻¹ + y.ρ⁻¹
    McL = factorise(M⁻¹, ρ)
    DiagPlus(x.R, x.C, x.N + y.N, McL, M⁻¹, ρ)
end

# This should be wrong
Base.rand(m::DiagPlus) = reshape(m.L * randn(length(m)), m.C, m.R)' |> collect

Base.:\(m::DiagPlus, x) = begin
    d = vec(x')
    dd = m.M⁻¹ .* d
    view(dd, 1:lastindex(d)-1) .+= m.ρ⁻¹ * view(d, 2:lastindex(d))
    view(dd, 2:lastindex(d)) .+= m.ρ⁻¹ * view(d, 1:lastindex(d)-1)
    reshape(dd, m.C, m.R)' |> collect
end

# This is likely wrong, shouldnt be negative if using L
LinearAlgebra.logabsdet(m::DiagPlus) = -2 * sum(log.(m.L.dv))