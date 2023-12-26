import StatsBase: StatsBase, sample

export Sample, sample

struct Sample{T <: AbstractVecOrMat}
    state::T
    ll::Float64
    accepted::Bool
    path::Vector{State}
end

function metropolis(s₀, s₁, θ)
    println(s₀)
    println(s₁)
    ΔH = sum(θ(s₀)) - sum(θ(s₁))
    return ΔH > 0.0 || exp(ΔH) > rand() ? q(s₁) : q(s₀)
end

function StatsBase.sample(q₀::AbstractVecOrMat, θ::Hamiltonian, ϕ::Integrator)
    p₀ = θ.mass.L * randn(q₀ |> size)
    path = ϕ(State(q₀, p₀), θ)
    q₁ = metropolis(first(path), rand(path), θ)
    return Sample(q₁, θ(q₁)[1], q₀ != q₁, path)
end

function StatsBase.sample(q::T, θ::Hamiltonian, ϕ::Integrator, n) where T <: AbstractVecOrMat
    samples = Vector{Sample{T}}(undef, n)
    for i in 1:n
        s = sample(q, θ, ϕ)
        samples[i] = s
        q = s.state
    end
    return samples
end
