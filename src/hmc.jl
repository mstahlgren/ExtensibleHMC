export HMC

struct HMC <: Sampler
    ϵ::Float64
    L::Int
end

function StatsBase.sample(ϕ::HMC, θ::Hamiltonian, q₀::AbstractVecOrMat)
    s₀ = s₁ = State(q₀, mass(θ) * randn(q₀ |> size), q₀ |> ∇(θ))
    for _ in 1:ϕ.L s₁ = leapfrog(θ, s₁, stepsize(ϕ)) end
    δH = sum(θ(s₀)) - sum(θ(s₁))
    accepted = δH > 0.0 || exp(δH) > rand()
    s = accepted ? s₁ : s₀
    return Sample(q(s), θ(s)[1], accepted, [s₀, s₁])
end