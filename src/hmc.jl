export HMC

struct HMC <: Sampler
    ϵ::Float64
    L::Int
end

function StatsBase.sample(ϕ::HMC, θ::Hamiltonian, q₀::AbstractVecOrMat)
    s₀ = s₁ = State(q₀, mass(θ) * randn(q₀ |> size), gradient(θ, q₀))
    for _ in 1:ϕ.L s₁ = leapfrog(θ, s₁, stepsize(ϕ)) end
    δE = energy(θ, s₁) - energy(θ, s₀)
    accepted = δE < 0.0 || exp(-δE) > rand()
    s = accepted ? s₁ : s₀
    return Sample(q(s), potential(θ, s), accepted)
end