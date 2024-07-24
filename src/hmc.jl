export HMC

struct HMC <: Sampler
    ϵ::Float64
    L::Int
end

function StatsBase.sample(ϕ::HMC, θ, q₀)
    s₀ = State(θ, q₀)
    for _ in 1:ϕ.L s₁ = leapfrog(θ, s₁, stepsize(ϕ)) end
    E₁ = energy(θ, s₁)
    δE = E₁ - energy(θ, s₀)
    d = -abs(δE) < - 1000
    accepted = !d && (δE < 0.0 || exp(-δE) > rand())
    s = accepted ? s₁ : s₀
    return Sample(q(s), potential(θ, s), accepted, d, false)
end