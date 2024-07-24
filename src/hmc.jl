export HMC

struct HMC <: Sampler
    ϵ::Float64
    L::Int
end

function StatsBase.sample(ϕ::HMC, θ, q₀; verbose = false)
    s₀, s₁ = State(θ, q₀), leapfrog(θ, s₀, stepsize(ϕ))
    for _ in 2:ϕ.L s₁ = leapfrog(θ, s₁, stepsize(ϕ)) end
    E₁ = energy(θ, s₁)
    δE = E₁ - energy(θ, s₀)
    d = -abs(δE) < - 1000
    accepted = !d && (δE < 0.0 || exp(-δE) > rand())
    s = accepted ? s₁ : s₀
    return Sample(q(s), potential(θ, s), accepted, d, false)
end