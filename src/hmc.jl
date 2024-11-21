# Plain Hybrid Monte Carlo sampler using a fixed steps and metropolis endpoint sampling.

struct HMC <: Sampler
    ϵ::Float64
    L::Int
    max_ΔE::Float64
end

HMC(ϵ, L) = HMC(ϵ, L, 1000)

function StatsBase.sample(ϕ::HMC, θ, q₀; verbose = false)
    s₀ = s₁ = State(θ, q₀)
    for _ in 1:ϕ.L s₁ = leapfrog(θ, s₁, ϕ.ϵ) end
    ΔE = energy(s₁) - energy(s₀)
    d = abs(ΔE) > ϕ.max_ΔE
    p = min(1.0, exp(-ΔE))
    accepted = !d && (ΔE < 0.0 || rand() < p)
    s = accepted ? s₁ : s₀
    return Sample(q(s), ll(s), ϕ.L, p, d, false)
end