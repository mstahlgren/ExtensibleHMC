# Original No U-turn sampler using slice sampling.

# NOTE: This Sampler currently does not return a proper acceptrate

struct NUTS <: Sampler
    ϵ::Float64
    max_depth::Int
    max_ΔE::Float64
end

NUTS(ϵ) = NUTS(ϵ, 12, 1000)

function uturn(s⁻, s⁺)
    δq = q(s⁺) - q(s⁻)
    return δq⋅p(s⁻) < 0 || δq⋅p(s⁺) < 0
end

function StatsBase.sample(ϕ::NUTS, θ, q₀; verbose = false)
    s₀ = State(θ, q₀)
    u = log(rand()) - energy(s₀)
    s⁻, s⁺, s₁, j, n, t, d = s₀, s₀, s₀, 0, 1, false, false
    while j < ϕ.max_depth
        if verbose print("New branch :: j ", j) end
        if rand(Bool) s⁻, _, s′, n′, t, d = buildleft(ϕ, θ, s⁻, u, j)
        else _, s⁺, s′, n′, t, d = buildright(ϕ, θ, s⁺, u, j) end
        if verbose println(" :: internal uturn ", t, " :: diverged ", d, " :: uturn ", uturn(s⁻, s⁺)) end
        if t || d break end
        if n′/n > rand() s₁ = s′ end
        if uturn(s⁻, s⁺) break end
        n += n′
        j += 1
    end
    return Sample(q(s₁), s₁.ll, 2^j, 0.0, d, j == ϕ.max_depth)
end

function buildleft(ϕ::NUTS, θ, s, u, j)
    if iszero(j) return buildleaf(ϕ, θ, s, u, -ϕ.ϵ) end
    s⁻, s⁺, s₁, n₁, t₁, d₁ = buildleft(ϕ, θ, s, u, j-1)
    if d₁ return s⁻, s⁺, s₁, n₁, t₁, true end
    s⁻, _, s₂, n₂, t₂, d₂ = buildleft(ϕ, θ, s⁻, u, j-1)
    n = n₂ / (n₁ + n₂) > rand()
    return s⁻, s⁺, n ? s₁ : s₂, n₁ + n₂, t₁ || t₂ || uturn(s⁻, s⁺), d₂
end

function buildright(ϕ::NUTS, θ, s, u, j)
    if iszero(j) return buildleaf(ϕ, θ, s, u, ϕ.ϵ) end
    s⁻, s⁺, s₁, n₁, t₁, d₁ = buildright(ϕ, θ, s, u, j-1)
    if d₁ return s⁻, s⁺, s₁, n₁, t₁, true end
    _, s⁺, s₂, n₂, t₂, d₂ = buildright(ϕ, θ, s⁻, u, j-1)
    n = n₂ / (n₁ + n₂) > rand()
    return s⁻, s⁺, n ? s₁ : s₂, n₁ + n₂, t₁ || t₂ || uturn(s⁻, s⁺), d₂
end

function buildleaf(ϕ::NUTS, θ, s, u, ϵ)
    s′ = leapfrog(θ, s, ϵ)
    E = energy(s′)
    n = Int(-E >= u)
    d = -E < u - ϕ.max_ΔE
    return s′, s′, s′, n, false, d
end