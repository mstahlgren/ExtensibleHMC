import LinearAlgebra: ⋅

export NUTS

# CONSIDER: Max tree depth

struct NUTS <: Sampler
    ϵ::Float64
end

function uturn(s⁻, s⁺)
    δq = q(s⁺) - q(s⁻)
    return δq⋅p(s⁻) < 0 || δq⋅p(s⁺) < 0
end

function StatsBase.sample(ϕ::NUTS, θ, q₀; verbose = false)
    s₀ = State(θ, q₀)
    u = log(rand()) - energy(θ, s₀)
    s⁻, s⁺, s₁, j, n, t, d = s₀, s₀, s₀, 0, 1, false, false
    while j <= 8
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
    return Sample(q(s₁), s₁.ll, q₀ != q(s₁), d, j == 8)
end

function buildleft(ϕ, θ, s, u, j)
    if iszero(j) return buildleaf(θ, s, u, -stepsize(ϕ)) end
    s⁻, s⁺, s₁, n₁, t₁, d₁, ll₁ = buildleft(ϕ, θ, s, u, j-1)
    if d₁ return s⁻, s⁺, s₁, n₁, t₁, true, Inf end
    s⁻, _, s₂, n₂, t₂, d₂, ll₂ = buildleft(ϕ, θ, s⁻, u, j-1)
    n = n₂ / (n₁ + n₂) > rand()
    return s⁻, s⁺, n ? s₁ : s₂, n₁ + n₂, t₁ || t₂ || uturn(s⁻, s⁺), d₂
end

function buildright(ϕ, θ, s, u, j)
    if iszero(j) return buildleaf(θ, s, u, stepsize(ϕ)) end
    s⁻, s⁺, s₁, n₁, t₁, d₁ = buildright(ϕ, θ, s, u, j-1)
    if d₁ return s⁻, s⁺, s₁, n₁, t₁, true, Inf end
    _, s⁺, s₂, n₂, t₂, d₂ = buildright(ϕ, θ, s⁻, u, j-1)
    n = n₂ / (n₁ + n₂) > rand()
    return s⁻, s⁺, n ? s₁ : s₂, n₁ + n₂, t₁ || t₂ || uturn(s⁻, s⁺), d₂
end

function buildleaf(θ, s, u, ϵ)
    s′ = leapfrog(θ, s, ϵ)
    E = energy(θ, s′)
    n = Int(-E >= u)
    d = -E < u - 1000
    return s′, s′, s′, n, false, d
end