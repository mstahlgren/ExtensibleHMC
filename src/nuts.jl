import LinearAlgebra: ⋅

export NUTS

# CONSIDER: Max tree depth
# CONSIDER: Use pullback instead of gradient and make use of ll. Measure time

struct NUTS <: Sampler
    ϵ::Float64
end

function uturn(s⁻, s⁺)
    δq = q(s⁺) - q(s⁻)
    return δq⋅p(s⁻) < 0 || δq⋅p(s⁺) < 0
end

function StatsBase.sample(ϕ::NUTS, θ, q₀)
    s₀ = State(q₀, mass(θ) * randn(q₀ |> size), gradient(θ, q₀)) # This suggests Chol of Mass
    u = rand() * exp(-energy(θ, s₀))
    s⁻, s⁺, s₁, j, n, t, d = s₀, s₀, s₀, 0, 1, false, false
    while true
        if rand(Bool) s⁻, _, s′, n′, t, d = buildleft(ϕ, θ, s⁻, u, j)
        else _, s⁺, s′, n′, t, d = buildright(ϕ, θ, s⁺, u, j) end
        if t || d break end
        if n′/n > rand() s₁ = s′ end
        if uturn(s⁻, s⁺) break end
        n += n′
        j += 1
    end
    return Sample(q(s₁), potential(θ, s₁), q₀ != q(s₁), d)
end

function buildleft(ϕ, θ, s, u, j)
    if iszero(j) return buildleaf(θ, s, u, -stepsize(ϕ)) end
    s⁻, s⁺, s₁, n₁, t₁, d₁ = buildleft(ϕ, θ, s, u, j-1)
    if d₁ return s⁻, s⁺, s₁, n₁, t₁, true end
    s⁻, _, s₂, n₂, t₂, d₂ = buildleft(ϕ, θ, s⁻, u, j-1)
    return s⁻, s⁺, n₂ / (n₁ + n₂) > rand() ? s₁ : s₂, n₁ + n₂, t₁ || t₂ || uturn(s⁻, s⁺), d₂
end

function buildright(ϕ, θ, s, u, j)
    if iszero(j) return buildleaf(θ, s, u, stepsize(ϕ)) end
    s⁻, s⁺, s₁, n₁, t₁, d₁ = buildright(ϕ, θ, s, u, j-1)
    if d₁ return s⁻, s⁺, s₁, n₁, t₁, true end
    _, s⁺, s₂, n₂, t₂, d₂ = buildright(ϕ, θ, s⁻, u, j-1)
    return s⁻, s⁺, n₂ / (n₁ + n₂) > rand() ? s₁ : s₂, n₁ + n₂, t₁ || t₂ || uturn(s⁻, s⁺), d₂
end

function buildleaf(θ, s, u, ϵ)
    s′ = leapfrog(θ, s, ϵ)
    E = energy(θ, s′)
    n = Int(exp(-E) >= u)
    d = -E < log(u) - 1000
    return s′, s′, s′, n, false, d
end