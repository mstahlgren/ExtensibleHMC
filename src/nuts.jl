import LinearAlgebra: ⋅

export NUTS

struct NUTS <: Sampler
    ϵ::Float64
end

function StatsBase.sample(ϕ::NUTS, θ, q₀)
    s₀ = State(q₀, mass(θ) * randn(q₀ |> size), gradient(θ, q₀))
    u = rand() * exp(-energy(θ, s₀))
    s⁻, s⁺, s₁, j, n, h = s₀, s₀, s₀, 0, 1, true
    while h
        if rand(Bool) s⁻, _, s′, n′, h′ = buildtree(θ, s⁻, u, -1, j, stepsize(ϕ))
        else _, s⁺, s′, n′, h′ = buildtree(θ, s⁺, u, 1, j, stepsize(ϕ)) end
        if h′ && min(1, n′/n) > rand() s₁ = s′ end
        δq = q(s⁺) - q(s⁻)
        h = h′ && δq⋅p(s⁻) >= 0 && δq⋅p(s⁺) >= 0
        n += n′
        j += 1
    end
    return Sample(q(s₁), potential(θ, s₁), q₀ != q(s₁))
end

function buildtree(θ, s, u, v, j, ϵ)
    if iszero(j) return buildleaf(θ, s, u, v, ϵ) end
    s⁻, s⁺, s′, n′, h′ = buildtree(θ, s, u, v, j-1, ϵ)
    if !h′ return s⁻, s⁺, s′, n′, h′ end
    if v == -1 s⁻, _, s″, n″, h″ = buildtree(θ, s⁻, u, v, j-1, ϵ)
    else _, s⁺, s″, n″, h″ = buildtree(θ, s⁺, u, v, j-1, ϵ) end
    if n″ / (n′ + n″) > rand() s′ = s″ end
    δq = q(s⁺) - q(s⁻)
    h′ = h″ && δq⋅p(s⁻) >= 0 && δq⋅p(s⁺) >= 0
    n′ += n″
    return s⁻, s⁺, s′, n′, h′
end

function buildleaf(θ, s, u, v, ϵ)
    s′ = leapfrog(θ, s, v * ϵ)
    E′ = energy(θ, s′)
    n′ = Int(u <= exp(-E′))
    h′ = -E′ > log(u) - 1000
    return s′, s′, s′, n′, h′
end