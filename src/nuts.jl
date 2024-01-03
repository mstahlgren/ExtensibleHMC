import LinearAlgebra: ⋅

export State2, p, q, a
export NUTS, stepsize
export buildtree, buildleaf,leapfrog

struct State2{T <: AbstractVecOrMat}
    q::T
    p::T
    a::T
end

q(s::State2) = s.q

p(s::State2) = s.p

a(s::State2) = s.a

Base.copy(s::State2) = State2(s.q |> copy, s.p |> copy, s.a |> copy)

struct NUTS <: Integrator
    ϵ::Float64
end

stepsize(x::Integrator) = x.ϵ

function sample(ϕ::NUTS, θ, q₀)
    s₀ = State2(q₀, mass(θ) * rand(length(q₀)), q₀ |> ∇(θ))
    u = rand() * exp(reduce(-, θ(s₀)))
    s⁻, s⁺, s₁, j, n, h = s₀, s₀, s₀, 0, 1, true
    while h
        v = rand([1, -1])
        if v == -1 s⁻, _, s′, n′, h′ = buildtree(θ, s⁻, u, v, j, stepsize(ϕ))
        else _, s⁺, s′, n′, h′ = buildtree(θ, s⁺, u, v, j, stepsize(ϕ)) end
        if h′ && min(1, n′/n) > rand() s₁ = s′ end
        δq = q(s⁺) - q(s⁻)
        h = h′ && δq⋅p(s⁻) >= 0 && δq⋅p(s⁺) >= 0
        n += n′
        j += 1
    end
    return Sample(q(s₁), θ(s₁)[1], q₀ != q(s₁), [State(q₀, q₀)])
end

function leapfrog(θ, s₀, ϵ)
    s₁ = copy(s₀)
    p(s₁) .-= 0.5 * ϵ * a(s₁)
    q(s₁) .+= ϵ * (θ.mass\p(s₁))
    a(s₁) .= q(s₁) |> ∇(θ)
    p(s₁) .-= 0.5 * ϵ * a(s₁)
    return s₁
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
    δE = reduce(-, θ(s′))
    n′ = Int(u <= exp(δE))
    h′ = δE > log(u) - 1000
    return s′, s′, s′, n′, h′
end