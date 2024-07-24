struct State{T <: AbstractVecOrMat}
    position::T
    momentum::T
    gradient::T
    ll::Float64
end

State(θ, q₀) = begin ll, Δll = θ(q₀); State(q₀, mass(θ) * randn(q₀ |> size), Δll[1], ll) end

Base.copy(s::State) = State(q(s) |> copy, p(s) |> copy, a(s) |> copy, s.ll)

q(s::State) = s.position

p(s::State) = s.momentum

a(s::State) = s.gradient

ll(s::State) = s.ll

function leapfrog(θ, s₀, ϵ)
    s₁ = copy(s₀)
    p(s₁) .+= 0.5 * ϵ * a(s₁)
    q(s₁) .+= ϵ * (θ.mass\p(s₁))
    ll, Δll = θ(q(s₁))
    a(s₁) .= Δll[1]
    p(s₁) .+= 0.5 * ϵ * a(s₁)
    return State(q(s₁), p(s₁), a(s₁), ll)
end