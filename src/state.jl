struct State{T <: AbstractVecOrMat}
    pos::T
    mom::T
    acc::T
end

State(θ, q₀) = State(q₀, mass(θ) * randn(q₀ |> size), gradient(θ, q₀)[2])

Base.copy(s::State) = State(q(s) |> copy, p(s) |> copy, a(s) |> copy)

q(s::State) = s.pos

p(s::State) = s.mom

a(s::State) = s.acc

function leapfrog(θ, s₀, ϵ)
    s₁ = copy(s₀)
    p(s₁) .+= 0.5 * ϵ * a(s₁)
    q(s₁) .+= ϵ * (θ.mass\p(s₁))
    ll, Δll = gradient(θ, q(s₁))
    a(s₁) .= Δll
    p(s₁) .+= 0.5 * ϵ * a(s₁)
    return ll, s₁
end