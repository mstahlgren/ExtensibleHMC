struct State{T <: AbstractVecOrMat}
    pos::T
    mom::T
    acc::T
end

Base.copy(s::State) = State(q(s) |> copy, p(s) |> copy, a(s) |> copy)

q(s::State) = s.pos

p(s::State) = s.mom

a(s::State) = s.acc

function leapfrog(θ, s₀, ϵ)
    s₁ = copy(s₀)
    p(s₁) .+= 0.5 * ϵ * a(s₁)
    q(s₁) .+= ϵ * (θ.mass\p(s₁))
    a(s₁) .= gradient(θ, q(s₁))
    p(s₁) .+= 0.5 * ϵ * a(s₁)
    return s₁
end