struct State{T <: AbstractVecOrMat}
    position::T
    momentum::T
    gradient::T
    ll::Float64
    ke::Float64
end

function State(θ, q₀)
    ll, Δll = θ(q₀)
    p = mass(θ) * randn(q₀ |> size)
    return State(copy(q₀), p, Δll[1], ll, kinetic(θ, p))
end

q(s::State) = s.position

p(s::State) = s.momentum

a(s::State) = s.gradient

ll(s::State) = s.ll

energy(s::State) = s.ke - s.ll

function leapfrog(θ, s₀, ϵ)
    p₁ = p(s₀) .+ 0.5 .* ϵ .* a(s₀)
    q₁ = q(s₀) .+ ϵ .* v(θ, p₁)
    ll, Δll = θ(q₁)
    a₁ = Δll[1]
    p₁ .+= 0.5 .* ϵ .* a₁
    return State(q₁, p₁, a₁, ll, kinetic(θ, p₁))
end

function leapdiff(θ, s₀, ϵ)
    s₁ = leapfrog(θ, s₀, ϵ)
    return s₁, energy(s₁) - energy(s₀)
end

function Base.show(io::IO, s::State)
    print("State(")
    print("U: ", round(-s.ll, digits = 3), " ")
    print("K: ", round(s.ke, digits = 3), " ")
    print("E: ", round(energy(s), digits = 3))
    print(")")
end