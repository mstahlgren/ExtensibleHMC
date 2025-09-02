struct State{T <: AbstractArray}
    position::T
    momentum::T
    gradient::T
    ll::Float64
    ke::Float64
end

function State(θ, q₀, buffer)
    println("hello")
    display(111111111)
    display(θ)
    p, ll, Δll = refresh(θ, buffer), θ(q₀, pop!(buffer))...
    return State(copy!(pop!(buffer), q₀), p, Δll, ll, kinetic(θ, p))
end

q(s::State) = s.position

p(s::State) = s.momentum

a(s::State) = s.gradient

ll(s::State) = s.ll

ke(s::State) = s.ke

energy(s::State) = s.ke - s.ll

function leapfrog(θ, s₀, ϵ, buffer)
    q₁, p₁ = pop!(buffer), pop!(buffer)
    p₁ .= p(s₀) .+ 0.5 .* ϵ .* a(s₀)
    q₁ .= q(s₀) .+ ϵ .* v(θ, p₁)
    ll, a₁ = θ(q₁, a(s₀))
    p₁ .+= 0.5 .* ϵ .* a₁
    return State(q₁, p₁, a₁, ll, kinetic(θ, p₁))
end

function Base.show(io::IO, s::State)
    print("State(")
    print("U: ", round(-s.ll, digits = 3), " ")
    print("K: ", round(s.ke, digits = 3), " ")
    print("E: ", round(energy(s), digits = 3))
    print(")")
end