import StatsBase: StatsBase, sample

export Sample, sample

# CONSIDER: Kinetic energy? Potential energy?
struct Sample
    q::Vector{Float64}
    accepted::Bool
    path::Vector{State}
end

function metropolis(path, θ)
    s₀ = first(path)
    s₁ = last(path)
    ΔH = sum(θ(s₀)) - sum(θ(s₁))
    return ΔH > 0.0 || exp(ΔH) > rand() ? q(s₁) : q(s₀)
end

function StatsBase.sample(q₀::Vector{Float64}, θ::Hamiltonian, ϕ::Integrator)
    p₀ = θ.mass.L * randn(size(q₀)) # This isn't very clear
    path = ϕ(State(q₀, p₀), θ)
    q₁ = metropolis(path, θ)
    return Sample(q₁, q₀ != q₁, path)
end

function StatsBase.sample(n, q::Vector{Float64}, θ::Hamiltonian, ϕ::Integrator)
    samples = Vector{Sample}(undef, n)
    for i in 1:n
        s = sample(q, θ, ϕ)
        samples[i] = s
        q = s.q
    end
    return samples
end

Base.Matrix(samples::Vector{Sample}) = reduce(hcat, [s.q for s in samples])'