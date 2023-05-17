export sample

function StatsBase.sample(q₀::Vector{Float64}, θ::Hamiltonian, ϕ::Integrator, σ::Sampler)
    p₀ = randn(size(q₀))
    path = ϕ(q₀, p₀, ∇(θ))
    return σ(path, θ)
end

function StatsBase.sample(n, q₀::Vector{Float64}, θ::Hamiltonian, ϕ::Integrator, σ::Sampler)
    samples = Vector{typeof(q₀)}(undef, n)
    for i in 1:n
        q₀ = sample(q₀, θ, ϕ, σ)
        samples[i] = q₀
    end
    return samples
end
# function plot_dlli(q, p, fun, tc, ϵ)
#     P₁ = dlli(q, p, fun, tc, ϵ)
#     q̂, p̂ = last(P₁)
#     P₂ = dlli(q̂, -p̂, fun, tc, ϵ)
#     d₁ = reduce(hcat, [vcat(x...) for x in P₁])'
#     d₂ = reduce(hcat, [vcat(x...) for x in P₂])'
#     plot(d₁[:,1], d₁[:,2])
#     scatter!(d₂[:,1], d₂[:,2])
# end
