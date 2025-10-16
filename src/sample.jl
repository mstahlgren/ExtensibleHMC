abstract type Sampler end

struct Sample{T <: AbstractArray}
    value::T
    ll::Float64
    nsteps::Int
    acceptrate::Float64
    diverged::Bool
    maxdepth::Bool
end

const Samples{T} = Vector{Sample{T}}

function sample(fun::Function, q::AbstractArray, n::Int, step = 0.05)
    return sample(MNUTS(step), Hamiltonian(fun, UnitMass(length(q))), q, n)
end

function sample(ϕ::Sampler, θ::Hamiltonian, q::T, n::Int; verbose = false) where T <: AbstractArray
    buffer = Buffer(ϕ, length(q))
    samples = Samples{T}(undef, n)
    for i in 1:n
        s = sample(ϕ, θ, q, buffer)
        q, samples[i] = s.value, s
        reset!(buffer)
        if !verbose continue end
        print("Sample $i Completed :: LL ", round(s.ll, digits = 4), " :: in $(s.nsteps) steps")
        if s.diverged println(" :: Diverged") 
        elseif !s.accepted println(" :: Rejected") 
        else println("") end
    end
    return samples
end

function adapt(ϕ::Sampler, θ::Hamiltonian, q, epochs, n, β)
    for e in 1:epochs
        S = sample(ϕ, θ, q, n)
        α, q = acceptrate(S), last(S).value
        ν = if α > 0.9 1.0/0.90 elseif α < 0.7 0.90 else 1.0 end
        ϕ = MNUTS(ϕ.ϵ * ν, ϕ.max_depth)
        θ = α > 0.4 ? θ(θ.mass(S, β)) : θ
        println("Epoch $e completed :: accept $α :: step $(ϕ.ϵ)")
    end
    return ϕ, θ, q
end

Base.show(io::IO, ::MIME"text/plain", x::Samples) = summary(x)
