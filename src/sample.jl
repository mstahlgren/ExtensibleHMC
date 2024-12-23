abstract type Sampler end

struct Sample{T <: AbstractVecOrMat}
    value::T
    ll::Float64
    nsteps::Int
    acceptrate::Float64
    diverged::Bool
    maxdepth::Bool
end

const Samples{T} = Vector{Sample{T}}

function StatsBase.sample(ϕ::Sampler, θ::Hamiltonian, q::AbstractVecOrMat, n::Int; verbose = false)
    samples = Vector{Sample{typeof(q)}}(undef, n)
    for i in 1:n
        if verbose println("Sample ", i) end
        s = sample(ϕ, θ, q)
        samples[i] = s
        q = s.value
        if verbose 
            print("Completed :: LL ", round(s.ll, digits = 4))
            print(" :: in $s.nsteps steps")
            if !s.accepted print(" :: Rejected") end
            if s.diverged print(" :: Diverged") end
            println("")
        end
    end
    return samples
end

Base.show(io::IO, ::MIME"text/plain", x::Samples) = summary(x)
