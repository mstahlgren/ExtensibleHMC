import StatsBase: StatsBase, sample

export Sample, sample

abstract type Sampler end

stepsize(x::Sampler) = x.ϵ

struct Sample{T <: AbstractVecOrMat}
    value::T
    ll::Float64
    accepted::Bool
    diverged::Bool
end

function StatsBase.sample(ϕ::Sampler, θ::Hamiltonian, q::T, n, verbose = false) where T <: AbstractVecOrMat
    samples = Vector{Sample{T}}(undef, n)
    for i in 1:n
        if verbose println("Sample ", i) end
        s = sample(ϕ, θ, q, verbose)
        samples[i] = s
        q = s.value
        if verbose 
            print("Completed :: LL ", round(s.ll, digits = 4))
            print(" :: ", s.accepted ? "Accepted" : "Rejected")
            println(s.diverged ? " :: Diverged" : "") 
        end
    end
    return samples
end

