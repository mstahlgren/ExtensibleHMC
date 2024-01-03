import StatsBase: StatsBase, sample

export Sample, sample

abstract type Sampler end

stepsize(x::Sampler) = x.ϵ

struct Sample{T <: AbstractVecOrMat}
    value::T
    ll::Float64
    accepted::Bool
    path::Vector{State{T}}
end

# CONSIDER: Time of samples?
# CONSIDER: Progress bar?
# CONSIDER: How to report divergences?
function StatsBase.sample(ϕ::Sampler, θ::Hamiltonian, q::T, n) where T <: AbstractVecOrMat
    samples = Vector{Sample{T}}(undef, n)
    for i in 1:n
        s = sample(ϕ, θ, q)
        samples[i] = s
        q = s.value
    end
    return samples
end

