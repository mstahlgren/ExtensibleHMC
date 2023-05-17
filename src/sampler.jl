import StatsBase: StatsBase, Weights, sample

export Debug, Metropolis, Proportional

abstract type Sampler end

struct Debug <: Sampler end

struct Metropolis <: Sampler end

struct Proportional <: Sampler end

function (::Debug)(path, θ) path end

function (::Metropolis)(path, θ)
    ΔH = θ(first(path)...) - θ(last(path)...)
    return ΔH > 0.0 || exp(ΔH) > rand() ? last(path)[1] : first(path)[1]
end

function (::Proportional)(path, θ)
    P = exp.([θ(qp...) for qp in path])
    return path[sample(Weights(P))][1]
end