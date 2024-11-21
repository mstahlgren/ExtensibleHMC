import StatsBase: mean, quantile, Histogram, fit, autocor
import RecipesBase: @recipe

export samples, summary, acceptrate, miness

samples(x::Samples, f = identity) = [f(s.value) for s in x]

StatsBase.autocor(x::Samples) = autocor(reduce(hcat, samples(x, vec))')

miness(x::Samples) = length(x) / (2 * maximum(sum(autocor(x), dims = 1)) - 1)

acceptrate(x::Samples) = sum(s.acceptrate * s.nsteps for s in x) / sum(s.nsteps for s in x)

ndivergences(x::Samples) = sum(s.diverged for s in x)

# Samples per nsteps
function Base.summary(S::Samples)
    qs = (0, 0.1, 0.5, 0.9, 1)
    rnd = x -> round(x; digits = 3)
    println("Number of free variables: ", length(S[1].value))
    println("Number of samples: ", length(S))
    println("Effective sample size: ", rnd(miness(S)))
    println("Acceptance rate: ", rnd(acceptrate(S)))
    println("Divergence rate: ", rnd(ndivergences(S) / length(S)))
    println("Neff per 1000 steps: ", rnd(1000*miness(S) / sum(s.nsteps for s in S)))
    println("Number of steps: ", rnd.(quantile([s.nsteps for s in S], qs)))
    println("Autocorrelation τ¹⁻³: ", rnd.([mean(s) for s in eachrow(autocor(S)[2:4])]))
    println("Posterior ll: ", rnd.(quantile([s.ll for s in S], qs)))
    return nothing
end

# QT plot
@recipe function plot(S::Samples, idx)
    label := nothing
    return samples(S, x->x[idx])
end

# QQ scatter plot
@recipe function scatter(S::Samples, id₁, id₂)
    q = reduce(hcat, samples(S, x->x[[id₁, id₂]]))'
    label := nothing
    alpha := 100/length(S)
    view(q,:,id₁), view(q,:,id₂)
end
