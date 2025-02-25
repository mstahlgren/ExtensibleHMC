samples(x::Samples, f = identity) = [f(s.value) for s in x]

#StatsBase.autocor(x::Samples) = autocor(reduce(hcat, samples(x, vec))')
autocor2(x::Samples) = autocor(reduce(hcat, samples(x, vec))')

ess(x::Samples, agg = maximum) = length(x) / (2 * agg(sum(autocor2(x), dims = 1)) - 1)

acceptrate(x::Samples) = sum(s.acceptrate * s.nsteps for s in x) / sum(s.nsteps for s in x)

ndivergences(x::Samples) = sum(s.diverged for s in x)

# Samples per nsteps
function Base.summary(S::Samples)
    qs = (0, 0.1, 0.5, 0.9, 1)
    rnd = x -> round(x; digits = 3)
    println("Number of free variables: ", length(S[1].value))
    println("Number of samples: ", length(S))
    println("Minimum ESS: ", rnd(ess(S)))
    println("Average ESS: ", rnd(ess(S, mean)))
    println("Highest autocorr: ", argmax(vec(sum(autocor2(S), dims = 1))))
    println("Acceptance rate: ", rnd(acceptrate(S)))
    println("Divergence rate: ", rnd(ndivergences(S) / length(S)))
    println("N minESS per 1000 steps: ", rnd(1000*ess(S) / sum(s.nsteps for s in S)))
    println("Number of steps: ", rnd.(quantile([s.nsteps for s in S], qs)))
    println("Autocorrelation τ¹⁻³: ", rnd.([mean(s) for s in eachrow(autocor2(S)[2:4])]))
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
