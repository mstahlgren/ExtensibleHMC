Base.values(x::Samples) = reduce(hcat, s.value for s in x)

StatsBase.autocor(x::Samples) = autocor(Base.values(x)')

ess(x::Matrix, N, agg) = agg(vec(N ./ (2 .* sum(x, dims = 1) .- 1)))

ess(x::Samples, agg) = ess(autocor(x), length(x), agg)

acceptrate(x::Samples) = sum(s.acceptrate * s.nsteps for s in x) / sum(s.nsteps for s in x)

ndivergences(x::Samples) = sum(s.diverged for s in x)

# Samples per nsteps
function Base.summary(S::Samples)
    rnd(x) = round(x; digits = 3)
    qs = (0, 0.1, 0.5, 0.9, 1)
    N, ac = length(S), autocor(S)
    println("Number of free variables: ", length(S[1].value))
    println("Number of samples: ", N)
    println("Minimum ESS: ", rnd(ess(ac, N, minimum)))
    println("Average ESS: ", rnd(ess(ac, N, mean)))
    println("Highest autocorr: ", ess(ac, N, argmax))
    println("Acceptance rate: ", rnd(acceptrate(S)))
    println("Divergences: ", rnd(ndivergences(S)))
    println("N minESS per 100 steps: ", rnd(100*ess(S, minimum) / sum(s.nsteps for s in S)))
    println("Number of steps: ", rnd.(quantile([s.nsteps for s in S], qs)))
    println("Autocorrelation τ¹⁻³: ", rnd.(vec(mean(ac[2:4, :], dims = 1))))
    print("Posterior ll: ", rnd.(quantile([s.ll for s in S], qs)))
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
