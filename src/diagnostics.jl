Base.values(x::Samples) = reduce(hcat, vec(s.value) for s in x)

Base.values(x::Samples, idx...) = reduce(hcat, s.value[[idx...]] for s in x)

StatsBase.autocor(x::Samples) = autocor(Base.values(x)')

ess(x::Matrix, N, agg) = agg(vec(N ./ (2 .* sum(abs.(x), dims = 1) .- 1)))

ess(x::Samples, agg) = ess(autocor(x), length(x), agg) 

acceptrate(x::Samples) = sum(s.acceptrate * s.nsteps for s in x) / sum(s.nsteps for s in x)

ndivergences(x::Samples) = sum(s.diverged for s in x)

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
    println("N minESS per 100 steps: ", rnd(100*ess(ac, N, minimum) / sum(s.nsteps for s in S)))
    println("Number of steps: ", rnd.(quantile([s.nsteps for s in S], qs)))
    println("Autocorrelation τ¹⁻³: ", rnd.(vec(mean(ac[2:4, :], dims = 2))))
    print("Posterior ll: ", rnd.(quantile([s.ll for s in S], qs)))
end

# QT plot
@recipe function f(S::Samples, id)
    seriestype := :line
    label := nothing
    return values(S, id)'
end

# QQ scatter plot
@recipe function f(S::Samples, id₁, id₂)
    seriestype := :scatter
    label := nothing
    alpha := 100 / length(S)
    q = values(S, id₁, id₂)'
    view(q, :, id₁), view(q, :, id₂)
end
