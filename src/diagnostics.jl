Base.values(x::Samples) = reduce(hcat, vec(s.value) for s in x)

Base.values(x::Samples, idx...) = reduce(hcat, s.value[[idx...]] for s in x)

# This implementation differs from standard ones, but gives much better (?) results for negative autocorrelations.
function ess(x)
    a = autocor(x)    
    l = length(a) + iseven(length(a)) - 1
    M = a[1:2:l] .+ a[2:2:l]
    i = findfirst(x->x<0, M)
    if isnothing(i) i = length(M)+1 end
    length(x) / sum(M[1:i-1])
end

ess(x::Samples) = [z |> ess for z in x |> values |> transpose |> eachcol]

acceptrate(x::Samples) = sum(s.acceptrate * s.nsteps for s in x) / sum(s.nsteps for s in x)

ndivergences(x::Samples) = sum(s.diverged for s in x)

function Base.summary(S::Samples)
    rnd(x) = round(x; digits = 3)
    qs, N, E = (0, 0.1, 0.5, 0.9, 1), length(S), ess(S)
    Emin = argmin(E)
    println("Number of free variables: ", length(S[1].value))
    println("Number of samples: ", N)
    println("Minimum ESS: ", E[Emin] |> rnd, "($Emin)")
    println("Average ESS: ", mean(E) |> rnd)
    println("Acceptance rate: ", acceptrate(S) |> rnd)
    println("Divergences: ", ndivergences(S) |> rnd)
    println("N minESS / steps: ", E[Emin] / sum([s.nsteps for s in S]))
    println("Number of steps: ", quantile([s.nsteps for s in S], qs) .|> rnd)
    println("Posterior ll: ", quantile([s.ll for s in S], qs) .|> rnd)
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
