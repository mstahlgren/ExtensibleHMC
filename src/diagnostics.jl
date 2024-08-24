import StatsBase: mean, quantile, Histogram, fit, autocor
import Plots: Plots, plot, scatter

export summary, acceptrate, effectivesize
export plot, scatter, density

const Samples{T} = Vector{Sample{T}}

StatsBase.autocor(x::Samples) = autocor(reduce(hcat, [s.value for s in x])')

effectivesize(x::Samples) = mean(length(x) ./ (2 .* [sum(a) for a in eachcol(autocor(x))] .- 1))

acceptrate(x) = sum(s.acceptrate * s.nsteps for s in x) / sum(s.nsteps for s in x)

ndivergences(x) = sum(s.diverged for s in x)

function Base.summary(S::Vector{Sample{T}}) where T
    qs = (0, 0.1, 0.5, 0.9, 1)
    rnd = x -> round(x; digits = 3)
    println("Number of free variables: ", length(S[1].value))
    println("Number of samples: ", length(S))
    println("Autocorrelation τ¹⁻³: ", rnd.([mean(s) for s in eachrow(autocor(S)[2:4])]))
    println("Acceptance rate: ", rnd(acceptrate(S)))
    println("Divergance rate: ", rnd(ndivergences(S) / length(S)))
    println("Number of steps: ", rnd.(quantile([s.nsteps for s in S], qs)))
    println("Posterior ll: ", rnd.(quantile([s.ll for s in S], qs)))
    return nothing
end

# QT plot
function Plots.plot(samples::Vector{Sample{T}}, idx...) where T
    [s.value[idx...] for s in samples] |> plot
end

# QQ scatter plot
function Plots.scatter(samples::Vector{Sample{T}}, id₁, id₂) where T
    q = reduce(hcat, [s.value[[id₁, id₂]] for s in samples])
    scatter(view(q,:,1), view(q,:,2), xlim = (-3,3), ylim = (-3,3), aspect_ratio = :equal)
end

function density(samples::Vector{Sample{T}}, idx) where T
    S = [s.value[idx] for s in samples]
    low = minimum(S)
    high = maximum(S)
    bins = low:(high-low)/32:high
    hist = fit(Histogram, S, bins)
    plot(hist, legend = false)
end
