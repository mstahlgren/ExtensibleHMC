import StatsBase: mean, quantile, Histogram, fit, autocor
import Plots: Plots, plot, scatter

export summary, acceptrate
export plot, scatter, density

const Samples{T} = Vector{Sample{T}}

extract(samples::Vector{Sample{T}}, idx...) where T = [s.value[idx...] for s in samples]

#StatsBase.autocor(samples::Vector{Sample}) = autocor(Matrix(samples))

#effective_size(x) = length(x) / (1 + 2 * sum(autocor(x)[2:end,:])/length(x[1].q))

#effective_size2(x) = length(x) / (1 + 2 * sum(abs.(autocor(x)[2:end,:]))/length(x[1].q))

naccepted(x) = sum(s.accepted for s in x)

ndivergences(x) = sum(s.diverged for s in x)

function Base.summary(samples::Vector{Sample{T}}) where T
    n = length(samples)
    qs = (0, 0.1, 0.5, 0.9, 1)
    rnd = x -> round(x; digits = 3)
    println("Number of free variables: ", length(samples[1].value))
    println("Number of samples: ", n)
    println("Autocoorelation τ¹⁻³: ", rnd.(autocor(extract(samples, 1))[2:4]))
    println("Acceptance rate: ", rnd(naccepted(samples) / n))
    println("Divergance rate: ", rnd(ndivergences(samples) / n))
    println("Number of steps: ", rnd.(quantile([s.nsteps for s in samples], qs))),
    println("Posterior ll: ", rnd.(quantile([s.ll for s in samples], qs)))
    return nothing
end

# QT plot
function Plots.plot(samples::Vector{Sample{T}}, idx...) where T
    extract(samples, idx...) |> plot
end

# QQ scatter plot
function Plots.scatter(samples::Vector{Sample{T}}, id₁, id₂) where T
    q₁, q₂ = extract(samples, id₁), extract(samples, id₂)
    scatter(q₁, q₂, xlim = (-3,3), ylim = (-3,3), aspect_ratio = :equal)
end

function density(samples::Vector{Sample{T}}, idx...) where T
    S = [s.value[idx...] for s in samples]
    low = minimum(S)
    high = maximum(S)
    bins = low:(high-low)/32:high
    hist = fit(Histogram, S, bins)
    plot(hist, legend = false)
end
