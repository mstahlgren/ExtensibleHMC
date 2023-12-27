import StatsBase: mean, quantile, Histogram, autocor, fit
import Plots: Plots, plot, scatter, quiver!

export effective_size, effective_size2, autocor, acceptance_rate
export plot, scatter, density

StatsBase.autocor(samples::Vector{Sample}) = autocor(Matrix(samples))

#effective_size(x) = length(x) / (1 + 2 * sum(autocor(x)[2:end,:])/length(x[1].q))

#effective_size2(x) = length(x) / (1 + 2 * sum(abs.(autocor(x)[2:end,:]))/length(x[1].q))

acceptance_rate(x) = mean([s.accepted for s in x])

function Base.summary(samples::Vector{Sample{T}}) where T
    println("Number of free variables: ", length(samples[1].state))
    println("Number of samples: ", length(samples))
    #println("Effective size: ", effective_size(samples))
    #println("Effective size 2: ", effective_size2(samples))
    println("Acceptance rate: ", acceptance_rate(samples))
    println("Trajectory length: ", quantile([length(s.path) for s in samples], (0, 0.1, 0.5, 0.9, 1)))
    println("Posterior ll: ", quantile([s.ll for s in samples], ((0, 0.1, 0.5, 0.9, 1))))
end

function Plots.plot(path::Vector{State}, id1 = 1, id2 = 2)
    Q = reduce(hcat, [s.q for s in path])'
    P = reduce(hcat, [s.p for s in path])'

    QS, PS = (Q'.+Q[1,:])', (P'.+P[1,:])'
    QD, PD = (Q'.-Q[1,:])', (P'.-P[1,:])'

    scatter(Q[:,id1], Q[:,id2], legend = false)
    quiver!(Q[:,id1], Q[:,id2], quiver = (P[:,id1], P[:,id2]))
    quiver!(Q[:,id1], Q[:,id2], quiver = (QD[:,id1], QD[:,id2]))
end

function Plots.plot(samples::Vector{Sample}, idx = 1)
    [s.state[idx...] for s in samples] |> plot
end

function Plots.scatter(samples::Vector{Sample}, id1 = 1, id2 = 2)
    M = Matrix(samples)
    scatter(M[:,id1], M[:,id2])
end

function density(samples::Vector{Sample}, idx = 1)
    S = [s.state[idx...] for s in samples]
    low = minimum(S)
    high = maximum(S)
    bins = low:(high-low)/32:high
    hist = fit(Histogram, S, bins)
    plot(hist, legend = false)
end