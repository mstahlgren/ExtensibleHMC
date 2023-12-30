import StatsBase: mean, quantile, Histogram, autocor, fit
import Plots: Plots, plot, scatter, quiver!, scatter!

export effective_size, effective_size2, autocor, acceptance_rate
export plot, scatter, path, density

extract(samples::Vector{Sample{T}}, idx...) where T = [s.state[idx...] for s in samples]

StatsBase.autocor(samples::Vector{Sample}) = autocor(Matrix(samples))

#effective_size(x) = length(x) / (1 + 2 * sum(autocor(x)[2:end,:])/length(x[1].q))

#effective_size2(x) = length(x) / (1 + 2 * sum(abs.(autocor(x)[2:end,:]))/length(x[1].q))

acceptance_rate(x) = mean([s.accepted for s in x])

function Base.summary(samples::Vector{Sample{T}}) where T
    println("Number of free variables: ", length(samples[1].state))
    println("Number of samples: ", length(samples))
    #println("Effective size: ", effective_size(samples))
    #println("Effective size 2: ", effective_size2(samples))
    println("Acceptance rate: ", round(acceptance_rate(samples); digits = 2))
    println("Trajectory length: ", round.(quantile([length(s.path) for s in samples], (0, 0.1, 0.5, 0.9, 1)); digits = 3))
    println("Posterior ll: ", round.(quantile([s.ll for s in samples], ((0, 0.1, 0.5, 0.9, 1))); digits = 3))
end

# QT plot
function Plots.plot(samples::Vector{Sample{T}}, idx...) where T
    extract(samples, idx...) |> plot
end

# QQ scatter plot
function Plots.scatter(samples::Vector{Sample{T}}, id₁, id₂) where T
    q₁, q₂ = extract(samples, id₁), extract(samples, id₂)
    scatter(q₁, q₂)
end

# OLDs
function Plots.scatter(path::Vector{State}, id₁, id₂)
    Q = reduce(hcat, [s.q[[id₁, id₂]] for s in path])'
    P = reduce(hcat, [s.p for s in path])'

    QS, PS = (Q'.+Q[1,:])', (P'.+P[1,:])'
    QD, PD = (Q'.-Q[1,:])', (P'.-P[1,:])'

    scatter(Q[:,id1], Q[:,id2], legend = false)
    quiver!(Q[:,id1], Q[:,id2], quiver = (P[:,id1], P[:,id2]))
    quiver!(Q[:,id1], Q[:,id2], quiver = (QD[:,id1], QD[:,id2]))
end

function density(samples::Vector{Sample}, idx...)
    S = [s.state[idx...] for s in samples]
    low = minimum(S)
    high = maximum(S)
    bins = low:(high-low)/32:high
    hist = fit(Histogram, S, bins)
    plot(hist, legend = false)
end

# QP path plot
function path(sample::Sample{T}, id) where T
    q, p = [s.q[id] for s in sample.path], [s.p[id] for s in sample.path]
    #δq, δp, Σp = q .- q[1], p .- p[1], p[1]
    scatter(q, p, legend = false, xlabel = "Position", ylabel = "Momentum")
   #scatter!([q[1]], [p[1]])
end

# Q₁Q₂ path plot
function path(sample::Sample{T}, id₁, id₂) where T
    d = reduce(hcat, [[s.q[id₁], s.q[id₂], s.p[id₁], s.p[id₂]] for s in sample.path])'
    δq, δp, Σp = d[:,1] .- d[:,2], d[:,3] .- d[:,4], d[:,3] + d[:,4]
    scatter(d[:,1], d[:,2], legend = false, xlabel = "Position", ylabel = "Position")
    scatter!([d[1,1]], [d[1,2]])
    quiver!(d[:,1], d[:,2], quiver = (d[:,3], d[:,4]))
end