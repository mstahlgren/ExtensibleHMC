import StatsBase: Histogram, autocor, fit
import Plots: Plots, plot, scatter, quiver!

export effective_size, effective_size2, autocor, acceptance_rate
export path_plot, plot, scatter, density

# Here I am assuming autocor operates on each column in a matrix
StatsBase.autocor(samples::Vector{Sample}) = autocor(Matrix(samples))

# Canonical. Suggests broken level negative correlation is good
effective_size(x) = length(x) / (1 + 2 * sum(autocor(x)[2:end,:]))

# Pretty effective. Sometimes seems too small?
effective_size2(x) = length(x) / (1 + 2 * sum(abs.(autocor(x)[2:end,:])))

acceptance_rate(x) = count([s.accepted for s in x]) / length(x)

function Base.summary(samples::Vector{Sample})
    M = Matrix(samples)
    println("Number of free variables: ", size(M)[2])
    println("Number of samples: ", length(samples))
    println("Effective size: ", effective_size(M))
    println("Effective size 2: ", effective_size2(M))
    println("Acceptance rate: ", count([s.accepted for s in samples])/length(samples))
    println("Average trajectory length: ", sum([length(s.path) for s in samples])/length(samples))
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

# TODO: Put each variable in separate plots
function Plots.plot(samples::Vector{Sample}, idx = 1)
    Matrix(samples)[:,idx] |> plot
end

function Plots.scatter(samples::Vector{Sample}, id1 = 1, id2 = 2)
    M = Matrix(samples)
    scatter(M[:,id1], M[:,id2])
end

function density(samples::Vector{Sample}, idx::Int = 1)
    M = Matrix(samples)[:,idx]
    low = minimum(M)
    high = maximum(M)
    bins = low:(high-low)/32:high
    hist = fit(Histogram, M, bins)
    plot(hist, legend = false)
end