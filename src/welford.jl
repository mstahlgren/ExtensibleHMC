import Statistics: Statistics

export Welford, update, var

struct Welford{T,S <: AbstractVecOrMat{T}}
    n::Int
    mean::S
    M2::S
end

Welford(v) = Welford{eltype(v), typeof(v)}(0, zeros(size(v)), zeros(size(v)))

Statistics.mean(w::Welford) = w.mean

Statistics.var(w::Welford) = w.M2 / (w.n - 1)

function update(w::Welford, values::Vector{T}) where T <: AbstractVecOrMat
    for v in values w = update(w, v) end
    return w
end

function update(w::Welford, value)
    n = w.n + 1
    delta = value - w.mean
    w.mean .= w.mean .+ delta ./ n
    w.M2 .+= delta .* (value .- mean)
    return w
end
