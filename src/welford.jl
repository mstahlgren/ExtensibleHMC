struct Welford{T, S <: AbstractVecOrMat{T}}
    n::Int
    mean::S
    M2::S
end

Welford(n) = Welford{Float64, Vector{Float64}}(0, zeros(n), zeros(n))

Welford(m, n) = Welford{Float64, Matrix{Float64}}(0, zeros(m, n), zeros(m, n))

Welford(n, m, m2) = Welford{Float64, typeof(m)}(n, m, m2)

Statistics.mean(w::Welford) = w.mean

Statistics.var(w::Welford) = w.M2 / (w.n - 1)

function update!(w::Welford, samples::Samples{T}) where T <: AbstractVecOrMat
    for s in samples w = update!(w, s.value) end
    return w
end

function update!(w::Welford, value)
    n = w.n + 1
    delta = value - w.mean
    w.mean .= w.mean .+ delta ./ n
    w.M2 .+= delta .* (value .- w.mean)
    return Welford(n, w.mean, w.M2)
end

function colsum(w::Welford)
    nb, na, m, m2 = w.n, w.n, w.mean[:, 1], w.M2[:, 1]
    for i in 2:size(w.mean)[2]
        δ = view(w.mean, :, i) - m
        m .+= δ .* nb ./ (na + nb)
        m2 .+= view(w.M2, :, i) .+ δ.^2 .* na * nb / (na + nb)
        na += nb
    end
    Welford{eltype(m), typeof(m)}(na, m, m2)
end