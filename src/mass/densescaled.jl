struct DenseScaled{S, M} <: Mass
    McL::S
    M⁻¹::M
    var::M
    N::Int
end

function DenseScaled(N, P)
    V, C = zeros(N, P), zeros(N, N)
end

function DenseScaled(S::Samples{T}, m = DenseScaled(size(S[1].value, 1))) where T <: Matrix
    ns, P = length(S), size(S[1].value, 2)
    V, C = fill(undef, ns, P), zeros(ns, ns)
    M = reduce(vcat, vec(s.value)' for s in S)
    for i in 1:P
        m = view(M, :, i*N-N+1:i*N)
        C .+= cor(m, corrected = false)
        V[:, i] .= std(m, corrected = false)
    end
    N += P
    c = cholesky!(Hermitian(cholesky(cov) |> inv)).L

    DenseMass{eltype(cov), typeof(c), typeof(cov)}(c, cov, sqrt.())
end

function DenseScaled(m::DenseScaled, s::Samples{T}) where T

end

Base.:*(m::Mass, x) = m.McL * x

Base.:\(m::Mass, x) =  m.M⁻¹ * x