struct DenseMass{S, M} <: Mass
    McL::S
    M⁻¹::M
end

function DenseMass(cov)
    c = cholesky(Hermitian(cholesky(cov) |> inv)).L
    DenseMass{eltype(cov), typeof(c), typeof(cov)}(c, cov)
end
