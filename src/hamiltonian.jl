import LinearAlgebra: Cholesky, cholesky

export Euclidean, ∇

abstract type Hamiltonian end

struct Euclidean{T} <: Hamiltonian
    U::Function # Potential energy | Negative PDF
    mass::Cholesky{T}
end

Euclidean(U, x::AbstractMatrix) = Euclidean(U, cholesky(x))

∇(H::Hamiltonian) = x->gradient(H.U, x)[1]

# CONSIDER: Should this be logpdf or the energy itself? pᵀM⁻¹p/2
# (H::Euclidean)(q, p) = H(q) - logpdf(IsoMvNormal(length(p)), p)
# (H::Euclidean)(q, p) = H(q) - p'*(H.mass\p)
(H::Euclidean)(s) = H.U(q(s)), 0.5*p(s)'*(H.mass\p(s))
