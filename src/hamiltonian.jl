import ComposableDistributions: IsoMvNormal, logpdf

export Euclidean, ∇

abstract type Hamiltonian end

struct Euclidean{T} <: Hamiltonian
    U::Function # Potential energy function (neg dens)
    m::T
end

(H::Hamiltonian)(q) = H.U(q)

(H::Euclidean)(q, p) = H(q) - logpdf(IsoMvNormal(length(p)), p)

∇(H::Hamiltonian) = x->gradient(H.U, x)[1]