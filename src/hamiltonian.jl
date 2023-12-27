import LinearAlgebra: Cholesky, cholesky

export Euclidean, ∇

abstract type Hamiltonian end

struct Euclidean{T} <: Hamiltonian
    U::Function # Potential energy | Negative PDF
    mass::T
end

∇(H::Hamiltonian) = x->gradient(H.U, x)[1]

(H::Euclidean)(s) = H.U(q(s)), 0.5*sum(p(s).*(H.mass\p(s)))
