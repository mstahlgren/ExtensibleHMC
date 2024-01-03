import Zygote: gradient
import LinearAlgebra: ⋅

export Euclidean, ∇, mass

abstract type Hamiltonian end

struct Euclidean{T <: Mass} <: Hamiltonian
    U::Function # Potential energy | Negative PDF
    mass::T # LLᵀ decomposition
end

∇(H::Hamiltonian) = x->gradient(H.U, x)[1]

(H::Hamiltonian)(s) = H.U(q(s)), 0.5 * p(s) ⋅ (H.mass\p(s))

mass(H::Hamiltonian) = H.mass