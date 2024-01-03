import Zygote: gradient
import LinearAlgebra: ⋅

export Hamiltonian, ∇, mass

struct Hamiltonian{T <: Mass}
    U::Function # Potential energy | Negative PDF
    mass::T # LLᵀ decomposition
end

∇(H::Hamiltonian) = x->gradient(H.U, x)[1]

(H::Hamiltonian)(s) = H.U(q(s)), 0.5 * p(s) ⋅ (H.mass\p(s))

mass(H::Hamiltonian) = H.mass