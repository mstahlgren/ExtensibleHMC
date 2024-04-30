import LinearAlgebra: ⋅
import Zygote: gradient

export Hamiltonian, ∇, mass, gradient

struct Hamiltonian{T <: Mass}
    U::Function # Unnormalized log posterior
    mass::T # LLᵀ decomposition
end

gradient(H::Hamiltonian, s::AbstractArray) = gradient(H.U, s)[1]

energy(H::Hamiltonian, s) = kinetic(H, s) + potential(H, s)

kinetic(H::Hamiltonian, s) = 0.5 * p(s)⋅(H.mass\p(s))

potential(H::Hamiltonian, s) = -H.U(q(s))

mass(H::Hamiltonian) = H.mass